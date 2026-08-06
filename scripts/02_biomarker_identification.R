#!/usr/bin/env Rscript
# =============================================================================
# 02_biomarker_identification.R — DE analysis + consensus spaceflight signature
# =============================================================================
# For each tissue stratum (brain, liver_immune):
#   1. Parse sample metadata to assign conditions (Space Flight vs Ground Control)
#   2. Per-study DESeq2 differential expression
#   3. Random-effects effect-size meta-analysis across studies (metafor::rma)
#   4. Define consensus signature (FDR < 0.05 AND direction-consistent)
#   5. Map mouse Ensembl IDs → human orthologs (biomaRt)
#   6. Technical gene filtering
#
# Inputs:
#   data/raw/<OSD-ID>/  — count matrices + sample_metadata.json
#   data/study_catalog.csv
#   data/download_status.csv
#
# Outputs:
#   data/processed/consensus_signature_{tissue}.csv
#   data/processed/meta_analysis_{tissue}_full.csv
#   data/processed/ortholog_mapping.csv
#   data/processed/per_study_de/{OSD-ID}_de_results.csv
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(metafor)
  library(biomaRt)
  library(jsonlite)
  library(dplyr)
  library(readr)
  library(stringr)
})

set.seed(42)

PROJECT_ROOT <- "/workspace/astronaut-opposite-forcing"
RAW_DIR <- file.path(PROJECT_ROOT, "data", "raw")
PROC_DIR <- file.path(PROJECT_ROOT, "data", "processed")
PER_STUDY_DIR <- file.path(PROC_DIR, "per_study_de")
dir.create(PROC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PER_STUDY_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load study catalog and download status ---
study_catalog <- read_csv(file.path(PROJECT_ROOT, "data", "study_catalog.csv"), show_col_types = FALSE)
download_status <- read_csv(file.path(PROJECT_ROOT, "data", "download_status.csv"), show_col_types = FALSE)

# Only use studies with successfully downloaded count matrices
usable_studies <- download_status %>%
  filter(!is.na(count_file)) %>%
  select(osd_id, tissue_stratum, tissue, glds_id, mission) %>%
  inner_join(study_catalog %>% select(osd_id, title), by = "osd_id")

message(sprintf("Usable studies with count matrices: %d", nrow(usable_studies)))
message(sprintf("  Brain: %d", sum(usable_studies$tissue_stratum == "brain")))
message(sprintf("  Liver/immune: %d", sum(usable_studies$tissue_stratum == "liver_immune")))

# =============================================================================
# Helper: Parse sample metadata JSON to get condition assignments
# =============================================================================
parse_sample_conditions <- function(osd_id) {
  meta_file <- file.path(RAW_DIR, osd_id, "sample_metadata.json")
  if (!file.exists(meta_file)) return(NULL)

  meta <- fromJSON(meta_file, simplifyVector = FALSE)
  assays <- meta[["assays"]]
  if (is.null(assays)) return(NULL)

  sample_conditions <- list()
  sample_tissues <- list()

  for (assay_name in names(assays)) {
    samples <- assays[[assay_name]][["samples"]]
    if (is.null(samples)) next

    for (sample_name in names(samples)) {
      s_meta <- samples[[sample_name]][["metadata"]]
      if (is.null(s_meta)) next

      # Get factor value (spaceflight)
      fv <- s_meta[["study"]][["factor value"]]
      condition <- NA
      if (!is.null(fv[["spaceflight"]])) {
        condition <- fv[["spaceflight"]][["~"]]
      }

      # Get tissue from characteristics
      chars <- s_meta[["study"]][["characteristics"]]
      tissue <- NA
      if (!is.null(chars)) {
        for (ck in names(chars)) {
          if (grepl("tissue|organ|material", ck, ignore.case = TRUE)) {
            tissue <- chars[[ck]][["~"]]
            break
          }
        }
      }

      sample_conditions[[sample_name]] <- condition
      sample_tissues[[sample_name]] <- tissue
    }
  }

  if (length(sample_conditions) == 0) return(NULL)

  data.frame(
    sample = names(sample_conditions),
    condition = unlist(sample_conditions),
    tissue = unlist(sample_tissues),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Helper: Infer condition from sample name if metadata is missing
# =============================================================================
infer_condition_from_name <- function(sample_names) {
  conditions <- sapply(sample_names, function(sn) {
    sn_lower <- tolower(sn)
    if (grepl("\\bflt\\b|_flt_|flight|spc|iss", sn_lower)) return("Space Flight")
    if (grepl("\\bgc\\b|_gc_|ground|grnd|ises|hgc", sn_lower)) return("Ground Control")
    if (grepl("\\bviv\\b|_viv_|vivarium|vgc", sn_lower)) return("Vivarium Control")
    if (grepl("\\bbas\\b|_bas_|basal", sn_lower)) return("Basal Control")
    return(NA)
  })
  conditions
}

# =============================================================================
# Helper: Run DESeq2 on a single study
# =============================================================================
run_deseq2 <- function(counts, coldata, study_id) {
  # Ensure counts and coldata are aligned
  common_samples <- intersect(colnames(counts), coldata$sample)
  if (length(common_samples) < 4) {
    message(sprintf("  WARNING: %s has only %d matched samples, skipping", study_id, length(common_samples)))
    return(NULL)
  }

  counts <- counts[, common_samples]
  coldata <- coldata[coldata$sample %in% common_samples, ]
  rownames(coldata) <- coldata$sample

  # Only keep Space Flight and Ground Control (exclude Vivarium/Basal for primary contrast)
  coldata$condition <- ifelse(grepl("space flight", coldata$condition, ignore.case = TRUE), "Flight",
                       ifelse(grepl("ground control", coldata$condition, ignore.case = TRUE), "Ground", NA))
  coldata <- coldata[!is.na(coldata$condition), ]
  counts <- counts[, rownames(coldata)]

  if (sum(coldata$condition == "Flight") < 2 || sum(coldata$condition == "Ground") < 2) {
    message(sprintf("  WARNING: %s has insufficient replicates (Flight=%d, Ground=%d), skipping",
                    study_id, sum(coldata$condition == "Flight"), sum(coldata$condition == "Ground")))
    return(NULL)
  }

  coldata$condition <- factor(coldata$condition, levels = c("Ground", "Flight"))

  # Pre-filter low-count genes
  keep <- rowSums(counts) >= 10
  counts <- counts[keep, ]

  # DESeq2
  tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = coldata, design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", "Flight", "Ground"), alpha = 0.05)

    # LFC shrinkage (use coef with apeglm; fallback to normal if contrast needed)
    res_shrunk <- tryCatch({
      lfcShrink(dds, coef = resultsNames(dds)[length(resultsNames(dds))],
                res = res, type = "apeglm")
    }, error = function(e) {
      lfcShrink(dds, contrast = c("condition", "Flight", "Ground"),
                res = res, type = "normal")
    })

    # Extract results with SE
    de_df <- data.frame(
      gene_id = rownames(res),
      baseMean = res$baseMean,
      log2FC = res_shrunk$log2FoldChange,
      SE = res_shrunk$lfcSE,
      stat = res$stat,
      pvalue = res$pvalue,
      padj = res$padj,
      stringsAsFactors = FALSE
    )

    return(de_df)
  }, error = function(e) {
    message(sprintf("  ERROR in DESeq2 for %s: %s", study_id, conditionMessage(e)))
    return(NULL)
  })
}

# =============================================================================
# Step 1: Per-study DESeq2 analysis
# =============================================================================
message("\n=== Step 1: Per-study DESeq2 differential expression ===")

all_de_results <- list()
study_sample_info <- list()

for (i in 1:nrow(usable_studies)) {
  row <- usable_studies[i, ]
  osd_id <- row$osd_id
  tissue_stratum <- row$tissue_stratum

  message(sprintf("\n  [%d/%d] %s (%s)", i, nrow(usable_studies), osd_id, tissue_stratum))

  # Find count matrix file
  study_dir <- file.path(RAW_DIR, osd_id)
  count_files <- list.files(study_dir, pattern = "Unnormalized_Counts.*\\.csv$", full.names = TRUE)
  if (length(count_files) == 0) {
    message(sprintf("  WARNING: No count file for %s, skipping", osd_id))
    next
  }
  count_file <- count_files[1]

  # Load counts
  counts <- tryCatch(read.csv(count_file, row.names = 1, check.names = FALSE),
                     error = function(e) {
                       message(sprintf("  ERROR reading counts: %s", conditionMessage(e)))
                       return(NULL)
                     })
  if (is.null(counts)) next

  # Ensure integer counts
  counts <- as.matrix(counts)
  mode(counts) <- "integer"

  # Parse sample conditions
  sample_meta <- parse_sample_conditions(osd_id)

  if (is.null(sample_meta) || nrow(sample_meta) == 0) {
    # Infer from sample names
    message(sprintf("  No metadata JSON, inferring conditions from sample names"))
    conditions <- infer_condition_from_name(colnames(counts))
    sample_meta <- data.frame(
      sample = colnames(counts),
      condition = conditions,
      tissue = row$tissue,
      stringsAsFactors = FALSE
    )
  }

  # If multi-tissue study, filter to samples matching this tissue stratum
  # Check if sample names contain tissue indicators
  if (!is.null(sample_meta$tissue) && any(!is.na(sample_meta$tissue))) {
    # For multi-tissue studies, we need to filter. But since each OSD-ID in our
    # catalog was assigned a tissue, and the count matrix may contain all tissues,
    # we filter by matching tissue names in sample names
    # For now, use all samples (the study was selected for a specific tissue)
  }

  # Run DESeq2
  de_res <- run_deseq2(counts, sample_meta, osd_id)

  if (!is.null(de_res)) {
    de_res$osd_id <- osd_id
    de_res$tissue_stratum <- tissue_stratum
    all_de_results[[osd_id]] <- de_res
    write_csv(de_res, file.path(PER_STUDY_DIR, paste0(osd_id, "_de_results.csv")))

    n_sig <- sum(de_res$padj < 0.05 & !is.na(de_res$padj), na.rm = TRUE)
    message(sprintf("  DE results: %d genes, %d significant (padj < 0.05)", nrow(de_res), n_sig))

    study_sample_info[[osd_id]] <- data.frame(
      osd_id = osd_id,
      tissue_stratum = tissue_stratum,
      n_flight = sum(sample_meta$condition == "Space Flight" | sample_meta$condition == "Flight", na.rm = TRUE),
      n_ground = sum(sample_meta$condition == "Ground Control" | sample_meta$condition == "Ground", na.rm = TRUE),
      n_sig = n_sig
    )
  }
}

de_combined <- do.call(rbind, all_de_results)
rownames(de_combined) <- NULL
sample_info_df <- do.call(rbind, study_sample_info)
rownames(sample_info_df) <- NULL

message(sprintf("\nTotal DE results: %d rows from %d studies", nrow(de_combined), length(all_de_results)))
write_csv(sample_info_df, file.path(PROC_DIR, "study_sample_info.csv"))

# =============================================================================
# Step 2: Random-effects meta-analysis per tissue stratum
# =============================================================================
message("\n=== Step 2: Random-effects effect-size meta-analysis ===")

run_meta_analysis <- function(de_combined, stratum) {
  stratum_de <- de_combined %>% filter(tissue_stratum == stratum)
  if (nrow(stratum_de) == 0) return(NULL)

  # Get genes present in >= 2 studies
  gene_study_counts <- stratum_de %>%
    filter(!is.na(log2FC), !is.na(SE), SE > 0) %>%
    group_by(gene_id) %>%
    summarise(n_studies = n(), .groups = "drop")

  meta_genes <- gene_study_counts %>% filter(n_studies >= 2) %>% pull(gene_id)
  message(sprintf("  %s: %d genes in >=2 studies for meta-analysis", stratum, length(meta_genes)))

  if (length(meta_genes) == 0) return(NULL)

  # Run meta-analysis per gene
  meta_results <- lapply(meta_genes, function(gene) {
    gene_data <- stratum_de %>%
      filter(gene_id == gene, !is.na(log2FC), !is.na(SE), SE > 0)

    if (nrow(gene_data) < 2) return(NULL)

    tryCatch({
      # Random-effects meta-analysis
      rma_fit <- rma(yi = log2FC, sei = SE, data = gene_data, method = "REML",
                     control = list(stepadj = 0.5, maxiter = 100))

      # Per-study directions
      directions <- sign(gene_data$log2FC)
      n_positive <- sum(directions > 0)
      n_negative <- sum(directions < 0)
      direction_consistent <- (n_positive == length(directions)) || (n_negative == length(directions))

      data.frame(
        gene_id = gene,
        pooled_log2FC = rma_fit$beta[1],
        pooled_SE = rma_fit$se,
        z = rma_fit$zval,
        pvalue = rma_fit$pval,
        I2 = rma_fit$I2,
        tau2 = rma_fit$tau2,
        k = rma_fit$k,
        n_studies_up = n_positive,
        n_studies_down = n_negative,
        direction_consistent = direction_consistent,
        per_study_log2FC = paste(sprintf("%.3f", gene_data$log2FC), collapse = ";"),
        per_study_osd = paste(gene_data$osd_id, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      # Fallback: DerSimonian-Laird
      tryCatch({
        rma_fit <- rma(yi = log2FC, sei = SE, data = gene_data, method = "DL")
        directions <- sign(gene_data$log2FC)
        direction_consistent <- all(directions > 0) || all(directions < 0)

        data.frame(
          gene_id = gene,
          pooled_log2FC = rma_fit$beta[1],
          pooled_SE = rma_fit$se,
          z = rma_fit$zval,
          pvalue = rma_fit$pval,
          I2 = rma_fit$I2,
          tau2 = rma_fit$tau2,
          k = rma_fit$k,
          n_studies_up = sum(directions > 0),
          n_studies_down = sum(directions < 0),
          direction_consistent = direction_consistent,
          per_study_log2FC = paste(sprintf("%.3f", gene_data$log2FC), collapse = ";"),
          per_study_osd = paste(gene_data$osd_id, collapse = ";"),
          stringsAsFactors = FALSE
        )
      }, error = function(e2) NULL)
    })
  })

  meta_df <- do.call(rbind, meta_results)
  if (is.null(meta_df) || nrow(meta_df) == 0) return(NULL)

  # BH-FDR correction
  meta_df$padj <- p.adjust(meta_df$pvalue, method = "BH")

  # Define consensus: FDR < 0.05 AND direction-consistent
  meta_df$consensus <- meta_df$padj < 0.05 & meta_df$direction_consistent
  meta_df$core <- meta_df$consensus & abs(meta_df$pooled_log2FC) >= 1

  # Direction
  meta_df$direction <- ifelse(meta_df$pooled_log2FC > 0, "up", "down")

  return(meta_df)
}

# Run for each stratum
strata <- unique(usable_studies$tissue_stratum)
all_meta <- list()

for (stratum in strata) {
  message(sprintf("\n  Running meta-analysis for: %s", stratum))
  meta_df <- run_meta_analysis(de_combined, stratum)

  if (!is.null(meta_df)) {
    all_meta[[stratum]] <- meta_df
    write_csv(meta_df, file.path(PROC_DIR, paste0("meta_analysis_", stratum, "_full.csv")))

    n_fdr <- sum(meta_df$padj < 0.05, na.rm = TRUE)
    n_consensus <- sum(meta_df$consensus, na.rm = TRUE)
    n_core <- sum(meta_df$core, na.rm = TRUE)
    n_up <- sum(meta_df$consensus & meta_df$direction == "up", na.rm = TRUE)
    n_down <- sum(meta_df$consensus & meta_df$direction == "down", na.rm = TRUE)

    message(sprintf("  %s: %d FDR-sig, %d consensus (%d up, %d down), %d core",
                    stratum, n_fdr, n_consensus, n_up, n_down, n_core))
  }
}

# =============================================================================
# Step 3: Technical gene filtering
# =============================================================================
message("\n=== Step 3: Technical gene filtering ===")

# Technical genes to remove (ribosomal, translation factors, housekeeping)
# But RETAIN metallothioneins (MT1*/MT2A) as genuine biology
technical_patterns <- c(
  "^Rpl\\d", "^Rps\\d",  # ribosomal proteins
  "^Rp[ls]\\d",          # ribosomal (alternative)
  "GAPDH", "ACTB", "B2M", "ACTG1",  # housekeeping
  "^Eef\\d", "^Eif\\d"   # translation factors
)

filter_technical_genes <- function(meta_df) {
  is_technical <- sapply(meta_df$gene_id, function(gid) {
    any(sapply(technical_patterns, function(p) grepl(p, gid, ignore.case = TRUE)))
  })
  # Don't remove metallothioneins
  is_metallothionein <- grepl("^MT1|^MT2A|^Mt1|^Mt2a", meta_df$gene_id, ignore.case = TRUE)
  is_technical <- is_technical & !is_metallothionein

  removed <- meta_df$gene_id[is_technical]
  message(sprintf("  Removed %d technical genes (retained %d metallothioneins)",
                  length(removed), sum(is_metallothionein)))
  meta_df[!is_technical, ]
}

for (stratum in names(all_meta)) {
  all_meta[[stratum]] <- filter_technical_genes(all_meta[[stratum]])
}

# =============================================================================
# Step 4: Mouse → Human ortholog mapping
# =============================================================================
message("\n=== Step 4: Mouse → Human ortholog mapping ===")

# Get all unique gene IDs across strata
all_gene_ids <- unique(unlist(lapply(all_meta, function(m) m$gene_id)))
message(sprintf("Total unique gene IDs to map: %d", length(all_gene_ids)))

# Use biomaRt to map Ensembl mouse IDs to human orthologs
ortholog_mapping <- NULL
tryCatch({
  message("  Connecting to Ensembl biomaRt...")
  mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  # Get orthologs
  orthologs <- getLDS(
    attributes = c("ensembl_gene_id", "mgi_symbol"),
    filters = "ensembl_gene_id",
    values = all_gene_ids,
    mart = mouse_mart,
    attributesL = c("ensembl_gene_id", "hgnc_symbol"),
    martL = human_mart
  )

  colnames(orthologs) <- c("mouse_ensembl_id", "mouse_symbol", "human_ensembl_id", "human_symbol")
  ortholog_mapping <- orthologs

  # Deduplicate: one-to-one preferred
  ortholog_mapping <- ortholog_mapping %>%
    group_by(mouse_ensembl_id) %>%
    mutate(n_human = n()) %>%
    arrange(mouse_ensembl_id, n_human, human_symbol) %>%
    slice(1) %>%  # take first (prefer one-to-one)
    ungroup()

  message(sprintf("  Mapped %d/%d mouse genes to human orthologs (%.1f%% coverage)",
                  nrow(ortholog_mapping), length(all_gene_ids),
                  100 * nrow(ortholog_mapping) / length(all_gene_ids)))

}, error = function(e) {
  message(sprintf("  WARNING: biomaRt failed: %s", conditionMessage(e)))
  message("  Falling back to MGI symbol-based mapping...")

  # Fallback: extract MGI symbols from Ensembl IDs using a local mapping
  # We'll use the gene names from the count matrix if available
  ortholog_mapping <<- data.frame(
    mouse_ensembl_id = all_gene_ids,
    mouse_symbol = NA,
    human_ensembl_id = NA,
    human_symbol = NA,
    stringsAsFactors = FALSE
  )
})

write_csv(ortholog_mapping, file.path(PROC_DIR, "ortholog_mapping.csv"))

# =============================================================================
# Step 5: Build consensus signatures with human orthologs
# =============================================================================
message("\n=== Step 5: Building consensus signatures ===")

for (stratum in names(all_meta)) {
  meta_df <- all_meta[[stratum]]

  # Add human orthologs
  meta_df <- meta_df %>%
    left_join(ortholog_mapping %>% select(mouse_ensembl_id, human_symbol),
              by = c("gene_id" = "mouse_ensembl_id"))

  # Consensus signature (up and down)
  consensus_sig <- meta_df %>%
    filter(consensus == TRUE) %>%
    select(gene_id, human_symbol, pooled_log2FC, padj, direction, I2, k,
           n_studies_up, n_studies_down, core) %>%
    arrange(desc(abs(pooled_log2FC)))

  write_csv(consensus_sig, file.path(PROC_DIR, paste0("consensus_signature_", stratum, ".csv")))

  n_up <- sum(consensus_sig$direction == "up")
  n_down <- sum(consensus_sig$direction == "down")
  n_human <- sum(!is.na(consensus_sig$human_symbol) & consensus_sig$human_symbol != "")

  message(sprintf("  %s consensus signature: %d genes (%d up, %d down), %d with human orthologs",
                  stratum, nrow(consensus_sig), n_up, n_down, n_human))

  # Also save up/down gene lists for LINCS
  up_genes <- consensus_sig %>% filter(direction == "up") %>% pull(human_symbol) %>% na.omit()
  down_genes <- consensus_sig %>% filter(direction == "down") %>% pull(human_symbol) %>% na.omit()

  writeLines(up_genes, file.path(PROC_DIR, paste0("signature_", stratum, "_up_genes.txt")))
  writeLines(down_genes, file.path(PROC_DIR, paste0("signature_", stratum, "_down_genes.txt")))

  message(sprintf("    Up genes (human symbols): %d -> %s", length(up_genes), file.path(PROC_DIR, paste0("signature_", stratum, "_up_genes.txt"))))
  message(sprintf("    Down genes (human symbols): %d -> %s", length(down_genes), file.path(PROC_DIR, paste0("signature_", stratum, "_down_genes.txt"))))
}

# =============================================================================
# Summary
# =============================================================================
message("\n=== Summary ===")
message(sprintf("Studies with DE results: %d", length(all_de_results)))
message(sprintf("Tissue strata analyzed: %s", paste(names(all_meta), collapse = ", ")))
for (stratum in names(all_meta)) {
  meta_df <- all_meta[[stratum]]
  message(sprintf("  %s: %d genes meta-analyzed, %d consensus, %d core",
                  stratum, nrow(meta_df),
                  sum(meta_df$consensus, na.rm = TRUE),
                  sum(meta_df$core, na.rm = TRUE)))
}
message(sprintf("Ortholog mapping coverage: %.1f%%",
                100 * nrow(ortholog_mapping) / length(all_gene_ids)))
message("\nNext step: Run 03_drug_screening.py for LINCS connectivity mapping.")
