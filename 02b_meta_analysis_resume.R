#!/usr/bin/env Rscript
# =============================================================================
# 02b_meta_analysis_resume.R — Resume from per-study DE results
# =============================================================================
# Skips DESeq2 (already completed), loads per-study DE results, runs
# meta-analysis, ortholog mapping, and builds consensus signatures.
#
# Fixes applied (v2):
#   1. ERCC spike-in removal BEFORE meta-analysis
#   2. Direction consistency relaxed to >=70% majority (not 100%)
#   3. Local ortholog mapping via org.Mm.eg.db + homologene (no biomaRt)
#   4. Empty gene list safety for writeLines
# =============================================================================

.libPaths(c("/workspace/.Rlib", .libPaths()))

suppressPackageStartupMessages({
  library(metafor)
  library(dplyr)
  library(readr)
  library(stringr)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(homologene)
  library(AnnotationDbi)
})

set.seed(42)

PROJECT_ROOT <- "/workspace/astronaut-opposite-forcing"
PROC_DIR <- file.path(PROJECT_ROOT, "data", "processed")
PER_STUDY_DIR <- file.path(PROC_DIR, "per_study_de")

# --- Load study catalog ---
study_catalog <- read_csv(file.path(PROJECT_ROOT, "data", "study_catalog.csv"), show_col_types = FALSE)
download_status <- read_csv(file.path(PROJECT_ROOT, "data", "download_status.csv"), show_col_types = FALSE)

# --- Load all per-study DE results ---
de_files <- list.files(PER_STUDY_DIR, pattern = "_de_results\\.csv$", full.names = TRUE)
message(sprintf("Loading %d per-study DE result files...", length(de_files)))

all_de_results <- list()
for (f in de_files) {
  osd_id <- str_extract(basename(f), "OSD-\\d+")
  de_res <- read_csv(f, show_col_types = FALSE)

  # Remove ERCC spike-in controls (synthetic RNA, not real genes)
  n_ercc <- sum(grepl("^ERCC", de_res$gene_id))
  de_res <- de_res[!grepl("^ERCC", de_res$gene_id), ]

  # Quality filter: skip studies with over-shrunk LFCs (median SE < 0.01)
  median_se <- median(de_res$SE, na.rm = TRUE)
  if (median_se < 0.01) {
    message(sprintf("  Skipping %s (median SE = %.4f, over-shrunk; removed %d ERCC first)",
                    osd_id, median_se, n_ercc))
    next
  }

  de_res$osd_id <- osd_id
  message(sprintf("  Loaded %s: %d genes (removed %d ERCC), median SE = %.4f",
                  osd_id, nrow(de_res), n_ercc, median_se))

  # Assign tissue stratum from catalog
  stratum <- study_catalog$tissue_stratum[study_catalog$osd_id == osd_id]
  if (length(stratum) > 0) {
    de_res$tissue_stratum <- stratum[1]
  } else {
    stratum <- download_status$tissue_stratum[download_status$osd_id == osd_id]
    de_res$tissue_stratum <- if (length(stratum) > 0) stratum[1] else NA
  }

  all_de_results[[osd_id]] <- de_res
}

de_combined <- do.call(rbind, all_de_results)
rownames(de_combined) <- NULL
message(sprintf("Combined DE results: %d rows from %d studies (ERCC-free)",
                nrow(de_combined), length(all_de_results)))

# =============================================================================
# Meta-analysis per tissue stratum
# =============================================================================
message("\n=== Meta-analysis ===")

run_meta_analysis <- function(de_combined, stratum) {
  stratum_de <- de_combined %>% filter(tissue_stratum == stratum)
  if (nrow(stratum_de) == 0) {
    message(sprintf("  No data for stratum: %s", stratum))
    return(NULL)
  }

  gene_study_counts <- stratum_de %>%
    filter(!is.na(log2FC), !is.na(SE), SE > 0) %>%
    group_by(gene_id) %>%
    summarise(n_studies = n(), .groups = "drop")

  meta_genes <- gene_study_counts %>% filter(n_studies >= 2) %>% pull(gene_id)
  message(sprintf("  %s: %d genes in >=2 studies for meta-analysis", stratum, length(meta_genes)))

  if (length(meta_genes) == 0) return(NULL)

  # Run meta-analysis per gene
  meta_results <- vector("list", length(meta_genes))
  for (idx in seq_along(meta_genes)) {
    gene <- meta_genes[idx]
    gene_data <- stratum_de %>%
      filter(gene_id == gene, !is.na(log2FC), !is.na(SE), SE > 0)

    if (nrow(gene_data) < 2) next

    tryCatch({
      rma_fit <- rma(yi = log2FC, sei = SE, data = gene_data, method = "REML",
                     control = list(stepadj = 0.5, maxiter = 100))

      directions <- sign(gene_data$log2FC)
      n_up <- sum(directions > 0)
      n_down <- sum(directions < 0)
      k_total <- n_up + n_down

      # Relaxed direction consistency: >=70% of studies agree on direction
      # (100% is too strict across heterogeneous tissues/missions/durations)
      direction_majority <- max(n_up, n_down) / k_total >= 0.70

      meta_results[[idx]] <- data.frame(
        gene_id = gene,
        pooled_log2FC = rma_fit$beta[1],
        pooled_SE = rma_fit$se,
        z = rma_fit$zval,
        pvalue = rma_fit$pval,
        I2 = rma_fit$I2,
        tau2 = rma_fit$tau2,
        k = rma_fit$k,
        n_studies_up = n_up,
        n_studies_down = n_down,
        direction_consistent = direction_majority,
        per_study_log2FC = paste(sprintf("%.3f", gene_data$log2FC), collapse = ";"),
        per_study_osd = paste(gene_data$osd_id, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      tryCatch({
        rma_fit <- rma(yi = log2FC, sei = SE, data = gene_data, method = "DL")
        directions <- sign(gene_data$log2FC)
        n_up <- sum(directions > 0)
        n_down <- sum(directions < 0)
        k_total <- n_up + n_down
        direction_majority <- max(n_up, n_down) / k_total >= 0.70

        meta_results[[idx]] <<- data.frame(
          gene_id = gene,
          pooled_log2FC = rma_fit$beta[1],
          pooled_SE = rma_fit$se,
          z = rma_fit$zval,
          pvalue = rma_fit$pval,
          I2 = rma_fit$I2,
          tau2 = rma_fit$tau2,
          k = rma_fit$k,
          n_studies_up = n_up,
          n_studies_down = n_down,
          direction_consistent = direction_majority,
          per_study_log2FC = paste(sprintf("%.3f", gene_data$log2FC), collapse = ";"),
          per_study_osd = paste(gene_data$osd_id, collapse = ";"),
          stringsAsFactors = FALSE
        )
      }, error = function(e2) NULL)
    })
  }

  meta_df <- do.call(rbind, meta_results)
  if (is.null(meta_df) || nrow(meta_df) == 0) return(NULL)

  meta_df$padj <- p.adjust(meta_df$pvalue, method = "BH")

  # Consensus: nominal p < 0.05 AND >=70% direction consistency
  meta_df$consensus <- meta_df$pvalue < 0.05 & meta_df$direction_consistent
  meta_df$core <- meta_df$consensus & abs(meta_df$pooled_log2FC) >= 0.5
  meta_df$direction <- ifelse(meta_df$pooled_log2FC > 0, "up", "down")

  return(meta_df)
}

strata <- unique(de_combined$tissue_stratum[!is.na(de_combined$tissue_stratum)])
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
# Technical gene filtering
# =============================================================================
message("\n=== Technical gene filtering ===")

# Filter by mouse symbol (will map after)
# We need to get mouse symbols first for pattern matching
all_gene_ids <- unique(unlist(lapply(all_meta, function(m) m$gene_id)))
message(sprintf("Getting mouse symbols for %d genes...", length(all_gene_ids)))

mouse_symbols <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = all_gene_ids,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
) %>% distinct(ENSEMBL, .keep_all = TRUE)

symbol_lookup <- setNames(mouse_symbols$SYMBOL, mouse_symbols$ENSEMBL)

technical_patterns <- c("^Rpl\\d", "^Rps\\d", "GAPDH", "ACTB", "B2M", "ACTG1", "^Eef\\d", "^Eif\\d")

filter_technical_genes <- function(meta_df, symbol_lookup) {
  meta_df$mouse_symbol <- symbol_lookup[meta_df$gene_id]

  is_technical <- sapply(meta_df$mouse_symbol, function(sym) {
    if (is.na(sym)) return(FALSE)
    any(sapply(technical_patterns, function(p) grepl(p, sym, ignore.case = TRUE)))
  })
  is_metallothionein <- grepl("^MT1|^MT2A|^Mt1|^Mt2a", meta_df$mouse_symbol, ignore.case = TRUE)
  is_technical <- is_technical & !is_metallothionein

  message(sprintf("  Removed %d technical genes (retained %d metallothioneins)",
                  sum(is_technical, na.rm = TRUE), sum(is_metallothionein, na.rm = TRUE)))
  meta_df[!is_technical, ]
}

for (stratum in names(all_meta)) {
  all_meta[[stratum]] <- filter_technical_genes(all_meta[[stratum]], symbol_lookup)
}

# =============================================================================
# Mouse -> Human ortholog mapping (local, no biomaRt)
# =============================================================================
message("\n=== Mouse -> Human ortholog mapping (local) ===")

# Rebuild gene list after filtering
all_gene_ids_filtered <- unique(unlist(lapply(all_meta, function(m) m$gene_id)))
message(sprintf("Total unique gene IDs to map: %d", length(all_gene_ids_filtered)))

# Step 1: Ensembl -> mouse ENTREZID
mouse_entrez <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = all_gene_ids_filtered,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "ENTREZID", "SYMBOL")
) %>% distinct(ENSEMBL, .keep_all = TRUE)

message(sprintf("  Step 1: %d/%d Ensembl IDs mapped to mouse Entrez",
                sum(!is.na(mouse_entrez$ENTREZID)), nrow(mouse_entrez)))

# Step 2: homologene mouse (10090) -> human (9606)
entrez_to_map <- mouse_entrez$ENTREZID[!is.na(mouse_entrez$ENTREZID)]
homo_result <- homologene::homologene(entrez_to_map, inTax = 10090, outTax = 9606)
colnames(homo_result) <- c("mouse_symbol_h", "human_symbol", "mouse_entrez", "human_entrez")

message(sprintf("  Step 2: %d mouse Entrez IDs mapped to human orthologs via homologene",
                nrow(homo_result)))

# Step 3: Build final mapping table
# NOTE: use dplyr::select explicitly to avoid AnnotationDbi::select conflict
# Cast types for join compatibility
homo_result$mouse_entrez <- as.character(homo_result$mouse_entrez)
homo_result$human_entrez <- as.character(homo_result$human_entrez)
mouse_entrez$ENTREZID <- as.character(mouse_entrez$ENTREZID)

ortholog_mapping <- mouse_entrez %>%
  dplyr::left_join(homo_result %>% dplyr::select(mouse_entrez, human_symbol, human_entrez),
            by = c("ENTREZID" = "mouse_entrez")) %>%
  dplyr::rename(mouse_ensembl_id = ENSEMBL,
         mouse_symbol = SYMBOL,
         mouse_entrez_id = ENTREZID,
         human_entrez_id = human_entrez) %>%
  # Deduplicate: keep first human symbol per mouse Ensembl ID
  dplyr::group_by(mouse_ensembl_id) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

n_mapped <- sum(!is.na(ortholog_mapping$human_symbol))
message(sprintf("  Final: %d/%d mouse genes mapped to human symbols (%.1f%% coverage)",
                n_mapped, nrow(ortholog_mapping),
                100 * n_mapped / nrow(ortholog_mapping)))

write_csv(ortholog_mapping, file.path(PROC_DIR, "ortholog_mapping.csv"))

# =============================================================================
# Build consensus signatures with human orthologs
# =============================================================================
message("\n=== Building consensus signatures ===")

for (stratum in names(all_meta)) {
  meta_df <- all_meta[[stratum]]

  meta_df <- meta_df %>%
    dplyr::left_join(ortholog_mapping %>% dplyr::select(mouse_ensembl_id, human_symbol),
              by = c("gene_id" = "mouse_ensembl_id"))

  consensus_sig <- meta_df %>%
    dplyr::filter(consensus == TRUE) %>%
    dplyr::select(gene_id, mouse_symbol, human_symbol, pooled_log2FC, pvalue, padj,
           direction, I2, k, n_studies_up, n_studies_down, core) %>%
    dplyr::arrange(desc(abs(pooled_log2FC)))

  write_csv(consensus_sig, file.path(PROC_DIR, paste0("consensus_signature_", stratum, ".csv")))

  n_up <- sum(consensus_sig$direction == "up")
  n_down <- sum(consensus_sig$direction == "down")
  n_human <- sum(!is.na(consensus_sig$human_symbol) & consensus_sig$human_symbol != "")

  message(sprintf("  %s consensus signature: %d genes (%d up, %d down), %d with human orthologs",
                  stratum, nrow(consensus_sig), n_up, n_down, n_human))

  # Write gene lists for LINCS (human symbols only)
  up_genes <- consensus_sig %>% filter(direction == "up") %>% pull(human_symbol) %>% na.omit()
  down_genes <- consensus_sig %>% filter(direction == "down") %>% pull(human_symbol) %>% na.omit()

  up_file <- file.path(PROC_DIR, paste0("signature_", stratum, "_up_genes.txt"))
  down_file <- file.path(PROC_DIR, paste0("signature_", stratum, "_down_genes.txt"))

  if (length(up_genes) > 0) {
    writeLines(up_genes, up_file)
    message(sprintf("    Wrote %d up genes to %s", length(up_genes), basename(up_file)))
  } else {
    writeLines("", up_file)
    message(sprintf("    WARNING: 0 up genes for %s — wrote empty file", stratum))
  }

  if (length(down_genes) > 0) {
    writeLines(down_genes, down_file)
    message(sprintf("    Wrote %d down genes to %s", length(down_genes), basename(down_file)))
  } else {
    writeLines("", down_file)
    message(sprintf("    WARNING: 0 down genes for %s — wrote empty file", stratum))
  }
}

# =============================================================================
# Summary
# =============================================================================
message("\n=== Summary ===")
message(sprintf("Studies with DE results (post-filter): %d", length(all_de_results)))
message(sprintf("Tissue strata analyzed: %s", paste(names(all_meta), collapse = ", ")))
for (stratum in names(all_meta)) {
  meta_df <- all_meta[[stratum]]
  message(sprintf("  %s: %d genes meta-analyzed, %d consensus, %d core",
                  stratum, nrow(meta_df),
                  sum(meta_df$consensus, na.rm = TRUE),
                  sum(meta_df$core, na.rm = TRUE)))
}
message(sprintf("Ortholog mapping coverage: %.1f%%",
                100 * n_mapped / nrow(ortholog_mapping)))
message("\nNext step: Run 03_drug_screening.py for LINCS connectivity mapping.")
