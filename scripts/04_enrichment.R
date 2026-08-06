#!/usr/bin/env Rscript
# =============================================================================
# 04_enrichment.R — Disease & functional enrichment of spaceflight signatures
# =============================================================================
# 1. GSEA on ranked meta-analysis statistics (Hallmark gene sets)
# 2. ORA on consensus up/down genes (GO:BP + Reactome)
# 3. JensenLab DISEASES enrichment (disease-gene associations)
#
# KEGG is excluded (commercial license). Hallmark + Reactome + GO:BP only.
# =============================================================================

.libPaths(c("/workspace/.Rlib", .libPaths()))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(AnnotationDbi)
  library(jsonlite)
})

set.seed(42)

PROJECT_ROOT <- "/workspace/astronaut-opposite-forcing"
PROC_DIR     <- file.path(PROJECT_ROOT, "data", "processed")
RESULTS_DIR  <- file.path(PROJECT_ROOT, "results", "tables")
FIGURES_DIR  <- file.path(PROJECT_ROOT, "results", "figures")
os.makedirs  <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
os.makedirs(RESULTS_DIR)
os.makedirs(FIGURES_DIR)

# =============================================================================
# Load meta-analysis results and build ranked lists
# =============================================================================
message("=== Loading meta-analysis results ===")

load_stratum <- function(stratum) {
  meta_file <- file.path(PROC_DIR, paste0("meta_analysis_", stratum, "_full.csv"))
  if (!file.exists(meta_file)) {
    message(sprintf("  WARNING: %s not found", meta_file))
    return(NULL)
  }
  meta_df <- read_csv(meta_file, show_col_types = FALSE)

  # Get human symbols
  ortho_file <- file.path(PROC_DIR, "ortholog_mapping.csv")
  ortho <- read_csv(ortho_file, show_col_types = FALSE)
  meta_df <- meta_df %>%
    dplyr::left_join(ortho %>% dplyr::select(mouse_ensembl_id, human_symbol),
              by = c("gene_id" = "mouse_ensembl_id"))

  # Filter to genes with human orthologs and valid stats
  meta_df <- meta_df %>%
    filter(!is.na(human_symbol), human_symbol != "",
           !is.na(pooled_log2FC), !is.na(pvalue))

  # Remove duplicates (keep most significant)
  meta_df <- meta_df %>%
    arrange(pvalue) %>%
    distinct(human_symbol, .keep_all = TRUE)

  # Ranked list for GSEA: sign(log2FC) * -log10(pvalue)
  meta_df$rank_score <- sign(meta_df$pooled_log2FC) * -log10(meta_df$pvalue + 1e-300)
  ranked_df <- meta_df %>%
    arrange(desc(rank_score)) %>%
    dplyr::select(human_symbol, rank_score)
  ranked <- setNames(ranked_df$rank_score, ranked_df$human_symbol)

  # Consensus up/down gene lists for ORA
  consensus <- meta_df %>% filter(consensus == TRUE)
  up_genes <- consensus %>% filter(direction == "up") %>% pull(human_symbol) %>% na.omit()
  dn_genes <- consensus %>% filter(direction == "down") %>% pull(human_symbol) %>% na.omit()

  message(sprintf("  %s: %d ranked genes, %d up, %d down consensus",
                  stratum, length(ranked), length(up_genes), length(dn_genes)))

  return(list(ranked = ranked, up_genes = up_genes, dn_genes = dn_genes,
              meta_df = meta_df, stratum = stratum))
}

strata_data <- list()
for (s in c("liver_immune", "brain")) {
  strata_data[[s]] <- load_stratum(s)
}

# =============================================================================
# Prepare gene sets from msigdbr
# =============================================================================
message("\n=== Preparing gene sets ===")

# Hallmark
h_df <- msigdbr(collection = "H")
hallmark_terms <- split(h_df$gene_symbol, h_df$gs_name)
message(sprintf("  Hallmark: %d gene sets", length(hallmark_terms)))

# Reactome
r_df <- msigdbr(collection = "C2", subcollection = "CP:REACTOME")
reactome_terms <- split(r_df$gene_symbol, r_df$gs_name)
message(sprintf("  Reactome: %d gene sets", length(reactome_terms)))

# GO:BP
go_df <- msigdbr(collection = "C5", subcollection = "GO:BP")
go_terms <- split(go_df$gene_symbol, go_df$gs_name)
message(sprintf("  GO:BP: %d gene sets", length(go_terms)))

# =============================================================================
# GSEA on Hallmark gene sets
# =============================================================================
message("\n=== GSEA (Hallmark) ===")

run_gsea <- function(ranked, term2gene, stratum, collection_name) {
  if (length(ranked) < 10) {
    message(sprintf("  Skipping GSEA for %s (%s): too few ranked genes", stratum, collection_name))
    return(NULL)
  }

  # Prepare term2gene for clusterProfiler
  t2g <- data.frame(
    term = names(term2gene),
    gene = unlist(lapply(term2gene, function(x) paste(x, collapse = ";")),
                  use.names = FALSE)
  )

  set.seed(42)
  gsea_res <- tryCatch({
    GSEA(geneList = ranked,
         TERM2GENE = t2g %>% tidyr::separate_rows(gene, sep = ";") %>%
           dplyr::rename(GENE = gene, TERM = term),
         pvalueCutoff = 0.25,  # relaxed for exploration
         pAdjustMethod = "BH",
         verbose = FALSE,
         seed = TRUE,
         minGSSize = 10,
         maxGSSize = 500)
  }, error = function(e) {
    message(sprintf("  GSEA error: %s", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
    message(sprintf("  No significant GSEA hits for %s (%s)", stratum, collection_name))
    return(NULL)
  }

  res_df <- as.data.frame(gsea_res)
  message(sprintf("  %s (%s): %d significant pathways (p.adj < 0.25)",
                  stratum, collection_name, nrow(res_df)))

  return(res_df)
}

# Run GSEA for each stratum × collection
all_gsea <- list()
for (s in names(strata_data)) {
  if (is.null(strata_data[[s]])) next
  ranked <- strata_data[[s]]$ranked

  all_gsea[[paste0(s, "_hallmark")]] <- run_gsea(ranked, hallmark_terms, s, "Hallmark")
  all_gsea[[paste0(s, "_reactome")]] <- run_gsea(ranked, reactome_terms, s, "Reactome")

  # Save
  for (coll in c("hallmark", "reactome")) {
    key <- paste0(s, "_", coll)
    if (!is.null(all_gsea[[key]])) {
      write_csv(all_gsea[[key]],
                file.path(RESULTS_DIR, paste0("gsea_", s, "_", coll, ".csv")))
    }
  }
}

# =============================================================================
# ORA on consensus up/down genes
# =============================================================================
message("\n=== ORA (GO:BP + Reactome) ===")

run_ora <- function(genes, term2gene, stratum, direction, collection_name) {
  if (length(genes) < 5) {
    message(sprintf("  Skipping ORA for %s %s (%s): too few genes (%d)",
                    stratum, direction, collection_name, length(genes)))
    return(NULL)
  }

  t2g <- data.frame(
    TERM = rep(names(term2gene), lengths(term2gene)),
    GENE = unlist(term2gene)
  )

  set.seed(42)
  ora_res <- tryCatch({
    enrichGO(gene = genes,
             OrgDb = org.Hs.eg.db,
             keyType = "SYMBOL",
             ont = "BP",
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             qvalueCutoff = 0.2)
  }, error = function(e) {
    # Fallback to generic enricher for non-GO collections
    tryCatch({
      enricher(gene = genes, TERM2GENE = t2g,
               pvalueCutoff = 0.05, pAdjustMethod = "BH")
    }, error = function(e2) {
      message(sprintf("  ORA error: %s", conditionMessage(e2)))
      return(NULL)
    })
  })

  if (is.null(ora_res) || nrow(as.data.frame(ora_res)) == 0) {
    message(sprintf("  No significant ORA hits for %s %s (%s)",
                    stratum, direction, collection_name))
    return(NULL)
  }

  res_df <- as.data.frame(ora_res)
  message(sprintf("  %s %s (%s): %d enriched terms (p.adj < 0.05)",
                  stratum, direction, collection_name, nrow(res_df)))
  return(res_df)
}

# Run ORA for each stratum × direction × collection
all_ora <- list()
for (s in names(strata_data)) {
  if (is.null(strata_data[[s]])) next

  for (dir in c("up", "dn")) {
    genes <- if (dir == "up") strata_data[[s]]$up_genes else strata_data[[s]]$dn_genes

    # GO:BP via enrichGO
    key_go <- paste0(s, "_", dir, "_gobp")
    all_ora[[key_go]] <- run_ora(genes, go_terms, s, dir, "GO:BP")
    if (!is.null(all_ora[[key_go]])) {
      write_csv(all_ora[[key_go]],
                file.path(RESULTS_DIR, paste0("ora_", s, "_", dir, "_gobp.csv")))
    }

    # Reactome via enricher
    key_react <- paste0(s, "_", dir, "_reactome")
    t2g_react <- data.frame(
      TERM = rep(names(reactome_terms), lengths(reactome_terms)),
      GENE = unlist(reactome_terms)
    )
    if (length(genes) >= 5) {
      set.seed(42)
      react_res <- tryCatch({
        enricher(gene = genes, TERM2GENE = t2g_react,
                 pvalueCutoff = 0.05, pAdjustMethod = "BH")
      }, error = function(e) NULL)
      all_ora[[key_react]] <- if (!is.null(react_res)) as.data.frame(react_res) else NULL
      if (!is.null(all_ora[[key_react]])) {
        write_csv(all_ora[[key_react]],
                  file.path(RESULTS_DIR, paste0("ora_", s, "_", dir, "_reactome.csv")))
        message(sprintf("  %s %s (Reactome): %d enriched terms",
                        s, dir, nrow(all_ora[[key_react]])))
      }
    }
  }
}

# =============================================================================
# JensenLab DISEASES enrichment
# =============================================================================
message("\n=== JensenLab DISEASES enrichment ===")

# Download JensenLab DISEASES knowledge channel (filtered, curated)
# URL: https://download.jensenlab.org/human_disease_knowledge_filtered.tsv
jensen_file <- file.path(PROJECT_ROOT, "data", "raw", "jensenlab_diseases_knowledge_filtered.tsv")
jensen_url <- "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv"

if (!file.exists(jensen_file)) {
  message(sprintf("  Downloading JensenLab DISEASES knowledge (filtered) from %s...", jensen_url))
  tryCatch({
    download.file(jensen_url, jensen_file, quiet = TRUE, timeout = 300)
    message(sprintf("  Downloaded: %.1f KB", file.size(jensen_file) / 1e3))
  }, error = function(e) {
    message(sprintf("  WARNING: Download failed: %s", conditionMessage(e)))
  })
}

run_jensen_enrichment <- function(up_genes, dn_genes, stratum) {
  if (!file.exists(jensen_file)) {
    message(sprintf("  JensenLab file not available for %s", stratum))
    return(NULL)
  }

  # Read JensenLab knowledge: gene_ensp, gene_symbol, disease_doid, disease_name, source, evidence_type, confidence
  jensen <- tryCatch({
    read_tsv(jensen_file,
             col_names = c("gene_id", "protein", "disease_id", "disease", "source", "evidence_type", "score"),
             show_col_types = FALSE, progress = FALSE)
  }, error = function(e) {
    message(sprintf("  Error reading JensenLab: %s", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(jensen)) return(NULL)
  message(sprintf("  JensenLab knowledge: %d disease-gene associations, %d diseases",
                  nrow(jensen), length(unique(jensen$disease))))

  # Build disease → gene set (using gene symbols, filter by confidence)
  disease_genes <- jensen %>%
    filter(score >= 3) %>%  # confidence filter (3+ = high confidence)
    group_by(disease) %>%
    summarise(genes = list(unique(protein)), n_genes = n(), .groups = "drop")

  # Enrichment via Fisher's exact for each disease
  all_sig_genes <- unique(c(up_genes, dn_genes))
  universe <- unique(jensen$protein)

  results <- list()
  for (i in seq_len(nrow(disease_genes))) {
    disease <- disease_genes$disease[i]
    disease_set <- unlist(disease_genes$genes[i])
    overlap <- intersect(all_sig_genes, disease_set)

    if (length(overlap) < 2) next

    # Fisher's exact test
    a <- length(overlap)
    b <- length(all_sig_genes) - a
    c <- length(disease_set) - a
    d <- length(universe) - length(disease_set) - b

    pval <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value

    if (pval < 0.01) {
      # Separate up/down contributions
      up_overlap <- intersect(up_genes, disease_set)
      dn_overlap <- intersect(dn_genes, disease_set)

      results[[length(results) + 1]] <- data.frame(
        disease = disease,
        n_overlap = length(overlap),
        n_disease_genes = length(disease_set),
        pvalue = pval,
        padj = NA,  # will adjust later
        up_genes = paste(up_overlap, collapse = ";"),
        n_up = length(up_overlap),
        down_genes = paste(dn_overlap, collapse = ";"),
        n_down = length(dn_overlap),
        direction = ifelse(length(up_overlap) > length(dn_overlap), "up_enriched", "down_enriched")
      )
    }
  }

  if (length(results) == 0) {
    message(sprintf("  No disease enrichment for %s", stratum))
    return(NULL)
  }

  res_df <- do.call(rbind, results)
  res_df$padj <- p.adjust(res_df$pvalue, method = "BH")
  res_df <- res_df %>% arrange(pvalue)
  res_df$stratum <- stratum

  message(sprintf("  %s: %d diseases enriched (p < 0.01, before FDR)", stratum, nrow(res_df)))
  return(res_df)
}

all_jensen <- list()
for (s in names(strata_data)) {
  if (is.null(strata_data[[s]])) next
  all_jensen[[s]] <- run_jensen_enrichment(
    strata_data[[s]]$up_genes,
    strata_data[[s]]$dn_genes,
    s
  )
  if (!is.null(all_jensen[[s]])) {
    write_csv(all_jensen[[s]],
              file.path(RESULTS_DIR, paste0("jensen_diseases_", s, ".csv")))
  }
}

# =============================================================================
# Summary
# =============================================================================
message("\n=== Enrichment Summary ===")
for (s in names(strata_data)) {
  if (is.null(strata_data[[s]])) next
  message(sprintf("\n  %s:", s))
  for (coll in c("hallmark", "reactome")) {
    key <- paste0(s, "_", coll)
    n <- if (!is.null(all_gsea[[key]])) nrow(all_gsea[[key]]) else 0
    message(sprintf("    GSEA %s: %d pathways", coll, n))
  }
  for (dir in c("up", "dn")) {
    for (coll in c("gobp", "reactome")) {
      key <- paste0(s, "_", dir, "_", coll)
      n <- if (!is.null(all_ora[[key]])) nrow(all_ora[[key]]) else 0
      message(sprintf("    ORA %s %s: %d terms", dir, coll, n))
    }
  }
  n_jensen <- if (!is.null(all_jensen[[s]])) nrow(all_jensen[[s]]) else 0
  message(sprintf("    JensenLab diseases: %d", n_jensen))
}

message("\nNext step: Run 05_nutrient_gene_mapping.py for food-compound mapping.")
