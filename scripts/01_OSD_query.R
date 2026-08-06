#!/usr/bin/env Rscript
# =============================================================================
# 01_OSD_query.R — Data retrieval & curation from NASA OSDR
# =============================================================================
# Queries the NASA OSDR Biological Data API (BDAPI) to enumerate all rodent
# spaceflight RNA-seq studies, selects studies for brain (neurological) and
# liver/immune (oncogenic) tissue strata, downloads processed count matrices
# and sample metadata, and cross-references human spaceflight datasets.
#
# Outputs:
#   data/raw/<OSD-ID>/          — count matrices + metadata per study
#   data/study_catalog.csv      — enumerated studies with tissue/stressor/sample annotations
#   data/human_study_catalog.csv — human spaceflight studies for validation
#
# Author: Opposite-Forcing Pipeline
# =============================================================================

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(readr)
  library(stringr)
})

# --- Configuration -----------------------------------------------------------
BDAPI_BASE <- "https://visualization.osdr.nasa.gov/biodata/api"
FILE_API_BASE <- "https://osdr.nasa.gov/osdr/data/osd/files"
PROJECT_ROOT <- "/workspace/astronaut-opposite-forcing"
RAW_DIR <- file.path(PROJECT_ROOT, "data", "raw")
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)

# Helper: robust JSON GET with retries
api_get <- function(url, params = list(), max_retries = 3) {
  for (attempt in 1:max_retries) {
    tryCatch({
      r <- GET(url, query = c(params, list(format = "json")),
               add_headers(Accept = "application/json"),
               timeout(120))
      if (status_code(r) == 200) {
        ct <- headers(r)$`content-type`
        if (!is.null(ct) && grepl("json", ct, ignore.case = TRUE)) {
          return(content(r, as = "parsed", type = "application/json"))
        }
      }
    }, error = function(e) {
      message(sprintf("  Attempt %d failed: %s", attempt, conditionMessage(e)))
    })
    if (attempt < max_retries) Sys.sleep(2^attempt)
  }
  return(NULL)
}

# =============================================================================
# Step 1: Enumerate ALL datasets in OSDR
# =============================================================================
message("\n=== Step 1: Enumerating all OSDR datasets ===")

all_datasets <- api_get(paste0(BDAPI_BASE, "/v2/datasets/"))
if (is.null(all_datasets)) stop("Failed to retrieve dataset listing from OSDR BDAPI")

osd_ids <- names(all_datasets)
message(sprintf("Total datasets in OSDR: %d", length(osd_ids)))

# =============================================================================
# Step 2: Retrieve metadata for each dataset and filter to rodent spaceflight RNA-seq
# =============================================================================
message("\n=== Step 2: Retrieving metadata and filtering to rodent spaceflight RNA-seq ===")

study_records <- list()

for (osd_id in osd_ids) {
  meta_url <- paste0(BDAPI_BASE, "/v2/dataset/", osd_id, "/metadata/")
  meta_resp <- api_get(meta_url)
  if (is.null(meta_resp)) next

  meta <- meta_resp[[osd_id]][["metadata"]]
  if (is.null(meta)) next

  # Extract organism — it's a top-level field in the metadata
  organism <- meta[["organism"]]
  if (is.null(organism)) organism <- NA
  # Handle list-valued organism
  if (is.list(organism)) organism <- paste(unlist(organism), collapse = "; ")

  # Study type / project type
  project_type <- meta[["project type"]]
  if (is.null(project_type)) project_type <- NA

  # Assay technology types
  assay_tech <- meta[["study assay technology type"]]
  if (is.null(assay_tech)) assay_tech <- list()
  has_rnaseq <- any(grepl("RNA[ -]?[Ss]eq|transcription profiling", unlist(assay_tech), ignore.case = TRUE))

  # Tissue / material type
  tissue <- meta[["material type"]]
  if (is.null(tissue)) tissue <- NA

  # Spaceflight factor
  factor_name <- meta[["study factor name"]]
  factor_type <- meta[["study factor type"]]
  # Handle list-valued factor_type
  factor_type_str <- if (is.list(factor_type)) paste(unlist(factor_type), collapse = "; ")
                    else if (is.null(factor_type)) NA else as.character(factor_type)
  is_spaceflight <- !is.na(factor_type_str) &&
    grepl("space|flight|microgravity|gravity", factor_type_str, ignore.case = TRUE)

  # Mission / project
  mission <- if (!is.null(meta[["mission"]][["name"]])) meta[["mission"]][["name"]] else NA
  project_id <- meta[["project identifier"]]
  if (is.null(project_id)) project_id <- NA

  # Title
  title <- meta[["study title"]]
  if (is.null(title)) title <- NA

  # GLDS identifier
  glds_id <- if (!is.null(meta[["identifiers"]]) && length(meta[["identifiers"]]) > 0) {
    meta[["identifiers"]][[1]]
  } else NA

  # Flight program
  flight_program <- meta[["flight program"]]
  if (is.null(flight_program)) flight_program <- NA

  # Record
  rec <- data.frame(
    osd_id = osd_id,
    glds_id = glds_id,
    title = as.character(title),
    organism = as.character(organism),
    tissue = as.character(tissue),
    project_type = as.character(project_type),
    project_id = as.character(project_id),
    mission = as.character(mission),
    flight_program = as.character(flight_program),
    factor_type = factor_type_str,
    has_rnaseq = has_rnaseq,
    is_spaceflight = is_spaceflight,
    stringsAsFactors = FALSE
  )
  study_records[[osd_id]] <- rec
}

study_df <- do.call(rbind, study_records)
rownames(study_df) <- NULL
message(sprintf("Retrieved metadata for %d datasets", nrow(study_df)))

# =============================================================================
# Step 3: Filter to rodent (mouse/rat) spaceflight RNA-seq studies
# =============================================================================
message("\n=== Step 3: Filtering to rodent spaceflight RNA-seq ===")

# Filter: organism is mouse or rat, has RNA-seq, is spaceflight
rodent_rnaseq <- study_df %>%
  filter(has_rnaseq == TRUE) %>%
  filter(grepl("mus musculus|mouse|rattus|rat", organism, ignore.case = TRUE) |
         grepl("mus musculus|mouse|rattus|rat", title, ignore.case = TRUE)) %>%
  filter(is_spaceflight == TRUE | grepl("spaceflight|space flight", project_type, ignore.case = TRUE))

message(sprintf("Rodent spaceflight RNA-seq studies: %d", nrow(rodent_rnaseq)))

# Also get human studies for validation
human_rnaseq <- study_df %>%
  filter(has_rnaseq == TRUE) %>%
  filter(grepl("homo sapiens|human", organism, ignore.case = TRUE) |
         grepl("human|astronaut", title, ignore.case = TRUE)) %>%
  filter(is_spaceflight == TRUE | grepl("spaceflight|space flight", project_type, ignore.case = TRUE))

message(sprintf("Human spaceflight RNA-seq studies: %d", nrow(human_rnaseq)))

# =============================================================================
# Step 4: Tissue-scope finalization
# =============================================================================
message("\n=== Step 4: Tissue-scope analysis ===")

# Group rodent studies by tissue
tissue_table <- rodent_rnaseq %>%
  group_by(tissue) %>%
  summarise(
    n_studies = n(),
    osd_ids = paste(osd_id, collapse = "; "),
    missions = paste(unique(mission[!is.na(mission)]), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_studies))

message("\nTissue distribution (rodent spaceflight RNA-seq):")
print(tissue_table, n = Inf)

# Define tissue strata based on availability
# Brain (neurological) and liver/immune (oncogenic) per plan
brain_tissues <- c("brain", "cerebrum", "cerebellum", "cerebral hemisphere",
                   "left cerebral hemisphere", "hippocampus", "cortex",
                   "frontal cortex", "hypothalamus")
liver_immune_tissues <- c("liver", "spleen", "thymus", "blood", "pbmc",
                          "leukocyte", "lymph node", "bone marrow",
                          "peripheral blood mononuclear cell")

# Classify each study into a stratum
classify_tissue <- function(tissue_name) {
  tn <- tolower(tissue_name)
  if (any(sapply(brain_tissues, function(bt) grepl(bt, tn, ignore.case = TRUE)))) {
    return("brain")
  }
  if (any(sapply(liver_immune_tissues, function(lt) grepl(lt, tn, ignore.case = TRUE)))) {
    return("liver_immune")
  }
  return("other")
}

rodent_rnaseq$tissue_stratum <- sapply(rodent_rnaseq$tissue, classify_tissue)

# Count studies per stratum
stratum_counts <- rodent_rnaseq %>%
  group_by(tissue_stratum) %>%
  summarise(n_studies = n(), .groups = "drop") %>%
  arrange(desc(n_studies))

message("\nStudies per tissue stratum:")
print(stratum_counts)

# Select studies: need >=2 per stratum for meta-analysis
selected_studies <- rodent_rnaseq %>%
  filter(tissue_stratum %in% c("brain", "liver_immune")) %>%
  group_by(tissue_stratum) %>%
  filter(n() >= 2) %>%
  ungroup()

message(sprintf("\nSelected studies for analysis: %d", nrow(selected_studies)))
if (nrow(selected_studies) == 0) {
  message("WARNING: No stratum has >=2 studies. Relaxing to include all rodent spaceflight RNA-seq.")
  selected_studies <- rodent_rnaseq
  selected_studies$tissue_stratum <- "pan_tissue"
}

# =============================================================================
# Step 5: Download count matrices and sample metadata for selected studies
# =============================================================================
message("\n=== Step 5: Downloading count matrices and metadata ===")

# For each selected study, find and download the STAR unnormalized counts file
# and the sample metadata (ISA-Tab or sample table)

download_study_data <- function(osd_id, glds_id, out_dir) {
  study_dir <- file.path(out_dir, osd_id)
  dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)

  # Get file listing
  files_url <- paste0(BDAPI_BASE, "/v2/dataset/", osd_id, "/files/")
  files_resp <- api_get(files_url)
  if (is.null(files_resp)) {
    message(sprintf("  WARNING: Could not retrieve files for %s", osd_id))
    return(list(osd_id = osd_id, count_file = NA, metadata_file = NA, n_samples = NA))
  }

  files <- files_resp[[osd_id]][["files"]]
  if (is.null(files)) {
    message(sprintf("  WARNING: No files found for %s", osd_id))
    return(list(osd_id = osd_id, count_file = NA, metadata_file = NA, n_samples = NA))
  }

  file_names <- names(files)

  # Find the STAR unnormalized counts file (preferred for DESeq2)
  # Pattern: *_STAR_Unnormalized_Counts_*.csv or *_RSEM_Unnormalized_Counts_*.csv
  count_file <- NULL
  star_pattern <- "STAR_Unnormalized_Counts"
  rsem_pattern <- "RSEM_Unnormalized_Counts_rRNArm"  # prefer rRNA-removed

  star_matches <- file_names[grepl(star_pattern, file_names, ignore.case = TRUE) & grepl("\\.csv$", file_names)]
  rsem_matches <- file_names[grepl(rsem_pattern, file_names, ignore.case = TRUE) & grepl("\\.csv$", file_names)]

  if (length(star_matches) > 0) {
    count_file <- star_matches[1]
  } else if (length(rsem_matches) > 0) {
    count_file <- rsem_matches[1]
  } else {
    # Fallback: any unnormalized counts
    any_counts <- file_names[grepl("Unnormalized_Counts", file_names, ignore.case = TRUE) & grepl("\\.csv$", file_names)]
    if (length(any_counts) > 0) count_file <- any_counts[1]
  }

  # Find metadata/sample table file
  # Pattern: *ISA* or *sample* or *metadata*
  meta_file <- NULL
  isa_matches <- file_names[grepl("ISA|investigation|sample_table|metadata", file_names, ignore.case = TRUE) &
                            grepl("\\.(txt|csv|tsv|zip)$", file_names)]
  if (length(isa_matches) > 0) meta_file <- isa_matches[1]

  # Download count file
  downloaded_count <- NA
  if (!is.null(count_file)) {
    file_info <- files[[count_file]]
    file_url <- file_info[["URL"]]
    if (is.null(file_url)) file_url <- file_info[["remote_url"]]
    if (!is.null(file_url)) {
      dest <- file.path(study_dir, count_file)
      tryCatch({
        download.file(file_url, dest, quiet = TRUE, mode = "wb")
        downloaded_count <- count_file
        message(sprintf("  Downloaded: %s (%s)", osd_id, count_file))
      }, error = function(e) {
        message(sprintf("  ERROR downloading %s: %s", count_file, conditionMessage(e)))
      })
    }
  } else {
    message(sprintf("  WARNING: No count matrix found for %s", osd_id))
  }

  # Download metadata file
  downloaded_meta <- NA
  if (!is.null(meta_file)) {
    file_info <- files[[meta_file]]
    file_url <- file_info[["URL"]]
    if (is.null(file_url)) file_url <- file_info[["remote_url"]]
    if (!is.null(file_url)) {
      dest <- file.path(study_dir, meta_file)
      tryCatch({
        download.file(file_url, dest, quiet = TRUE, mode = "wb")
        downloaded_meta <- meta_file
      }, error = function(e) {
        message(sprintf("  ERROR downloading metadata %s: %s", meta_file, conditionMessage(e)))
      })
    }
  }

  # Also get sample-level metadata via API
  assay_url <- paste0(BDAPI_BASE, "/v2/dataset/", osd_id, "/assay/*/sample/*/")
  assay_resp <- api_get(assay_url)
  n_samples <- NA
  if (!is.null(assay_resp) && !is.null(assay_resp[[osd_id]])) {
    samples_data <- assay_resp[[osd_id]]
    # Save sample metadata as JSON
    write_json(samples_data, file.path(study_dir, "sample_metadata.json"), auto_unbox = TRUE, pretty = TRUE)
    # Count samples
    if ("samples" %in% names(samples_data)) {
      n_samples <- length(samples_data[["samples"]])
    }
  }

  return(list(osd_id = osd_id, count_file = downloaded_count,
              metadata_file = downloaded_meta, n_samples = n_samples))
}

# Download data for each selected study
download_results <- list()
for (i in 1:nrow(selected_studies)) {
  row <- selected_studies[i, ]
  message(sprintf("\n  [%d/%d] Processing %s (%s, tissue: %s)",
                  i, nrow(selected_studies), row$osd_id, row$tissue_stratum, row$tissue))
  result <- download_study_data(row$osd_id, row$glds_id, RAW_DIR)
  result$tissue_stratum <- row$tissue_stratum
  result$tissue <- row$tissue
  result$glds_id <- row$glds_id
  result$mission <- row$mission
  download_results[[row$osd_id]] <- result
}

download_df <- do.call(rbind, lapply(download_results, as.data.frame, stringsAsFactors = FALSE))
rownames(download_df) <- NULL

# =============================================================================
# Step 6: Save study catalogs
# =============================================================================
message("\n=== Step 6: Saving study catalogs ===")

# Full study catalog (all rodent spaceflight RNA-seq)
catalog <- rodent_rnaseq %>%
  select(osd_id, glds_id, title, organism, tissue, tissue_stratum,
         project_type, project_id, mission, flight_program, factor_type) %>%
  arrange(tissue_stratum, osd_id)
write_csv(catalog, file.path(PROJECT_ROOT, "data", "study_catalog.csv"))
message(sprintf("Saved study_catalog.csv (%d studies)", nrow(catalog)))

# Download status
write_csv(download_df, file.path(PROJECT_ROOT, "data", "download_status.csv"))
message(sprintf("Saved download_status.csv (%d studies)", nrow(download_df)))

# Human study catalog for validation
if (nrow(human_rnaseq) > 0) {
  human_catalog <- human_rnaseq %>%
    select(osd_id, glds_id, title, organism, tissue,
           project_type, project_id, mission, flight_program)
  write_csv(human_catalog, file.path(PROJECT_ROOT, "data", "human_study_catalog.csv"))
  message(sprintf("Saved human_study_catalog.csv (%d studies)", nrow(human_catalog)))
} else {
  message("No human spaceflight RNA-seq studies found for validation.")
}

# =============================================================================
# Summary
# =============================================================================
message("\n=== Summary ===")
message(sprintf("Total OSDR datasets scanned: %d", length(osd_ids)))
message(sprintf("Rodent spaceflight RNA-seq studies found: %d", nrow(rodent_rnaseq)))
message(sprintf("Studies selected for analysis: %d", nrow(selected_studies)))
message(sprintf("  Brain stratum: %d studies", sum(selected_studies$tissue_stratum == "brain")))
message(sprintf("  Liver/immune stratum: %d studies", sum(selected_studies$tissue_stratum == "liver_immune")))
message(sprintf("Human validation studies: %d", nrow(human_rnaseq)))
message(sprintf("Count matrices downloaded: %d", sum(!is.na(download_df$count_file))))
message("\nNext step: Run 02_biomarker_identification.R for DE analysis and consensus signature.")
