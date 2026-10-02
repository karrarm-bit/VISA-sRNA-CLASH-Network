# ============================================================
# BacRegRNA Project
# Script 13B — RsaOI Known-Target Resolution
#
# PURPOSE
# ------------------------------------------------------------
# 1. Resolve known RsaOI targets using JKD6008 annotation.
# 2. Search CLASH edges using:
#       gene symbol
#       locus tag
#       RefSeq protein accession
#       available aliases
# 3. Specifically benchmark:
#       RsaOI -> atl
#       RsaOI -> HPr
# 4. Determine whether known interactions are:
#       A) present in CLASH and retained by DEG filter
#       B) present in CLASH but excluded by DEG filter
#       C) absent from this CLASH edge table
#
# IMPORTANT:
# This script does NOT change thresholds.
# It does NOT modify Scripts 13 or 14.
# ============================================================


# ------------------------------------------------------------
# 1. PROJECT PATHS
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

clash_file <- file.path(
  project_folder,
  "13_CLASH_network",
  "tables",
  "GSE254532_sRNA_mRNA_edges_collapsed.csv"
)

kb_file <- file.path(
  project_folder,
  "12_master_knowledgebase",
  "tables",
  "JKD6008_master_gene_knowledgebase.csv"
)

output_folder <- file.path(
  project_folder,
  "13B_RsaOI_target_resolution"
)

table_folder <- file.path(
  output_folder,
  "tables"
)

dir.create(
  table_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "tidyr",
  "tibble"
)

installed_names <- rownames(
  installed.packages()
)

missing_packages <- setdiff(
  required_packages,
  installed_names
)

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)


# ------------------------------------------------------------
# 3. VERIFY INPUT FILES
# ------------------------------------------------------------

if (!file.exists(clash_file)) {
  stop(
    paste0(
      "CLASH file not found:\n",
      clash_file
    )
  )
}

if (!file.exists(kb_file)) {
  stop(
    paste0(
      "Knowledgebase file not found:\n",
      kb_file
    )
  )
}


# ------------------------------------------------------------
# 4. READ DATA
# ------------------------------------------------------------

edges <- readr::read_csv(
  clash_file,
  show_col_types = FALSE
)

kb <- readr::read_csv(
  kb_file,
  show_col_types = FALSE
)


cat(
  "\n========================================\n",
  "SCRIPT 13B — RsaOI TARGET RESOLUTION\n",
  "========================================\n"
)

cat(
  "\nCLASH edges:",
  nrow(edges),
  "\n"
)

cat(
  "Knowledgebase genes:",
  nrow(kb),
  "\n"
)


# ------------------------------------------------------------
# 5. SAFE CONVERSION FUNCTIONS
# ------------------------------------------------------------

safe_chr <- function(x) {
  
  x <- as.character(x)
  
  x[is.na(x)] <- ""
  
  x
}


safe_num <- function(x) {
  
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


safe_logical <- function(x) {
  
  if (is.logical(x)) {
    
    x[is.na(x)] <- FALSE
    
    return(x)
  }
  
  x <- str_to_lower(
    str_squish(
      as.character(x)
    )
  )
  
  x %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}


# ------------------------------------------------------------
# 6. FIND RELEVANT COLUMN NAMES AUTOMATICALLY
# ------------------------------------------------------------

cat(
  "\nKnowledgebase columns:\n"
)

print(
  names(kb)
)


cat(
  "\nCLASH columns:\n"
)

print(
  names(edges)
)


# ------------------------------------------------------------
# 7. BUILD SEARCHABLE KNOWLEDGEBASE TEXT
#
# This deliberately searches across ALL columns.
# This prevents failure caused by different annotation columns.
# ------------------------------------------------------------

kb_search <- apply(
  kb,
  1,
  function(row) {
    
    paste(
      safe_chr(row),
      collapse = " | "
    )
  }
)

kb$all_annotation_text <- str_to_lower(
  kb_search
)


# ------------------------------------------------------------
# 8. DEFINE KNOWN TARGET IDENTIFIERS
#
# These identifiers were resolved from the JKD6008
# knowledgebase.
#
# atl:
#   SAA6008_RS05315
#   WP_001074557.1
#
# HPr:
#   SAA6008_RS05475
#   WP_000437472.1
# ------------------------------------------------------------

known_targets <- tibble(
  
  benchmark_target = c(
    "atl",
    "HPr"
  ),
  
  locus_tag = c(
    "SAA6008_RS05315",
    "SAA6008_RS05475"
  ),
  
  refseq_protein = c(
    "WP_001074557.1",
    "WP_000437472.1"
  ),
  
  biological_role = c(
    "major autolysin",
    "phosphocarrier protein HPr"
  )
)


print(
  known_targets
)


# ------------------------------------------------------------
# 9. VERIFY TARGETS AGAINST KNOWLEDGEBASE
# ------------------------------------------------------------

atl_kb <- kb |>
  filter(
    str_detect(
      all_annotation_text,
      regex(
        "SAA6008_RS05315|WP_001074557\\.1|\\batl\\b",
        ignore_case = TRUE
      )
    )
  )


hpr_kb <- kb |>
  filter(
    str_detect(
      all_annotation_text,
      regex(
        "SAA6008_RS05475|WP_000437472\\.1|phosphocarrier protein HPr",
        ignore_case = TRUE
      )
    )
  )


cat(
  "\n----------------------------------------\n",
  "KNOWLEDGEBASE VERIFICATION\n",
  "----------------------------------------\n"
)

cat(
  "atl annotation rows:",
  nrow(atl_kb),
  "\n"
)

cat(
  "HPr annotation rows:",
  nrow(hpr_kb),
  "\n"
)


readr::write_csv(
  atl_kb,
  file.path(
    table_folder,
    "atl_knowledgebase_annotation.csv"
  )
)


readr::write_csv(
  hpr_kb,
  file.path(
    table_folder,
    "HPr_knowledgebase_annotation.csv"
  )
)


# ------------------------------------------------------------
# 10. BUILD SEARCHABLE CLASH TEXT
#
# Search ALL columns rather than relying on one target field.
# ------------------------------------------------------------

edge_search <- apply(
  edges,
  1,
  function(row) {
    
    paste(
      safe_chr(row),
      collapse = " | "
    )
  }
)

edges$all_edge_text <- str_to_lower(
  edge_search
)


# ------------------------------------------------------------
# 11. IDENTIFY RsaOI EDGES
# ------------------------------------------------------------

rsaoi_edges <- edges |>
  filter(
    str_detect(
      all_edge_text,
      regex(
        "rsaoi|sau[-_ ]?6477",
        ignore_case = TRUE
      )
    )
  )


cat(
  "\n----------------------------------------\n",
  "RsaOI RECOVERY\n",
  "----------------------------------------\n"
)

cat(
  "Total RsaOI-associated CLASH edges:",
  nrow(rsaoi_edges),
  "\n"
)


# ------------------------------------------------------------
# 12. SEARCH atl ACROSS ALL CLASH EDGES
# ------------------------------------------------------------

atl_regex <- paste0(
  
  "SAA6008_RS05315",
  "|",
  "WP_001074557\\.1",
  "|",
  "\\batl\\b"
)


atl_all_edges <- edges |>
  filter(
    str_detect(
      all_edge_text,
      regex(
        atl_regex,
        ignore_case = TRUE
      )
    )
  )


# ------------------------------------------------------------
# 13. SEARCH HPr ACROSS ALL CLASH EDGES
# ------------------------------------------------------------

hpr_regex <- paste0(
  
  "SAA6008_RS05475",
  "|",
  "WP_000437472\\.1",
  "|",
  "\\bhpr\\b",
  "|",
  "phosphocarrier protein hpr"
)


hpr_all_edges <- edges |>
  filter(
    str_detect(
      all_edge_text,
      regex(
        hpr_regex,
        ignore_case = TRUE
      )
    )
  )


cat(
  "\n----------------------------------------\n",
  "KNOWN TARGET RECOVERY — ALL CLASH\n",
  "----------------------------------------\n"
)

cat(
  "atl-associated edges:",
  nrow(atl_all_edges),
  "\n"
)

cat(
  "HPr-associated edges:",
  nrow(hpr_all_edges),
  "\n"
)


# ------------------------------------------------------------
# 14. RsaOI -> atl DIRECT SEARCH
# ------------------------------------------------------------

rsaoi_atl <- rsaoi_edges |>
  filter(
    str_detect(
      all_edge_text,
      regex(
        atl_regex,
        ignore_case = TRUE
      )
    )
  )


# ------------------------------------------------------------
# 15. RsaOI -> HPr DIRECT SEARCH
# ------------------------------------------------------------

rsaoi_hpr <- rsaoi_edges |>
  filter(
    str_detect(
      all_edge_text,
      regex(
        hpr_regex,
        ignore_case = TRUE
      )
    )
  )


cat(
  "\n----------------------------------------\n",
  "DIRECT POSITIVE-CONTROL TEST\n",
  "----------------------------------------\n"
)

cat(
  "RsaOI -> atl:",
  nrow(rsaoi_atl),
  "\n"
)

cat(
  "RsaOI -> HPr:",
  nrow(rsaoi_hpr),
  "\n"
)


# ------------------------------------------------------------
# 16. IDENTIFY THE ACTUAL CURRENT RsaOI TARGET(S)
# ------------------------------------------------------------

rsaoi_target_columns <- intersect(
  
  c(
    "source_node",
    "source_label",
    "source_raw_names",
    
    "target_node",
    "target_label",
    "target_raw_names",
    
    "number_of_supporting_rows",
    "total_hybrid_count",
    
    "maximum_number_of_experiments",
    
    "best_raw_p_value",
    "best_adjusted_p_value",
    
    "best_connection_score",
    
    "target_log2_fold_change",
    "target_deseq_adjusted_p",
    
    "target_is_significant_deg",
    "target_strong_response"
  ),
  
  names(rsaoi_edges)
)


rsaoi_compact <- rsaoi_edges |>
  select(
    all_of(
      rsaoi_target_columns
    )
  )


readr::write_csv(
  rsaoi_compact,
  file.path(
    table_folder,
    "RsaOI_actual_CLASH_targets.csv"
  )
)


# ------------------------------------------------------------
# 17. VERIFY CURRENT RsaOI TARGET WP_000825929.1
#     AGAINST KNOWLEDGEBASE
# ------------------------------------------------------------

current_rsaoi_target_kb <- kb |>
  filter(
    str_detect(
      all_annotation_text,
      fixed(
        "wp_000825929.1"
      )
    )
  )


readr::write_csv(
  current_rsaoi_target_kb,
  file.path(
    table_folder,
    "RsaOI_WP_000825929_annotation.csv"
  )
)


cat(
  "\nWP_000825929.1 annotation rows:",
  nrow(current_rsaoi_target_kb),
  "\n"
)


# ------------------------------------------------------------
# 18. CHECK DEG STATUS SAFELY
# ------------------------------------------------------------

get_logical_column <- function(
    data,
    column_name
) {
  
  if (!column_name %in% names(data)) {
    
    return(
      rep(
        NA,
        nrow(data)
      )
    )
  }
  
  safe_logical(
    data[[column_name]]
  )
}


get_numeric_column <- function(
    data,
    column_name
) {
  
  if (!column_name %in% names(data)) {
    
    return(
      rep(
        NA_real_,
        nrow(data)
      )
    )
  }
  
  safe_num(
    data[[column_name]]
  )
}


rsaoi_edges$benchmark_target_significant <-
  get_logical_column(
    rsaoi_edges,
    "target_is_significant_deg"
  )


rsaoi_edges$benchmark_target_log2FC <-
  get_numeric_column(
    rsaoi_edges,
    "target_log2_fold_change"
  )


# ------------------------------------------------------------
# 19. CREATE KNOWN-TARGET BENCHMARK SUMMARY
# ------------------------------------------------------------

benchmark_summary <- tibble(
  
  benchmark = c(
    "RsaOI detected in collapsed CLASH",
    "atl detected anywhere in collapsed CLASH",
    "HPr detected anywhere in collapsed CLASH",
    "RsaOI-atl interaction detected",
    "RsaOI-HPr interaction detected"
  ),
  
  result = c(
    
    nrow(rsaoi_edges) > 0,
    
    nrow(atl_all_edges) > 0,
    
    nrow(hpr_all_edges) > 0,
    
    nrow(rsaoi_atl) > 0,
    
    nrow(rsaoi_hpr) > 0
  ),
  
  number_of_edges = c(
    
    nrow(rsaoi_edges),
    
    nrow(atl_all_edges),
    
    nrow(hpr_all_edges),
    
    nrow(rsaoi_atl),
    
    nrow(rsaoi_hpr)
  )
)


readr::write_csv(
  benchmark_summary,
  file.path(
    table_folder,
    "RsaOI_known_target_benchmark_summary.csv"
  )
)


print(
  benchmark_summary
)


# ------------------------------------------------------------
# 20. SAVE ALL TARGET SEARCH RESULTS
# ------------------------------------------------------------

readr::write_csv(
  atl_all_edges,
  file.path(
    table_folder,
    "atl_all_CLASH_edges.csv"
  )
)


readr::write_csv(
  hpr_all_edges,
  file.path(
    table_folder,
    "HPr_all_CLASH_edges.csv"
  )
)


readr::write_csv(
  rsaoi_atl,
  file.path(
    table_folder,
    "RsaOI_atl_direct_CLASH_edges.csv"
  )
)


readr::write_csv(
  rsaoi_hpr,
  file.path(
    table_folder,
    "RsaOI_HPr_direct_CLASH_edges.csv"
  )
)


# ------------------------------------------------------------
# 21. CLASSIFY BENCHMARK OUTCOME
# ------------------------------------------------------------

classification <- case_when(
  
  nrow(rsaoi_atl) > 0 |
    nrow(rsaoi_hpr) > 0 ~
    
    paste0(
      "KNOWN TARGET RECOVERED: ",
      "At least one RsaOI positive-control target ",
      "is present in the collapsed CLASH network."
    ),
  
  
  nrow(atl_all_edges) > 0 |
    nrow(hpr_all_edges) > 0 ~
    
    paste0(
      "KNOWN TARGET PRESENT ELSEWHERE: ",
      "atl and/or HPr is represented in the collapsed ",
      "CLASH network but was not recovered as an ",
      "RsaOI interaction."
    ),
  
  
  TRUE ~
    
    paste0(
      "KNOWN TARGETS NOT RECOVERED IN COLLAPSED TABLE: ",
      "RsaOI is present, but neither the RsaOI-atl nor ",
      "RsaOI-HPr positive-control interaction was found ",
      "using gene symbol, locus tag, or RefSeq protein ",
      "accession. This should be interpreted as a ",
      "dataset/pipeline recovery limitation, not evidence ",
      "that the experimentally reported interactions ",
      "do not exist."
    )
)


cat(
  "\n========================================\n",
  "BENCHMARK CLASSIFICATION\n",
  "========================================\n\n"
)

cat(
  classification,
  "\n"
)


# ------------------------------------------------------------
# 22. WRITE AUDIT REPORT
# ------------------------------------------------------------

report_file <- file.path(
  output_folder,
  "RsaOI_Known_Target_Audit_Report.txt"
)


sink(
  report_file
)


cat(
  "BacRegRNA Project\n"
)

cat(
  "RsaOI Known-Target Resolution Audit\n"
)

cat(
  "====================================\n\n"
)


cat(
  "Total collapsed CLASH edges: ",
  nrow(edges),
  "\n",
  sep = ""
)


cat(
  "RsaOI-associated edges: ",
  nrow(rsaoi_edges),
  "\n\n",
  sep = ""
)


cat(
  "atl benchmark\n",
  "-------------\n"
)

cat(
  "Locus tag: SAA6008_RS05315\n"
)

cat(
  "RefSeq protein: WP_001074557.1\n"
)

cat(
  "All CLASH edges containing atl identifiers: ",
  nrow(atl_all_edges),
  "\n",
  sep = ""
)

cat(
  "Direct RsaOI-atl edges: ",
  nrow(rsaoi_atl),
  "\n\n",
  sep = ""
)


cat(
  "HPr benchmark\n",
  "-------------\n"
)

cat(
  "Locus tag: SAA6008_RS05475\n"
)

cat(
  "RefSeq protein: WP_000437472.1\n"
)

cat(
  "All CLASH edges containing HPr identifiers: ",
  nrow(hpr_all_edges),
  "\n",
  sep = ""
)

cat(
  "Direct RsaOI-HPr edges: ",
  nrow(rsaoi_hpr),
  "\n\n",
  sep = ""
)


cat(
  "Current RsaOI target accession\n",
  "------------------------------\n"
)

cat(
  "WP_000825929.1 knowledgebase rows: ",
  nrow(current_rsaoi_target_kb),
  "\n\n",
  sep = ""
)


cat(
  "Final classification\n",
  "--------------------\n"
)

cat(
  classification,
  "\n\n"
)


cat(
  "Interpretation\n",
  "--------------\n"
)

cat(
  paste0(
    "This benchmark separates CLASH recovery from ",
    "transcriptomic filtering. Failure to recover a ",
    "known RsaOI-target pair in this collapsed table ",
    "must not be interpreted as biological absence. ",
    "It instead identifies a limitation in interaction ",
    "recovery, identifier mapping, dataset coverage, ",
    "or analytical selection."
  )
)


sink()


# ------------------------------------------------------------
# 23. SAVE SESSION INFORMATION
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_folder,
    "sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 24. FINAL MESSAGE
# ------------------------------------------------------------

cat(
  "\n========================================\n"
)

cat(
  "SCRIPT 13B COMPLETED SUCCESSFULLY\n"
)

cat(
  "========================================\n"
)

cat(
  "\nOutput folder:\n",
  output_folder,
  "\n"
)

cat(
  "\nMOST IMPORTANT FILE:\n"
)

cat(
  "tables/RsaOI_known_target_benchmark_summary.csv\n"
)

cat(
  "\nAlso inspect:\n"
)

cat(
  "tables/RsaOI_actual_CLASH_targets.csv\n"
)

cat(
  "tables/atl_all_CLASH_edges.csv\n"
)

cat(
  "tables/HPr_all_CLASH_edges.csv\n"
)

cat(
  "tables/RsaOI_atl_direct_CLASH_edges.csv\n"
)

cat(
  "tables/RsaOI_HPr_direct_CLASH_edges.csv\n"
)

cat(
  "\nDO NOT MODIFY SCRIPT 14 YET.\n"
)