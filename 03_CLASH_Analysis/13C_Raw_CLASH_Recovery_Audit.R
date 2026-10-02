# ============================================================
# BacRegRNA Project
# Script 13C — RAW CLASH RECOVERY AUDIT
#
# PURPOSE
# ------------------------------------------------------------
# Trace selected benchmark/candidate RNAs through:
#
# RAW GSE254532 CLASH
#       ↓
# Script 13 standardized interactions
#       ↓
# sRNA-mRNA all rows
#       ↓
# collapsed sRNA-mRNA edges
#       ↓
# vancomycin-responsive target subset
#
# Benchmark entities:
#   RsaOI
#   SprA2 / RsaJ
#   atl
#   HPr
#   nrdF
#   3UTR-RS04120
#   3UTR-RS05585
#
# IMPORTANT:
# This is an AUDIT ONLY.
# It does not modify any previous result.
# ============================================================


# ------------------------------------------------------------
# 1. PROJECT PATHS
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

raw_file <- file.path(
  project_folder,
  "07_processed_data",
  "downloaded_files",
  "GSE254532",
  "GSE254532_Final_RATT_VISACLASH.txt.gz"
)

script13_folder <- file.path(
  project_folder,
  "13_CLASH_network",
  "tables"
)

standardized_file <- file.path(
  script13_folder,
  "GSE254532_CLASH_standardized_all_interactions.csv"
)

all_srna_mrna_file <- file.path(
  script13_folder,
  "GSE254532_sRNA_mRNA_edges_all_rows.csv"
)

collapsed_file <- file.path(
  script13_folder,
  "GSE254532_sRNA_mRNA_edges_collapsed.csv"
)

responsive_file <- file.path(
  script13_folder,
  "GSE254532_vancomycin_responsive_CLASH_targets.csv"
)


output_folder <- file.path(
  project_folder,
  "13C_raw_CLASH_recovery_audit"
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
# 3. VERIFY FILES
# ------------------------------------------------------------

required_files <- c(
  raw_file,
  standardized_file,
  all_srna_mrna_file,
  collapsed_file,
  responsive_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat(
    "\nMissing required files:\n"
  )
  
  print(missing_files)
  
  stop(
    "Run Script 13 first and verify project paths."
  )
}


# ------------------------------------------------------------
# 4. READ FILES
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "SCRIPT 13C — RAW CLASH RECOVERY AUDIT\n",
  "========================================\n"
)


raw <- readr::read_tsv(
  raw_file,
  show_col_types = FALSE,
  progress = FALSE
)


standardized <- readr::read_csv(
  standardized_file,
  show_col_types = FALSE
)


all_srna_mrna <- readr::read_csv(
  all_srna_mrna_file,
  show_col_types = FALSE
)


collapsed <- readr::read_csv(
  collapsed_file,
  show_col_types = FALSE
)


responsive <- readr::read_csv(
  responsive_file,
  show_col_types = FALSE
)


cat(
  "\nDataset sizes:\n"
)

cat(
  "RAW CLASH rows: ",
  nrow(raw),
  "\n",
  sep = ""
)

cat(
  "Standardized rows: ",
  nrow(standardized),
  "\n",
  sep = ""
)

cat(
  "sRNA-mRNA all rows: ",
  nrow(all_srna_mrna),
  "\n",
  sep = ""
)

cat(
  "Collapsed sRNA-mRNA edges: ",
  nrow(collapsed),
  "\n",
  sep = ""
)

cat(
  "Responsive-target edges: ",
  nrow(responsive),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. HELPER FUNCTIONS
# ------------------------------------------------------------

safe_chr <- function(x) {
  
  x <- as.character(x)
  
  x[is.na(x)] <- ""
  
  x
}


build_search_text <- function(data) {
  
  apply(
    data,
    1,
    function(row) {
      
      paste(
        safe_chr(row),
        collapse = " | "
      )
    }
  ) |>
    stringr::str_to_lower()
}


# ------------------------------------------------------------
# 6. BUILD SEARCHABLE REPRESENTATION
# ------------------------------------------------------------

raw$AUDIT_SEARCH_TEXT <-
  build_search_text(raw)

standardized$AUDIT_SEARCH_TEXT <-
  build_search_text(standardized)

all_srna_mrna$AUDIT_SEARCH_TEXT <-
  build_search_text(all_srna_mrna)

collapsed$AUDIT_SEARCH_TEXT <-
  build_search_text(collapsed)

responsive$AUDIT_SEARCH_TEXT <-
  build_search_text(responsive)


# ------------------------------------------------------------
# 7. DEFINE BENCHMARK ENTITIES
#
# Patterns deliberately include names + known accessions.
# ------------------------------------------------------------

benchmark_dictionary <- tibble::tribble(
  
  ~entity,
  ~pattern,
  
  "RsaOI",
  "rsaoi|sau[-_ ]?6477",
  
  "SprA2_RsaJ",
  "spra2|rsaj",
  
  "atl",
  "wp_001074557\\.1|saa6008_rs05315|\\batl\\b",
  
  "HPr",
  paste0(
    "wp_000437472\\.1",
    "|saa6008_rs05475",
    "|phosphocarrier protein hpr",
    "|\\bhpr\\b"
  ),
  
  "nrdF",
  "\\bnrdf\\b",
  
  "3UTR_RS04120",
  "rs04120|3utr[-_].*04120",
  
  "3UTR_RS05585",
  "rs05585|3utr[-_].*05585"
)


# ------------------------------------------------------------
# 8. FUNCTION: COUNT ENTITY AT EACH PIPELINE STAGE
# ------------------------------------------------------------

count_entity <- function(
    data,
    pattern
) {
  
  sum(
    str_detect(
      data$AUDIT_SEARCH_TEXT,
      regex(
        pattern,
        ignore_case = TRUE
      )
    ),
    na.rm = TRUE
  )
}


# ------------------------------------------------------------
# 9. PIPELINE RECOVERY TABLE
# ------------------------------------------------------------

recovery_table <- benchmark_dictionary |>
  rowwise() |>
  mutate(
    
    raw_rows =
      count_entity(
        raw,
        pattern
      ),
    
    standardized_rows =
      count_entity(
        standardized,
        pattern
      ),
    
    srna_mrna_rows =
      count_entity(
        all_srna_mrna,
        pattern
      ),
    
    collapsed_edges =
      count_entity(
        collapsed,
        pattern
      ),
    
    responsive_edges =
      count_entity(
        responsive,
        pattern
      )
  ) |>
  ungroup()


# ------------------------------------------------------------
# 10. DERIVE RECOVERY FLAGS
# ------------------------------------------------------------

recovery_table <- recovery_table |>
  mutate(
    
    detected_raw =
      raw_rows > 0,
    
    survives_standardization =
      standardized_rows > 0,
    
    classified_srna_mrna =
      srna_mrna_rows > 0,
    
    survives_collapse =
      collapsed_edges > 0,
    
    survives_DE_focus =
      responsive_edges > 0
  )


# ------------------------------------------------------------
# 11. DETERMINE FIRST APPARENT LOSS STAGE
# ------------------------------------------------------------

recovery_table <- recovery_table |>
  mutate(
    
    first_apparent_loss =
      case_when(
        
        raw_rows == 0 ~
          "Not detected in raw CLASH",
        
        standardized_rows == 0 ~
          "Lost during standardization/mapping",
        
        srna_mrna_rows == 0 ~
          "Excluded during sRNA-mRNA classification",
        
        collapsed_edges == 0 ~
          "Lost during edge collapsing",
        
        responsive_edges == 0 ~
          "Excluded by vancomycin-responsive target filter",
        
        TRUE ~
          "Retained through all audited stages"
      )
  )


# ------------------------------------------------------------
# 12. SAVE RECOVERY TABLE
# ------------------------------------------------------------

write_csv(
  recovery_table,
  file.path(
    table_folder,
    "CLASH_pipeline_recovery_summary.csv"
  )
)


cat(
  "\n----------------------------------------\n",
  "PIPELINE RECOVERY SUMMARY\n",
  "----------------------------------------\n"
)

print(
  recovery_table
)


# ------------------------------------------------------------
# 13. EXTRACT ENTITY-SPECIFIC ROWS
# ------------------------------------------------------------

extract_entity <- function(
    data,
    pattern
) {
  
  data |>
    filter(
      str_detect(
        AUDIT_SEARCH_TEXT,
        regex(
          pattern,
          ignore_case = TRUE
        )
      )
    ) |>
    select(
      -AUDIT_SEARCH_TEXT
    )
}


for (i in seq_len(nrow(benchmark_dictionary))) {
  
  entity_name <-
    benchmark_dictionary$entity[i]
  
  entity_pattern <-
    benchmark_dictionary$pattern[i]
  
  
  raw_hits <- extract_entity(
    raw,
    entity_pattern
  )
  
  
  standardized_hits <- extract_entity(
    standardized,
    entity_pattern
  )
  
  
  srna_mrna_hits <- extract_entity(
    all_srna_mrna,
    entity_pattern
  )
  
  
  collapsed_hits <- extract_entity(
    collapsed,
    entity_pattern
  )
  
  
  responsive_hits <- extract_entity(
    responsive,
    entity_pattern
  )
  
  
  write_csv(
    raw_hits,
    file.path(
      table_folder,
      paste0(
        entity_name,
        "_01_RAW.csv"
      )
    )
  )
  
  
  write_csv(
    standardized_hits,
    file.path(
      table_folder,
      paste0(
        entity_name,
        "_02_STANDARDIZED.csv"
      )
    )
  )
  
  
  write_csv(
    srna_mrna_hits,
    file.path(
      table_folder,
      paste0(
        entity_name,
        "_03_sRNA_mRNA.csv"
      )
    )
  )
  
  
  write_csv(
    collapsed_hits,
    file.path(
      table_folder,
      paste0(
        entity_name,
        "_04_COLLAPSED.csv"
      )
    )
  )
  
  
  write_csv(
    responsive_hits,
    file.path(
      table_folder,
      paste0(
        entity_name,
        "_05_RESPONSIVE.csv"
      )
    )
  )
}


# ------------------------------------------------------------
# 14. PAIR-LEVEL SEARCH FUNCTION
# ------------------------------------------------------------

find_pair <- function(
    data,
    source_pattern,
    target_pattern
) {
  
  data |>
    filter(
      
      str_detect(
        AUDIT_SEARCH_TEXT,
        regex(
          source_pattern,
          ignore_case = TRUE
        )
      ),
      
      str_detect(
        AUDIT_SEARCH_TEXT,
        regex(
          target_pattern,
          ignore_case = TRUE
        )
      )
    ) |>
    select(
      -AUDIT_SEARCH_TEXT
    )
}


# ------------------------------------------------------------
# 15. DEFINE CRITICAL PAIRS
# ------------------------------------------------------------

rsaoi_pattern <-
  "rsaoi|sau[-_ ]?6477"

atl_pattern <-
  "wp_001074557\\.1|saa6008_rs05315|\\batl\\b"

hpr_pattern <-
  paste0(
    "wp_000437472\\.1",
    "|saa6008_rs05475",
    "|phosphocarrier protein hpr",
    "|\\bhpr\\b"
  )

nrdf_pattern <-
  "\\bnrdf\\b"

utr04120_pattern <-
  "rs04120|3utr[-_].*04120"

utr05585_pattern <-
  "rs05585|3utr[-_].*05585"


# ------------------------------------------------------------
# 16. SEARCH CRITICAL PAIRS IN RAW DATA
# ------------------------------------------------------------

raw_rsaoi_atl <- find_pair(
  raw,
  rsaoi_pattern,
  atl_pattern
)

raw_rsaoi_hpr <- find_pair(
  raw,
  rsaoi_pattern,
  hpr_pattern
)

raw_04120_nrdf <- find_pair(
  raw,
  utr04120_pattern,
  nrdf_pattern
)

raw_05585_nrdf <- find_pair(
  raw,
  utr05585_pattern,
  nrdf_pattern
)


# ------------------------------------------------------------
# 17. SEARCH CRITICAL PAIRS AFTER STANDARDIZATION
# ------------------------------------------------------------

std_rsaoi_atl <- find_pair(
  standardized,
  rsaoi_pattern,
  atl_pattern
)

std_rsaoi_hpr <- find_pair(
  standardized,
  rsaoi_pattern,
  hpr_pattern
)

std_04120_nrdf <- find_pair(
  standardized,
  utr04120_pattern,
  nrdf_pattern
)

std_05585_nrdf <- find_pair(
  standardized,
  utr05585_pattern,
  nrdf_pattern
)


# ------------------------------------------------------------
# 18. SEARCH CRITICAL PAIRS IN COLLAPSED NETWORK
# ------------------------------------------------------------

collapsed_rsaoi_atl <- find_pair(
  collapsed,
  rsaoi_pattern,
  atl_pattern
)

collapsed_rsaoi_hpr <- find_pair(
  collapsed,
  rsaoi_pattern,
  hpr_pattern
)

collapsed_04120_nrdf <- find_pair(
  collapsed,
  utr04120_pattern,
  nrdf_pattern
)

collapsed_05585_nrdf <- find_pair(
  collapsed,
  utr05585_pattern,
  nrdf_pattern
)


# ------------------------------------------------------------
# 19. CRITICAL PAIR SUMMARY
# ------------------------------------------------------------

pair_summary <- tibble(
  
  interaction = c(
    "RsaOI -> atl",
    "RsaOI -> HPr",
    "3UTR-RS04120 -> nrdF",
    "3UTR-RS05585 -> nrdF"
  ),
  
  raw_rows = c(
    nrow(raw_rsaoi_atl),
    nrow(raw_rsaoi_hpr),
    nrow(raw_04120_nrdf),
    nrow(raw_05585_nrdf)
  ),
  
  standardized_rows = c(
    nrow(std_rsaoi_atl),
    nrow(std_rsaoi_hpr),
    nrow(std_04120_nrdf),
    nrow(std_05585_nrdf)
  ),
  
  collapsed_edges = c(
    nrow(collapsed_rsaoi_atl),
    nrow(collapsed_rsaoi_hpr),
    nrow(collapsed_04120_nrdf),
    nrow(collapsed_05585_nrdf)
  )
)


pair_summary <- pair_summary |>
  mutate(
    
    interpretation =
      case_when(
        
        raw_rows == 0 ~
          paste0(
            "Pair not recovered in the raw CLASH table; ",
            "do not interpret as biological absence."
          ),
        
        raw_rows > 0 &
          standardized_rows == 0 ~
          "Present in raw CLASH but lost during standardization.",
        
        standardized_rows > 0 &
          collapsed_edges == 0 ~
          "Present after standardization but absent after network collapse/classification.",
        
        collapsed_edges > 0 ~
          "Recovered in collapsed CLASH network.",
        
        TRUE ~
          "Requires manual review."
      )
  )


write_csv(
  pair_summary,
  file.path(
    table_folder,
    "Critical_pair_recovery_summary.csv"
  )
)


cat(
  "\n----------------------------------------\n",
  "CRITICAL PAIR RECOVERY\n",
  "----------------------------------------\n"
)

print(
  pair_summary
)


# ------------------------------------------------------------
# 20. SAVE RAW CRITICAL PAIRS
# ------------------------------------------------------------

write_csv(
  raw_rsaoi_atl,
  file.path(
    table_folder,
    "RAW_RsaOI_atl.csv"
  )
)

write_csv(
  raw_rsaoi_hpr,
  file.path(
    table_folder,
    "RAW_RsaOI_HPr.csv"
  )
)

write_csv(
  raw_04120_nrdf,
  file.path(
    table_folder,
    "RAW_3UTR_RS04120_nrdF.csv"
  )
)

write_csv(
  raw_05585_nrdf,
  file.path(
    table_folder,
    "RAW_3UTR_RS05585_nrdF.csv"
  )
)


# ------------------------------------------------------------
# 21. SPECIAL AUDIT — ALL RAW RsaOI INTERACTIONS
# ------------------------------------------------------------

raw_rsaoi_all <- extract_entity(
  raw,
  rsaoi_pattern
)

write_csv(
  raw_rsaoi_all,
  file.path(
    table_folder,
    "RAW_RsaOI_ALL_interactions.csv"
  )
)


# ------------------------------------------------------------
# 22. SPECIAL AUDIT — ALL RAW atl / HPr
# ------------------------------------------------------------

raw_atl_all <- extract_entity(
  raw,
  atl_pattern
)

raw_hpr_all <- extract_entity(
  raw,
  hpr_pattern
)

write_csv(
  raw_atl_all,
  file.path(
    table_folder,
    "RAW_atl_ALL_interactions.csv"
  )
)

write_csv(
  raw_hpr_all,
  file.path(
    table_folder,
    "RAW_HPr_ALL_interactions.csv"
  )
)


# ------------------------------------------------------------
# 23. CREATE LONG-FORM FLOW TABLE
# ------------------------------------------------------------

flow_table <- recovery_table |>
  select(
    entity,
    raw_rows,
    standardized_rows,
    srna_mrna_rows,
    collapsed_edges,
    responsive_edges
  ) |>
  pivot_longer(
    
    cols = -entity,
    
    names_to = "pipeline_stage",
    
    values_to = "number_detected"
  )


write_csv(
  flow_table,
  file.path(
    table_folder,
    "CLASH_pipeline_flow_long.csv"
  )
)


# ------------------------------------------------------------
# 24. HUMAN-READABLE AUDIT REPORT
# ------------------------------------------------------------

report_file <- file.path(
  output_folder,
  "Raw_CLASH_Recovery_Audit_Report.txt"
)


sink(
  report_file
)


cat(
  "BacRegRNA Project\n"
)

cat(
  "Raw CLASH Recovery Audit\n"
)

cat(
  "========================================\n\n"
)


cat(
  "PIPELINE SIZE\n"
)

cat(
  "-------------\n"
)

cat(
  "Raw CLASH rows: ",
  nrow(raw),
  "\n",
  sep = ""
)

cat(
  "Standardized rows: ",
  nrow(standardized),
  "\n",
  sep = ""
)

cat(
  "sRNA-mRNA rows: ",
  nrow(all_srna_mrna),
  "\n",
  sep = ""
)

cat(
  "Collapsed edges: ",
  nrow(collapsed),
  "\n",
  sep = ""
)

cat(
  "Responsive-target edges: ",
  nrow(responsive),
  "\n\n",
  sep = ""
)


cat(
  "ENTITY RECOVERY\n"
)

cat(
  "---------------\n"
)

print(
  recovery_table
)


cat(
  "\n\nCRITICAL PAIR RECOVERY\n"
)

cat(
  "----------------------\n"
)

print(
  pair_summary
)


cat(
  "\n\nINTERPRETATION PRINCIPLE\n"
)

cat(
  "------------------------\n"
)

cat(
  paste0(
    "This audit distinguishes failure to recover an ",
    "interaction computationally from evidence against ",
    "the interaction biologically. A pair absent from ",
    "the analyzed raw CLASH table cannot be created by ",
    "changing downstream prioritization thresholds. ",
    "Likewise, an interaction present in CLASH but ",
    "excluded by transcriptomic filtering remains ",
    "CLASH-supported and should be represented separately ",
    "from vancomycin-responsive target evidence."
  )
)


sink()


# ------------------------------------------------------------
# 25. SESSION INFORMATION
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_folder,
    "sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 26. FINAL CONSOLE OUTPUT
# ------------------------------------------------------------

cat(
  "\n========================================\n"
)

cat(
  "SCRIPT 13C COMPLETED SUCCESSFULLY\n"
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
  "\nPLEASE INSPECT THESE TWO FILES FIRST:\n"
)

cat(
  "1. tables/CLASH_pipeline_recovery_summary.csv\n"
)

cat(
  "2. tables/Critical_pair_recovery_summary.csv\n"
)

cat(
  "\nAdditional audit files:\n"
)

cat(
  "3. tables/RAW_RsaOI_ALL_interactions.csv\n"
)

cat(
  "4. tables/RAW_atl_ALL_interactions.csv\n"
)

cat(
  "5. tables/RAW_HPr_ALL_interactions.csv\n"
)

cat(
  "6. Raw_CLASH_Recovery_Audit_Report.txt\n"
)

cat(
  "\nDO NOT CHANGE SCRIPT 14 YET.\n"
)