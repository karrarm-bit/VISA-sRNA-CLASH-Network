# ============================================================
# BacRegRNA Project
# Script 06B: Apply confirmed manual metadata corrections
# Fixed version: replicate_curated kept numeric
# ============================================================

# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

input_file <- file.path(
  project_folder,
  "06_curated_design",
  "BacRegRNA_curated_sample_design.csv"
)

output_folder <- file.path(
  project_folder,
  "06_curated_design"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Install and load packages
# ------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "stringr"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(readr)
library(dplyr)
library(stringr)

# ------------------------------------------------------------
# 3. Check and read input file
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop(
    "Input file was not found:\n",
    input_file
  )
}

design <- read_csv(
  input_file,
  show_col_types = FALSE
)

# Ensure replicate column is numeric before case_when()
design <- design |>
  mutate(
    replicate_curated = suppressWarnings(
      as.numeric(replicate_curated)
    )
  )

# ------------------------------------------------------------
# 4. Apply confirmed corrections
# ------------------------------------------------------------

final_design <- design |>
  mutate(
    title_lower = str_to_lower(
      coalesce(sample_title, "")
    ),
    
    # --------------------------------------------------------
    # Treatment corrections
    # --------------------------------------------------------
    
    treatment_curated = case_when(
      
      # GSE254530:
      # Typographical error "vancomcycin" means vancomycin
      dataset == "GSE254530" &
        str_detect(
          title_lower,
          "vancomcycin_30mins_3"
        ) ~ "Vancomycin",
      
      # GSE254531:
      # Both Ribo-seq samples are untreated baseline controls
      dataset == "GSE254531" &
        str_detect(
          title_lower,
          "ribo_seq_1|ribo_seq_2|riboseq_cont_1|riboseq_cont_2"
        ) ~ "Untreated",
      
      # GSE254532:
      # Explicit no-treatment CLASH samples
      dataset == "GSE254532" &
        str_detect(
          title_lower,
          "rnaseiii-clash_no_treatment_1|rnaseiii-clash_no_treatment_2"
        ) ~ "Untreated",
      
      TRUE ~ treatment_curated
    ),
    
    # --------------------------------------------------------
    # Vancomycin status
    # --------------------------------------------------------
    
    vancomycin_status = case_when(
      
      treatment_curated == "Vancomycin" ~ "Treated",
      
      treatment_curated %in% c(
        "Untreated",
        "No antibiotic"
      ) ~ "Untreated",
      
      treatment_curated == "Background control" ~
        "Technical control",
      
      TRUE ~ vancomycin_status
    ),
    
    # --------------------------------------------------------
    # Concentration corrections
    # --------------------------------------------------------
    
    concentration_curated = case_when(
      
      dataset == "GSE254530" &
        str_detect(
          title_lower,
          "vancomcycin_30mins_3"
        ) ~ "8 ug/mL",
      
      TRUE ~ concentration_curated
    ),
    
    # --------------------------------------------------------
    # Exposure-time corrections
    # --------------------------------------------------------
    
    exposure_time_curated = case_when(
      
      dataset == "GSE254530" &
        str_detect(
          title_lower,
          "vancomcycin_30mins_3"
        ) ~ "30 min",
      
      TRUE ~ exposure_time_curated
    ),
    
    # --------------------------------------------------------
    # Replicate corrections
    # Values are numeric, not character strings
    # --------------------------------------------------------
    
    replicate_curated = case_when(
      
      dataset == "GSE254530" &
        str_detect(
          title_lower,
          "vancomcycin_30mins_3"
        ) ~ 3,
      
      dataset == "GSE254531" &
        str_detect(
          title_lower,
          "ribo_seq_1|riboseq_cont_1"
        ) ~ 1,
      
      dataset == "GSE254531" &
        str_detect(
          title_lower,
          "ribo_seq_2|riboseq_cont_2"
        ) ~ 2,
      
      dataset == "GSE254532" &
        str_detect(
          title_lower,
          "rnaseiii-clash_no_treatment_1|no_treatment_1"
        ) ~ 1,
      
      dataset == "GSE254532" &
        str_detect(
          title_lower,
          "rnaseiii-clash_no_treatment_2|no_treatment_2"
        ) ~ 2,
      
      TRUE ~ replicate_curated
    ),
    
    # --------------------------------------------------------
    # Analysis-module corrections
    # --------------------------------------------------------
    
    analysis_module = case_when(
      
      dataset == "GSE254530" &
        assay_curated == "RNA-seq" ~
        "Discovery_RNAseq_vancomycin_response",
      
      dataset == "GSE254531" &
        assay_curated == "Ribo-seq" ~
        "Discovery_Riboseq_baseline_only",
      
      dataset == "GSE254532" &
        assay_curated == "RNase III-CLASH" &
        vancomycin_status %in% c(
          "Treated",
          "Untreated"
        ) ~
        "Discovery_CLASH_vancomycin_response",
      
      dataset == "GSE254532" &
        vancomycin_status == "Technical control" ~
        "CLASH_technical_control",
      
      TRUE ~ analysis_module
    ),
    
    # --------------------------------------------------------
    # Record the reason for every manual correction
    # --------------------------------------------------------
    
    correction_note = case_when(
      
      dataset == "GSE254530" &
        str_detect(
          title_lower,
          "vancomcycin_30mins_3"
        ) ~ paste0(
          "Corrected typographical error: ",
          "vancomycin-treated replicate 3, ",
          "8 ug/mL for 30 min"
        ),
      
      dataset == "GSE254531" &
        str_detect(
          title_lower,
          "ribo_seq_1|ribo_seq_2|riboseq_cont_1|riboseq_cont_2"
        ) ~ paste0(
          "Confirmed as untreated baseline Ribo-seq control; ",
          "no vancomycin-treated Ribo-seq counterpart identified"
        ),
      
      dataset == "GSE254532" &
        str_detect(
          title_lower,
          "rnaseiii-clash_no_treatment_1|rnaseiii-clash_no_treatment_2"
        ) ~
        "Confirmed as untreated RNase III-CLASH sample",
      
      TRUE ~ NA_character_
    )
  ) |>
  select(-title_lower) |>
  arrange(
    dataset,
    analysis_module,
    vancomycin_status,
    replicate_curated
  )

# ------------------------------------------------------------
# 5. Build final summary
# ------------------------------------------------------------

final_summary <- final_design |>
  count(
    dataset,
    analysis_module,
    assay_curated,
    genotype_curated,
    vancomycin_status,
    concentration_curated,
    exposure_time_curated,
    name = "number_of_samples"
  ) |>
  arrange(
    dataset,
    analysis_module,
    genotype_curated,
    vancomycin_status
  )

# ------------------------------------------------------------
# 6. Identify remaining unclear samples
# ------------------------------------------------------------

remaining_unclear <- final_design |>
  filter(
    vancomycin_status == "Unclear" |
      analysis_module == "Exclude_or_review"
  )

# ------------------------------------------------------------
# 7. Create corrections log
# ------------------------------------------------------------

corrected_samples <- final_design |>
  filter(
    !is.na(correction_note),
    correction_note != ""
  ) |>
  select(
    dataset,
    sample_id,
    gsm_accession,
    sample_title,
    assay_curated,
    genotype_curated,
    treatment_curated,
    vancomycin_status,
    concentration_curated,
    exposure_time_curated,
    replicate_curated,
    analysis_module,
    correction_note
  )

# ------------------------------------------------------------
# 8. Basic validation checks
# ------------------------------------------------------------

validation_checks <- tibble(
  check = c(
    "GSE254530 untreated RNA-seq samples",
    "GSE254530 treated RNA-seq samples",
    "GSE254531 baseline Ribo-seq samples",
    "GSE254532 untreated CLASH samples",
    "GSE254532 treated CLASH samples",
    "Remaining unclear samples"
  ),
  
  observed = c(
    
    sum(
      final_design$dataset == "GSE254530" &
        final_design$assay_curated == "RNA-seq" &
        final_design$vancomycin_status == "Untreated",
      na.rm = TRUE
    ),
    
    sum(
      final_design$dataset == "GSE254530" &
        final_design$assay_curated == "RNA-seq" &
        final_design$vancomycin_status == "Treated",
      na.rm = TRUE
    ),
    
    sum(
      final_design$dataset == "GSE254531" &
        final_design$assay_curated == "Ribo-seq" &
        final_design$vancomycin_status == "Untreated",
      na.rm = TRUE
    ),
    
    sum(
      final_design$dataset == "GSE254532" &
        final_design$assay_curated == "RNase III-CLASH" &
        final_design$vancomycin_status == "Untreated",
      na.rm = TRUE
    ),
    
    sum(
      final_design$dataset == "GSE254532" &
        final_design$assay_curated == "RNase III-CLASH" &
        final_design$vancomycin_status == "Treated",
      na.rm = TRUE
    ),
    
    nrow(remaining_unclear)
  ),
  
  expected = c(
    3,
    3,
    2,
    2,
    4,
    0
  )
) |>
  mutate(
    status = if_else(
      observed == expected,
      "PASS",
      "REVIEW"
    )
  )

# ------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------

write_csv(
  final_design,
  file.path(
    output_folder,
    "BacRegRNA_final_curated_sample_design.csv"
  )
)

write_csv(
  final_summary,
  file.path(
    output_folder,
    "BacRegRNA_final_curated_design_summary.csv"
  )
)

write_csv(
  corrected_samples,
  file.path(
    output_folder,
    "BacRegRNA_manual_corrections_log.csv"
  )
)

write_csv(
  remaining_unclear,
  file.path(
    output_folder,
    "BacRegRNA_remaining_unclear_samples.csv"
  )
)

write_csv(
  validation_checks,
  file.path(
    output_folder,
    "BacRegRNA_design_validation_checks.csv"
  )
)

# ------------------------------------------------------------
# 10. Final console report
# ------------------------------------------------------------

cat("\n============================================\n")
cat("MANUAL CORRECTIONS APPLIED SUCCESSFULLY\n")
cat("============================================\n")

cat(
  "Total samples:",
  nrow(final_design),
  "\n"
)

cat(
  "Corrected samples:",
  nrow(corrected_samples),
  "\n"
)

cat(
  "Remaining unclear samples:",
  nrow(remaining_unclear),
  "\n\n"
)

cat("Validation checks:\n")

print(
  validation_checks,
  n = Inf
)

cat("\nFinal design summary:\n")

print(
  final_summary,
  n = Inf
)

cat(
  "\nFiles saved inside:\n",
  output_folder,
  "\n"
)