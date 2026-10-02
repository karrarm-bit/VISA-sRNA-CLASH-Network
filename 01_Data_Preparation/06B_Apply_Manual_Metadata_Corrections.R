# ============================================================
# BacRegRNA Project
# Script 06B: Apply manual metadata corrections
# ============================================================

library(readr)
library(dplyr)
library(stringr)

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

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    input_file
  )
}

design <- read_csv(
  input_file,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Apply confirmed corrections
# ------------------------------------------------------------

final_design <- design |>
  mutate(
    title_lower = str_to_lower(sample_title),
    
    # --------------------------------------------------------
    # 1. GSE254530:
    # vancomcycin_30mins_3 is treated replicate 3
    # --------------------------------------------------------
    
    treatment_curated = case_when(
      
      dataset == "GSE254530" &
        str_detect(
          title_lower,
          "vancomcycin_30mins_3"
        ) ~ "Vancomycin",
      
      # ------------------------------------------------------
      # 2. GSE254531:
      # both Ribo-seq samples are untreated controls
      # ------------------------------------------------------
      
      dataset == "GSE254531" &
        str_detect(
          title_lower,
          "ribo_seq_1|ribo_seq_2|riboseq_cont"
        ) ~ "Untreated",
      
      # ------------------------------------------------------
      # 3. GSE254532:
      # explicit no-treatment CLASH samples
      # ------------------------------------------------------
      
      dataset == "GSE254532" &
        str_detect(
          title_lower,
          "rnaseiii-clash_no_treatment_1|rnaseiii-clash_no_treatment_2"
        ) ~ "Untreated",
      
      TRUE ~ treatment_curated
    ),
    
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
    
    concentration_curated = case_when(
      
      dataset == "GSE254530" &
        str_detect(
          title_lower,
          "vancomcycin_30mins_3"
        ) ~ "8 ug/mL",
      
      TRUE ~ concentration_curated
    ),
    
    exposure_time_curated = case_when(
      
      dataset == "GSE254530" &
        str_detect(
          title_lower,
          "vancomcycin_30mins_3"
        ) ~ "30 min",
      
      TRUE ~ exposure_time_curated
    ),
    
    replicate_curated = case_when(
      
      dataset == "GSE254530" &
        str_detect(
          title_lower,
          "vancomcycin_30mins_3"
        ) ~ "3",
      
      dataset == "GSE254531" &
        str_detect(
          title_lower,
          "ribo_seq_1|riboseq_cont_1"
        ) ~ "1",
      
      dataset == "GSE254531" &
        str_detect(
          title_lower,
          "ribo_seq_2|riboseq_cont_2"
        ) ~ "2",
      
      dataset == "GSE254532" &
        str_detect(
          title_lower,
          "no_treatment_1"
        ) ~ "1",
      
      dataset == "GSE254532" &
        str_detect(
          title_lower,
          "no_treatment_2"
        ) ~ "2",
      
      TRUE ~ replicate_curated
    ),
    
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
    
    correction_note = case_when(
      
      dataset == "GSE254530" &
        str_detect(
          title_lower,
          "vancomcycin_30mins_3"
        ) ~
        "Corrected typo: treated replicate 3, 8 ug/mL, 30 min",
      
      dataset == "GSE254531" &
        str_detect(
          title_lower,
          "ribo_seq_1|ribo_seq_2|riboseq_cont"
        ) ~
        "Confirmed as untreated baseline Ribo-seq control",
      
      dataset == "GSE254532" &
        str_detect(
          title_lower,
          "rnaseiii-clash_no_treatment_1|rnaseiii-clash_no_treatment_2"
        ) ~
        "Confirmed as untreated CLASH sample",
      
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
# Final validation summary
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

remaining_unclear <- final_design |>
  filter(
    vancomycin_status == "Unclear" |
      analysis_module == "Exclude_or_review"
  )

corrected_samples <- final_design |>
  filter(!is.na(correction_note)) |>
  select(
    dataset,
    sample_id,
    gsm_accession,
    sample_title,
    treatment_curated,
    vancomycin_status,
    concentration_curated,
    exposure_time_curated,
    replicate_curated,
    analysis_module,
    correction_note
  )

# ------------------------------------------------------------
# Save outputs
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

# ------------------------------------------------------------
# Final report
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

print(
  final_summary,
  n = Inf
)

cat(
  "\nFiles saved inside:\n",
  output_folder,
  "\n"
)