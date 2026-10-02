# ============================================================
# BacRegRNA Project
# Script 06: Curate the experimental design
# ============================================================

library(readr)
library(dplyr)
library(stringr)
library(tidyr)

project_folder <- "D:/Bac-sRNA"

input_file <- file.path(
  project_folder,
  "05_parsed_metadata",
  "BacRegRNA_all_sample_metadata.csv"
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

metadata <- read_csv(
  input_file,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 1. Normalize sample titles
# ------------------------------------------------------------

curated <- metadata |>
  mutate(
    title_lower = str_to_lower(sample_title),
    
    # Extract SRX directly from the relations field
    srx_accession = str_extract(
      relations,
      "SRX[0-9]+"
    ),
    
    biosample_accession = str_extract(
      relations,
      "SAMN[0-9]+"
    )
  )

# ------------------------------------------------------------
# 2. Infer the actual assay from each sample
# ------------------------------------------------------------

curated <- curated |>
  mutate(
    assay_curated = case_when(
      
      str_detect(
        title_lower,
        "rnase iii-clash|rnaseiii-clash|rnc-htf"
      ) ~ "RNase III-CLASH",
      
      str_detect(
        title_lower,
        "ribo_seq|ribo-seq|riboseq"
      ) ~ "Ribo-seq",
      
      str_detect(
        title_lower,
        "drna-seq"
      ) ~ "dRNA-seq",
      
      str_detect(
        title_lower,
        "term-seq"
      ) ~ "Term-seq",
      
      str_detect(
        title_lower,
        "rna-seq"
      ) ~ "RNA-seq",
      
      TRUE ~ assay_expected
    )
  )

# ------------------------------------------------------------
# 3. Curate treatment using sample-specific information
# Priority is given to explicit no-treatment/control wording
# ------------------------------------------------------------

curated <- curated |>
  mutate(
    treatment_curated = case_when(
      
      # Explicit vancomycin treatment in sample title
      str_detect(
        title_lower,
        "vancomycin|vancomcyin"
      ) &
        !str_detect(
          title_lower,
          "no vancomycin"
        ) ~ "Vancomycin",
      
      # Explicit untreated groups
      str_detect(
        title_lower,
        "no treatment|no vancomycin|control_30mins"
      ) ~ "Untreated",
      
      # CLASH background controls
      str_detect(
        title_lower,
        "negative control|untagged_bkg"
      ) ~ "Background control",
      
      # GSE158830 CLASH and mapping assays without antibiotic
      dataset == "GSE158830" &
        assay_curated %in% c(
          "RNase III-CLASH",
          "dRNA-seq",
          "Term-seq"
        ) ~ "No antibiotic",
      
      TRUE ~ "Needs manual review"
    ),
    
    vancomycin_status = case_when(
      treatment_curated == "Vancomycin" ~ "Treated",
      treatment_curated %in% c(
        "Untreated",
        "No antibiotic"
      ) ~ "Untreated",
      treatment_curated == "Background control" ~ "Technical control",
      TRUE ~ "Unclear"
    )
  )

# ------------------------------------------------------------
# 4. Curate concentration and exposure time
# ------------------------------------------------------------

curated <- curated |>
  mutate(
    concentration_curated = case_when(
      
      dataset %in% c(
        "GSE254530",
        "GSE254532"
      ) &
        vancomycin_status == "Treated" ~ "8 ug/mL",
      
      dataset == "GSE158830" &
        vancomycin_status == "Treated" ~ "3 ug/mL",
      
      TRUE ~ NA_character_
    ),
    
    exposure_time_curated = case_when(
      
      dataset %in% c(
        "GSE254530",
        "GSE254532"
      ) &
        vancomycin_status == "Treated" ~ "30 min",
      
      dataset == "GSE158830" &
        vancomycin_status == "Treated" ~ "10 min",
      
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------
# 5. Extract biological replicate from title
# ------------------------------------------------------------

curated <- curated |>
  mutate(
    replicate_curated = case_when(
      
      str_detect(
        title_lower,
        "replicate\\s*[0-9]+"
      ) ~ str_extract(
        title_lower,
        "(?<=replicate\\s)[0-9]+"
      ),
      
      str_detect(
        title_lower,
        "_[1-9]$"
      ) ~ str_extract(
        title_lower,
        "[1-9]$"
      ),
      
      str_detect(
        title_lower,
        "30mins_[1-9]"
      ) ~ str_extract(
        title_lower,
        "(?<=30mins_)[1-9]"
      ),
      
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------
# 6. Curate genotype / experimental condition
# ------------------------------------------------------------

curated <- curated |>
  mutate(
    genotype_curated = case_when(
      
      str_detect(
        title_lower,
        "∆srna275::srna275|repaired"
      ) ~ "vigR_repaired",
      
      str_detect(
        title_lower,
        "∆srna275|∆vigr"
      ) ~ "vigR_deletion",
      
      str_detect(
        title_lower,
        "crispri knock-down|psd-1::srna275"
      ) ~ "vigR_CRISPRi",
      
      str_detect(
        title_lower,
        "empty vector"
      ) ~ "empty_vector",
      
      str_detect(
        title_lower,
        "\\bwt\\b"
      ) ~ "WT",
      
      str_detect(
        title_lower,
        "rnc-htf"
      ) ~ "RNaseIII_tagged",
      
      str_detect(
        title_lower,
        "untagged"
      ) ~ "untagged_control",
      
      TRUE ~ "Not specified"
    )
  )

# ------------------------------------------------------------
# 7. Define which samples serve each planned analysis
# ------------------------------------------------------------

curated <- curated |>
  mutate(
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
      
      dataset == "GSE158830" &
        assay_curated == "RNA-seq" &
        genotype_curated %in% c(
          "WT",
          "vigR_deletion",
          "vigR_CRISPRi",
          "empty_vector",
          "vigR_repaired"
        ) ~
        "Validation_vigR_RNAseq",
      
      dataset == "GSE158830" &
        assay_curated == "RNase III-CLASH" ~
        "Validation_CLASH_network",
      
      dataset == "GSE158830" &
        assay_curated %in% c(
          "dRNA-seq",
          "Term-seq"
        ) ~
        "Transcript_annotation_support",
      
      TRUE ~ "Exclude_or_review"
    )
  )

# ------------------------------------------------------------
# 8. Create summaries
# ------------------------------------------------------------

sample_design <- curated |>
  select(
    dataset,
    sample_id,
    gsm_accession,
    srx_accession,
    biosample_accession,
    sample_title,
    organism,
    strain_inferred,
    assay_curated,
    genotype_curated,
    treatment_curated,
    vancomycin_status,
    concentration_curated,
    exposure_time_curated,
    replicate_curated,
    analysis_module,
    supplementary_files
  ) |>
  arrange(
    dataset,
    analysis_module,
    genotype_curated,
    vancomycin_status,
    replicate_curated
  )

design_summary <- sample_design |>
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

manual_review <- sample_design |>
  filter(
    vancomycin_status == "Unclear" |
      analysis_module == "Exclude_or_review"
  )

# ------------------------------------------------------------
# 9. Save files
# ------------------------------------------------------------

write_csv(
  sample_design,
  file.path(
    output_folder,
    "BacRegRNA_curated_sample_design.csv"
  )
)

write_csv(
  design_summary,
  file.path(
    output_folder,
    "BacRegRNA_curated_design_summary.csv"
  )
)

write_csv(
  manual_review,
  file.path(
    output_folder,
    "BacRegRNA_samples_for_manual_review.csv"
  )
)

cat("\n============================================\n")
cat("STEP 06 COMPLETED\n")
cat("============================================\n")

cat(
  "Total samples:",
  nrow(sample_design),
  "\n"
)

cat(
  "Samples requiring manual review:",
  nrow(manual_review),
  "\n\n"
)

print(
  design_summary,
  n = Inf
)

cat(
  "\nFiles saved inside:\n",
  output_folder,
  "\n"
)