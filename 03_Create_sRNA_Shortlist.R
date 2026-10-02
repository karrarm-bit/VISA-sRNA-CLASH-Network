# ============================================================
# BacRegRNA Project
# Script 03: Create a shortlist of relevant sRNA studies
# ============================================================

library(readr)
library(dplyr)
library(stringr)

project_folder <- "D:/Bac-sRNA"

input_file <- file.path(
  project_folder,
  "02_screening",
  "Staphylococcus_aureus_metadata_repaired.csv"
)

output_folder <- file.path(
  project_folder,
  "03_selected_datasets"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

records <- read_csv(
  input_file,
  show_col_types = FALSE
)

shortlist <- records |>
  mutate(
    full_text = str_to_lower(
      paste(
        title_new,
        metadata_text,
        srna_terms,
        resistance_terms
      )
    ),
    
    study_accession = str_extract(
      paste(accession_new, metadata_text),
      "\\bGSE[0-9]+\\b"
    ),
    
    has_srna_evidence = str_detect(
      full_text,
      paste0(
        "small regulatory rna|",
        "small rna interactome|",
        "\\bsrna\\b|",
        "rnase iii.?clash|",
        "rna.?rna interaction|",
        "regulatory rna"
      )
    ),
    
    has_amr_evidence = str_detect(
      full_text,
      paste0(
        "vancomycin|",
        "methicillin|",
        "antibiotic resistance|",
        "antimicrobial resistance|",
        "drug resistance|",
        "tolerance|",
        "\\bmrsa\\b"
      )
    ),
    
    priority = case_when(
      has_srna_evidence & has_amr_evidence ~ "Priority 1",
      has_srna_evidence ~ "Priority 2",
      has_amr_evidence ~ "Priority 3",
      TRUE ~ "Exclude"
    )
  ) |>
  filter(priority != "Exclude") |>
  arrange(
    priority,
    desc(times_retrieved)
  )

# Unique GEO studies only
geo_shortlist <- shortlist |>
  filter(
    database == "gds",
    !is.na(study_accession)
  ) |>
  distinct(
    study_accession,
    .keep_all = TRUE
  ) |>
  select(
    priority,
    study_accession,
    title_new,
    date_new,
    srna_relevance,
    amr_relevance,
    resistance_category,
    times_retrieved
  )

write_csv(
  shortlist,
  file.path(
    output_folder,
    "Staphylococcus_sRNA_all_candidate_records.csv"
  )
)

write_csv(
  geo_shortlist,
  file.path(
    output_folder,
    "Staphylococcus_sRNA_unique_GEO_shortlist.csv"
  )
)

cat("\nStep 3 completed successfully.\n")
cat("Candidate records:", nrow(shortlist), "\n")
cat("Unique GEO studies:", nrow(geo_shortlist), "\n\n")

print(geo_shortlist)