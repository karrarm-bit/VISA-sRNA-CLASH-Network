# ============================================================
# BacRegRNA Project
# Script 04: Build GEO sample-level metadata
# Datasets:
# GSE254530 = RNA-seq
# GSE254531 = Ribo-seq
# GSE254532 = RNase III-CLASH
# GSE158830 = External validation
# ============================================================

# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

output_folder <- file.path(
  project_folder,
  "04_metadata"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

setwd(project_folder)

# ------------------------------------------------------------
# 2. Install required packages
# ------------------------------------------------------------

cran_packages <- c(
  "dplyr",
  "purrr",
  "stringr",
  "tidyr",
  "tibble",
  "readr"
)

missing_cran <- cran_packages[
  !cran_packages %in% rownames(installed.packages())
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install(
    "GEOquery",
    ask = FALSE,
    update = FALSE
  )
}

# ------------------------------------------------------------
# 3. Load packages
# ------------------------------------------------------------

library(GEOquery)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(tibble)
library(readr)

# ------------------------------------------------------------
# 4. Define datasets and assays
# ------------------------------------------------------------

dataset_table <- tribble(
  ~dataset,    ~assay,             ~role,
  "GSE254530", "RNA-seq",          "Discovery",
  "GSE254531", "Ribo-seq",         "Discovery",
  "GSE254532", "RNase III-CLASH",  "Discovery",
  "GSE158830", "RNase III-CLASH",  "Validation"
)

write_csv(
  dataset_table,
  file.path(
    output_folder,
    "BacRegRNA_dataset_list.csv"
  )
)

# ------------------------------------------------------------
# 5. Helper: safely collapse GEO metadata fields
# ------------------------------------------------------------

collapse_field <- function(x) {
  
  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }
  
  value <- paste(
    as.character(x),
    collapse = " | "
  )
  
  value <- str_squish(value)
  
  if (is.na(value) || value == "") {
    return(NA_character_)
  }
  
  value
}

# ------------------------------------------------------------
# 6. Helper: extract value from characteristic text
# Example:
# "treatment: vancomycin" -> "vancomycin"
# ------------------------------------------------------------

extract_characteristic <- function(text, patterns) {
  
  if (
    is.na(text) ||
    text == ""
  ) {
    return(NA_character_)
  }
  
  pattern <- paste(
    patterns,
    collapse = "|"
  )
  
  pieces <- unlist(
    str_split(
      text,
      "\\s*\\|\\s*|;|\\n"
    )
  )
  
  matched <- pieces[
    str_detect(
      str_to_lower(pieces),
      pattern
    )
  ]
  
  if (length(matched) == 0) {
    return(NA_character_)
  }
  
  value <- str_remove(
    matched[1],
    "^[^:=]+\\s*[:=]\\s*"
  )
  
  str_squish(value)
}

# ------------------------------------------------------------
# 7. Retrieve one GEO Series and its samples
# GSEMatrix = FALSE preserves GSM-level metadata
# ------------------------------------------------------------

fetch_geo_samples <- function(
    gse_accession,
    assay_name,
    study_role
) {
  
  cat(
    "\n============================================\n",
    "Downloading: ", gse_accession, "\n",
    "Assay: ", assay_name, "\n",
    "============================================\n",
    sep = ""
  )
  
  gse_object <- tryCatch(
    getGEO(
      GEO = gse_accession,
      GSEMatrix = FALSE,
      getGPL = FALSE,
      destdir = output_folder
    ),
    error = function(e) {
      
      message(
        "Failed to download ",
        gse_accession,
        ": ",
        conditionMessage(e)
      )
      
      return(NULL)
    }
  )
  
  if (is.null(gse_object)) {
    return(tibble())
  }
  
  gsm_list <- GSMList(gse_object)
  
  if (length(gsm_list) == 0) {
    
    message(
      "No GSM samples were found for ",
      gse_accession
    )
    
    return(tibble())
  }
  
  cat(
    "Number of GSM samples found:",
    length(gsm_list),
    "\n"
  )
  
  sample_table <- imap_dfr(
    gsm_list,
    function(gsm_object, gsm_name) {
      
      metadata <- Meta(gsm_object)
      
      characteristics <- collapse_field(
        metadata$characteristics_ch1
      )
      
      source_name <- collapse_field(
        metadata$source_name_ch1
      )
      
      title <- collapse_field(
        metadata$title
      )
      
      description <- collapse_field(
        metadata$description
      )
      
      molecule <- collapse_field(
        metadata$molecule_ch1
      )
      
      extract_protocol <- collapse_field(
        metadata$extract_protocol_ch1
      )
      
      library_strategy <- collapse_field(
        metadata$library_strategy
      )
      
      library_source <- collapse_field(
        metadata$library_source
      )
      
      library_selection <- collapse_field(
        metadata$library_selection
      )
      
      organism <- collapse_field(
        metadata$organism_ch1
      )
      
      platform <- collapse_field(
        metadata$platform_id
      )
      
      relation <- collapse_field(
        metadata$relation
      )
      
      supplementary_files <- collapse_field(
        metadata$supplementary_file
      )
      
      all_text <- str_to_lower(
        paste(
          title,
          source_name,
          characteristics,
          description,
          extract_protocol,
          library_strategy,
          library_source,
          library_selection
        )
      )
      
      treatment <- extract_characteristic(
        characteristics,
        c(
          "\\btreatment\\b",
          "\\bcondition\\b",
          "\\bexposure\\b",
          "\\bdrug\\b",
          "\\bvancomycin\\b"
        )
      )
      
      strain <- extract_characteristic(
        characteristics,
        c(
          "\\bstrain\\b",
          "\\bisolate\\b",
          "\\bgenotype\\b"
        )
      )
      
      replicate <- extract_characteristic(
        characteristics,
        c(
          "\\breplicate\\b",
          "\\bbiological replicate\\b",
          "\\brep\\b"
        )
      )
      
      concentration <- extract_characteristic(
        characteristics,
        c(
          "\\bconcentration\\b",
          "\\bdose\\b",
          "µg/ml",
          "ug/ml",
          "mg/l"
        )
      )
      
      exposure_time <- extract_characteristic(
        characteristics,
        c(
          "\\btime\\b",
          "\\bduration\\b",
          "\\bminutes?\\b",
          "\\bhours?\\b"
        )
      )
      
      inferred_group <- case_when(
        
        str_detect(
          all_text,
          "vancomycin|vanc-treated|van-treated"
        ) ~ "Vancomycin",
        
        str_detect(
          all_text,
          "untreated|no treatment|control|vehicle"
        ) ~ "Control",
        
        TRUE ~ "Unclear"
      )
      
      inferred_assay <- case_when(
        
        str_detect(
          all_text,
          "clash|rnase iii"
        ) ~ "RNase III-CLASH",
        
        str_detect(
          all_text,
          "ribo.?seq|ribosome profiling"
        ) ~ "Ribo-seq",
        
        str_detect(
          all_text,
          "rna.?seq|transcriptom"
        ) ~ "RNA-seq",
        
        TRUE ~ assay_name
      )
      
      raw_data_available <- case_when(
        
        !is.na(relation) &&
          str_detect(
            str_to_lower(relation),
            "sra|srx|srr"
          ) ~ "Yes",
        
        TRUE ~ "Unclear"
      )
      
      processed_data_available <- case_when(
        
        !is.na(supplementary_files) &&
          supplementary_files != "" ~ "Yes",
        
        TRUE ~ "No or unclear"
      )
      
      tibble(
        dataset = gse_accession,
        study_role = study_role,
        sample_id = gsm_name,
        sample_title = title,
        assay_expected = assay_name,
        assay_inferred = inferred_assay,
        organism = organism,
        strain = strain,
        treatment_reported = treatment,
        treatment_group_inferred = inferred_group,
        concentration = concentration,
        exposure_time = exposure_time,
        biological_replicate = replicate,
        source_name = source_name,
        characteristics = characteristics,
        molecule = molecule,
        library_strategy = library_strategy,
        library_source = library_source,
        library_selection = library_selection,
        platform = platform,
        raw_data_available = raw_data_available,
        processed_data_available = processed_data_available,
        sra_relation = relation,
        supplementary_files = supplementary_files,
        description = description
      )
    }
  )
  
  sample_table
}

# ------------------------------------------------------------
# 8. Download all four datasets
# ------------------------------------------------------------

all_metadata <- pmap_dfr(
  dataset_table,
  function(dataset, assay, role) {
    
    Sys.sleep(1)
    
    fetch_geo_samples(
      gse_accession = dataset,
      assay_name = assay,
      study_role = role
    )
  }
)

# ------------------------------------------------------------
# 9. Stop if nothing was retrieved
# ------------------------------------------------------------

if (nrow(all_metadata) == 0) {
  
  stop(
    "No GEO sample metadata were retrieved. ",
    "Check internet connection and GEO accessions."
  )
}

# ------------------------------------------------------------
# 10. Add manual-review flags
# ------------------------------------------------------------

all_metadata <- all_metadata |>
  mutate(
    
    needs_group_review = case_when(
      treatment_group_inferred == "Unclear" ~ "Yes",
      TRUE ~ "No"
    ),
    
    needs_strain_review = case_when(
      is.na(strain) | strain == "" ~ "Yes",
      TRUE ~ "No"
    ),
    
    needs_replicate_review = case_when(
      is.na(biological_replicate) |
        biological_replicate == "" ~ "Yes",
      TRUE ~ "No"
    ),
    
    metadata_complete = case_when(
      needs_group_review == "No" &
        needs_strain_review == "No" &
        needs_replicate_review == "No" ~ "Yes",
      
      TRUE ~ "Incomplete"
    )
  ) |>
  arrange(
    dataset,
    treatment_group_inferred,
    sample_id
  )

# ------------------------------------------------------------
# 11. Build dataset-level summary
# ------------------------------------------------------------

dataset_summary <- all_metadata |>
  group_by(
    dataset,
    study_role,
    assay_expected
  ) |>
  summarise(
    total_samples = n(),
    control_samples = sum(
      treatment_group_inferred == "Control",
      na.rm = TRUE
    ),
    vancomycin_samples = sum(
      treatment_group_inferred == "Vancomycin",
      na.rm = TRUE
    ),
    unclear_samples = sum(
      treatment_group_inferred == "Unclear",
      na.rm = TRUE
    ),
    unique_strains = paste(
      sort(
        unique(
          na.omit(strain)
        )
      ),
      collapse = "; "
    ),
    raw_data_yes = sum(
      raw_data_available == "Yes",
      na.rm = TRUE
    ),
    processed_data_yes = sum(
      processed_data_available == "Yes",
      na.rm = TRUE
    ),
    incomplete_metadata = sum(
      metadata_complete == "Incomplete",
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 12. Build group-level design summary
# ------------------------------------------------------------

design_summary <- all_metadata |>
  count(
    dataset,
    assay_expected,
    treatment_group_inferred,
    strain,
    concentration,
    exposure_time,
    name = "number_of_samples"
  ) |>
  arrange(
    dataset,
    treatment_group_inferred
  )

# ------------------------------------------------------------
# 13. Extract potential SRA accessions
# ------------------------------------------------------------

sra_accessions <- all_metadata |>
  mutate(
    sra_accession = str_extract(
      sra_relation,
      "\\b(SRP|SRR|SRX|SRS)[0-9]+\\b"
    )
  ) |>
  filter(!is.na(sra_accession)) |>
  distinct(
    dataset,
    sample_id,
    sra_accession,
    .keep_all = TRUE
  ) |>
  select(
    dataset,
    sample_id,
    sample_title,
    assay_expected,
    treatment_group_inferred,
    sra_accession,
    sra_relation
  )

# ------------------------------------------------------------
# 14. Save output files
# ------------------------------------------------------------

write_csv(
  all_metadata,
  file.path(
    output_folder,
    "BacRegRNA_all_GEO_sample_metadata.csv"
  )
)

write_csv(
  dataset_summary,
  file.path(
    output_folder,
    "BacRegRNA_dataset_summary.csv"
  )
)

write_csv(
  design_summary,
  file.path(
    output_folder,
    "BacRegRNA_experimental_design_summary.csv"
  )
)

write_csv(
  sra_accessions,
  file.path(
    output_folder,
    "BacRegRNA_SRA_accessions.csv"
  )
)

# Save records requiring manual review
manual_review <- all_metadata |>
  filter(
    metadata_complete == "Incomplete"
  )

write_csv(
  manual_review,
  file.path(
    output_folder,
    "BacRegRNA_metadata_manual_review.csv"
  )
)

# ------------------------------------------------------------
# 15. Print results
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 04 COMPLETED\n")
cat("============================================\n")

cat(
  "Total samples retrieved:",
  nrow(all_metadata),
  "\n"
)

cat(
  "Datasets retrieved:",
  paste(
    unique(all_metadata$dataset),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Samples needing manual review:",
  nrow(manual_review),
  "\n"
)

cat(
  "\nFiles saved inside:\n",
  output_folder,
  "\n\n"
)

print(dataset_summary)

cat("\nExperimental design summary:\n")
print(design_summary)

cat("\nSRA accessions detected:", nrow(sra_accessions), "\n")