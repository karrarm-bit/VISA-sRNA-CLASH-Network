# ============================================================
# BacRegRNA Project
# Script 04B: Build GEO metadata using Series Matrix
# ============================================================

project_folder <- "D:/Bac-sRNA"
output_folder <- file.path(project_folder, "04_metadata")

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

setwd(project_folder)

# ------------------------------------------------------------
# Install packages
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
# Load packages
# ------------------------------------------------------------

library(GEOquery)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(tibble)
library(readr)

# ------------------------------------------------------------
# Dataset definitions
# ------------------------------------------------------------

dataset_table <- tribble(
  ~dataset,    ~assay,            ~study_role,
  "GSE254530", "RNA-seq",         "Discovery",
  "GSE254531", "Ribo-seq",        "Discovery",
  "GSE254532", "RNase III-CLASH", "Discovery",
  "GSE158830", "RNase III-CLASH", "Validation"
)

# ------------------------------------------------------------
# Function to retrieve one GEO dataset
# ------------------------------------------------------------

fetch_series_matrix <- function(dataset, assay, study_role) {
  
  cat("\n--------------------------------------------\n")
  cat("Downloading:", dataset, "\n")
  cat("Assay:", assay, "\n")
  cat("--------------------------------------------\n")
  
  gse_list <- tryCatch(
    getGEO(
      dataset,
      GSEMatrix = TRUE,
      getGPL = FALSE,
      AnnotGPL = FALSE,
      destdir = output_folder
    ),
    error = function(e) {
      message(
        "Download failed for ",
        dataset,
        ": ",
        conditionMessage(e)
      )
      return(NULL)
    }
  )
  
  if (is.null(gse_list)) {
    return(tibble())
  }
  
  # GEOquery usually returns a list of ExpressionSet objects
  if (!is.list(gse_list)) {
    gse_list <- list(gse_list)
  }
  
  map_dfr(
    seq_along(gse_list),
    function(i) {
      
      eset <- gse_list[[i]]
      
      metadata <- Biobase::pData(eset) |>
        as.data.frame() |>
        rownames_to_column("sample_id") |>
        as_tibble()
      
      if (nrow(metadata) == 0) {
        return(tibble())
      }
      
      metadata |>
        mutate(
          dataset = dataset,
          assay_expected = assay,
          study_role = study_role,
          platform_index = i,
          .before = 1
        )
    }
  )
}

# ------------------------------------------------------------
# Download all datasets
# ------------------------------------------------------------

all_raw_metadata <- pmap_dfr(
  dataset_table,
  function(dataset, assay, study_role) {
    
    Sys.sleep(1)
    
    fetch_series_matrix(
      dataset = dataset,
      assay = assay,
      study_role = study_role
    )
  }
)

# ------------------------------------------------------------
# Verify retrieval
# ------------------------------------------------------------

if (nrow(all_raw_metadata) == 0) {
  stop(
    "No GEO metadata were retrieved. ",
    "Run the single-dataset test shown below."
  )
}

cat(
  "\nMetadata successfully retrieved for",
  nrow(all_raw_metadata),
  "samples.\n"
)

# ------------------------------------------------------------
# Combine all textual metadata columns
# ------------------------------------------------------------

text_columns <- names(all_raw_metadata)[
  vapply(all_raw_metadata, is.character, logical(1))
]

all_metadata <- all_raw_metadata |>
  mutate(
    combined_metadata = apply(
      select(., all_of(text_columns)),
      1,
      function(x) {
        str_squish(
          paste(
            x[!is.na(x)],
            collapse = " | "
          )
        )
      }
    ),
    
    metadata_lower = str_to_lower(combined_metadata)
  )

# ------------------------------------------------------------
# Infer treatment group
# ------------------------------------------------------------

all_metadata <- all_metadata |>
  mutate(
    treatment_group = case_when(
      
      str_detect(
        metadata_lower,
        "untreated|un-treated|without vancomycin|no vancomycin|control"
      ) ~ "Control",
      
      str_detect(
        metadata_lower,
        "vancomycin|vanco"
      ) ~ "Vancomycin",
      
      TRUE ~ "Unclear"
    ),
    
    strain_inferred = case_when(
      
      str_detect(metadata_lower, "jkd6008") ~ "JKD6008",
      
      str_detect(metadata_lower, "jkd6009") ~ "JKD6009",
      
      str_detect(metadata_lower, "usa300") ~ "USA300",
      
      str_detect(metadata_lower, "newman") ~ "Newman",
      
      str_detect(metadata_lower, "hg003") ~ "HG003",
      
      TRUE ~ NA_character_
    ),
    
    assay_inferred = case_when(
      
      str_detect(
        metadata_lower,
        "rnase iii.?clash|rnaseiii.?clash|clash"
      ) ~ "RNase III-CLASH",
      
      str_detect(
        metadata_lower,
        "ribo.?seq|ribosome profiling|ribosome footprint"
      ) ~ "Ribo-seq",
      
      str_detect(
        metadata_lower,
        "rna.?seq|transcriptom"
      ) ~ "RNA-seq",
      
      TRUE ~ assay_expected
    ),
    
    concentration_inferred = str_extract(
      metadata_lower,
      "[0-9.]+\\s*(µg/ml|ug/ml|μg/ml|mg/l)"
    ),
    
    exposure_time_inferred = str_extract(
      metadata_lower,
      "[0-9.]+\\s*(min|mins|minute|minutes|hour|hours|hr|hrs)"
    ),
    
    replicate_inferred = str_extract(
      metadata_lower,
      "(biological\\s*)?(replicate|rep|sample)\\s*[:=_-]?\\s*[0-9]+"
    )
  )

# ------------------------------------------------------------
# Summary tables
# ------------------------------------------------------------

dataset_summary <- all_metadata |>
  count(
    dataset,
    study_role,
    assay_expected,
    name = "total_samples"
  ) |>
  arrange(dataset)

group_summary <- all_metadata |>
  count(
    dataset,
    assay_expected,
    treatment_group,
    strain_inferred,
    concentration_inferred,
    exposure_time_inferred,
    name = "number_of_samples"
  ) |>
  arrange(
    dataset,
    treatment_group
  )

manual_review <- all_metadata |>
  filter(
    treatment_group == "Unclear" |
      is.na(strain_inferred) |
      is.na(replicate_inferred)
  )

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

write_csv(
  all_raw_metadata,
  file.path(
    output_folder,
    "BacRegRNA_raw_GEO_metadata.csv"
  )
)

write_csv(
  all_metadata,
  file.path(
    output_folder,
    "BacRegRNA_cleaned_GEO_metadata.csv"
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
  group_summary,
  file.path(
    output_folder,
    "BacRegRNA_group_summary.csv"
  )
)

write_csv(
  manual_review,
  file.path(
    output_folder,
    "BacRegRNA_metadata_manual_review.csv"
  )
)

# ------------------------------------------------------------
# Final console report
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 04B COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat("Total retrieved samples:", nrow(all_metadata), "\n")
cat(
  "Datasets:",
  paste(unique(all_metadata$dataset), collapse = ", "),
  "\n"
)
cat(
  "Samples requiring manual review:",
  nrow(manual_review),
  "\n\n"
)

print(dataset_summary)

cat("\nTreatment-group summary:\n")
print(group_summary)

cat(
  "\nFiles saved inside:\n",
  output_folder,
  "\n"
)