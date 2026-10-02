# ============================================================
# BacRegRNA Project
# Script 04D: Download GEO metadata without FTP
# ============================================================

project_folder <- "D:/Bac-sRNA"
output_folder <- file.path(project_folder, "04_metadata")
soft_folder <- file.path(output_folder, "SOFT_files")

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(soft_folder, recursive = TRUE, showWarnings = FALSE)

setwd(project_folder)
options(timeout = 1200)

# ------------------------------------------------------------
# Install packages
# ------------------------------------------------------------

cran_packages <- c(
  "httr2",
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

library(httr2)
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
  ~dataset,    ~assay_expected,   ~study_role,
  "GSE254530", "RNA-seq",         "Discovery",
  "GSE254531", "Ribo-seq",        "Discovery",
  "GSE254532", "RNase III-CLASH", "Discovery",
  "GSE158830", "RNase III-CLASH", "Validation"
)

# ------------------------------------------------------------
# Download SOFT file directly from GEO query page
# ------------------------------------------------------------

download_geo_soft <- function(gse) {
  
  soft_file <- file.path(
    soft_folder,
    paste0(gse, "_family.soft")
  )
  
  if (file.exists(soft_file) && file.info(soft_file)$size > 1000) {
    cat("Existing SOFT file found:", gse, "\n")
    return(soft_file)
  }
  
  url <- paste0(
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?",
    "acc=", gse,
    "&targ=self",
    "&form=text",
    "&view=full"
  )
  
  cat("\nDownloading:", gse, "\n")
  cat("From GEO webpage instead of FTP.\n")
  
  success <- tryCatch({
    
    response <- request(url) |>
      req_user_agent(
        "Mozilla/5.0 BacRegRNA academic research"
      ) |>
      req_headers(
        Accept = "text/plain"
      ) |>
      req_timeout(300) |>
      req_retry(max_tries = 3) |>
      req_perform()
    
    if (resp_status(response) != 200) {
      stop(
        "HTTP status: ",
        resp_status(response)
      )
    }
    
    content <- resp_body_raw(response)
    
    writeBin(
      content,
      soft_file
    )
    
    TRUE
    
  }, error = function(e) {
    
    message(
      "Download failed for ",
      gse,
      ": ",
      conditionMessage(e)
    )
    
    FALSE
  })
  
  if (
    !success ||
    !file.exists(soft_file) ||
    file.info(soft_file)$size < 1000
  ) {
    return(NA_character_)
  }
  
  cat(
    "Downloaded successfully:",
    basename(soft_file),
    "\n"
  )
  
  soft_file
}

# ------------------------------------------------------------
# Parse SOFT file using GEOquery
# ------------------------------------------------------------

parse_geo_soft <- function(
    soft_file,
    dataset,
    assay_expected,
    study_role
) {
  
  if (is.na(soft_file) || !file.exists(soft_file)) {
    return(tibble())
  }
  
  cat("Parsing:", dataset, "\n")
  
  gse_object <- tryCatch(
    
    getGEO(
      filename = soft_file,
      GSEMatrix = FALSE,
      getGPL = FALSE
    ),
    
    error = function(e) {
      
      message(
        "Parsing failed for ",
        dataset,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(gse_object)) {
    return(tibble())
  }
  
  gsm_list <- GSMList(gse_object)
  
  if (length(gsm_list) == 0) {
    message("No GSM samples found for ", dataset)
    return(tibble())
  }
  
  cat(
    "Samples detected:",
    length(gsm_list),
    "\n"
  )
  
  imap_dfr(
    gsm_list,
    function(gsm_object, gsm_id) {
      
      metadata <- Meta(gsm_object)
      
      collapse_value <- function(x) {
        
        if (is.null(x) || length(x) == 0) {
          return(NA_character_)
        }
        
        value <- paste(
          as.character(x),
          collapse = " | "
        )
        
        str_squish(value)
      }
      
      title <- collapse_value(metadata$title)
      source <- collapse_value(metadata$source_name_ch1)
      characteristics <- collapse_value(
        metadata$characteristics_ch1
      )
      treatment_protocol <- collapse_value(
        metadata$treatment_protocol_ch1
      )
      growth_protocol <- collapse_value(
        metadata$growth_protocol_ch1
      )
      extract_protocol <- collapse_value(
        metadata$extract_protocol_ch1
      )
      organism <- collapse_value(
        metadata$organism_ch1
      )
      platform <- collapse_value(
        metadata$platform_id
      )
      molecule <- collapse_value(
        metadata$molecule_ch1
      )
      relation <- collapse_value(
        metadata$relation
      )
      supplementary <- collapse_value(
        metadata$supplementary_file
      )
      
      combined_text <- str_to_lower(
        str_squish(
          paste(
            title,
            source,
            characteristics,
            treatment_protocol,
            growth_protocol,
            extract_protocol,
            relation
          )
        )
      )
      
      treatment_group <- case_when(
        
        str_detect(
          combined_text,
          "untreated|un-treated|without vancomycin|no vancomycin|control"
        ) ~ "Control",
        
        str_detect(
          combined_text,
          "vancomycin|vanco"
        ) ~ "Vancomycin",
        
        TRUE ~ "Unclear"
      )
      
      strain_inferred <- case_when(
        
        str_detect(combined_text, "jkd6008") ~ "JKD6008",
        str_detect(combined_text, "jkd6009") ~ "JKD6009",
        str_detect(combined_text, "usa300") ~ "USA300",
        str_detect(combined_text, "newman") ~ "Newman",
        str_detect(combined_text, "hg003") ~ "HG003",
        
        TRUE ~ str_match(
          combined_text,
          "strain\\s*[:=]\\s*([^|;,]+)"
        )[, 2]
      )
      
      concentration_inferred <- str_extract(
        combined_text,
        "[0-9.]+\\s*(µg/ml|μg/ml|ug/ml|mg/l)"
      )
      
      exposure_time_inferred <- str_extract(
        combined_text,
        "[0-9.]+\\s*(min|mins|minute|minutes|hour|hours|hr|hrs)"
      )
      
      replicate_inferred <- str_extract(
        combined_text,
        "(biological\\s*)?(replicate|rep)\\s*[:=_-]?\\s*[0-9]+"
      )
      
      sra_accession <- str_extract(
        combined_text,
        "\\b(SRR|SRX|SRS|SRP)[0-9]+\\b"
      )
      
      tibble(
        dataset = dataset,
        study_role = study_role,
        assay_expected = assay_expected,
        sample_id = gsm_id,
        sample_title = title,
        organism = organism,
        source = source,
        strain_inferred = strain_inferred,
        treatment_group = treatment_group,
        concentration_inferred = concentration_inferred,
        exposure_time_inferred = exposure_time_inferred,
        replicate_inferred = replicate_inferred,
        molecule = molecule,
        platform = platform,
        characteristics = characteristics,
        treatment_protocol = treatment_protocol,
        growth_protocol = growth_protocol,
        extract_protocol = extract_protocol,
        sra_accession = sra_accession,
        relations = relation,
        supplementary_files = supplementary,
        combined_metadata = combined_text
      )
    }
  )
}

# ------------------------------------------------------------
# Download and parse all datasets
# ------------------------------------------------------------

all_metadata <- pmap_dfr(
  dataset_table,
  function(dataset, assay_expected, study_role) {
    
    soft_file <- download_geo_soft(dataset)
    
    Sys.sleep(1)
    
    parse_geo_soft(
      soft_file = soft_file,
      dataset = dataset,
      assay_expected = assay_expected,
      study_role = study_role
    )
  }
)

# ------------------------------------------------------------
# Stop if retrieval failed
# ------------------------------------------------------------

if (nrow(all_metadata) == 0) {
  stop(
    "No sample metadata were retrieved. ",
    "Check whether the downloaded SOFT files contain GEO text."
  )
}

# ------------------------------------------------------------
# Add review flags
# ------------------------------------------------------------

all_metadata <- all_metadata |>
  mutate(
    metadata_status = case_when(
      
      treatment_group != "Unclear" &
        !is.na(strain_inferred) ~ "Partly complete",
      
      TRUE ~ "Needs review"
    )
  ) |>
  arrange(
    dataset,
    treatment_group,
    sample_id
  )

# ------------------------------------------------------------
# Create summaries
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
      treatment_group == "Control",
      na.rm = TRUE
    ),
    
    vancomycin_samples = sum(
      treatment_group == "Vancomycin",
      na.rm = TRUE
    ),
    
    unclear_samples = sum(
      treatment_group == "Unclear",
      na.rm = TRUE
    ),
    
    detected_strains = paste(
      sort(unique(na.omit(strain_inferred))),
      collapse = "; "
    ),
    
    .groups = "drop"
  )

design_summary <- all_metadata |>
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
  filter(metadata_status == "Needs review")

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

write_csv(
  all_metadata,
  file.path(
    output_folder,
    "BacRegRNA_SOFT_all_sample_metadata.csv"
  )
)

write_csv(
  dataset_summary,
  file.path(
    output_folder,
    "BacRegRNA_SOFT_dataset_summary.csv"
  )
)

write_csv(
  design_summary,
  file.path(
    output_folder,
    "BacRegRNA_SOFT_design_summary.csv"
  )
)

write_csv(
  manual_review,
  file.path(
    output_folder,
    "BacRegRNA_SOFT_manual_review.csv"
  )
)

# ------------------------------------------------------------
# Final console output
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 04D COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat(
  "Total samples:",
  nrow(all_metadata),
  "\n"
)

print(dataset_summary)

cat("\nExperimental design:\n")
print(design_summary)

cat(
  "\nFiles saved in:\n",
  output_folder,
  "\n"
)
install.packages("curl")
library(curl)

dir.create(
  "D:/Bac-sRNA/04_metadata/manual_download",
  recursive = TRUE,
  showWarnings = FALSE
)

url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE254nnn/GSE254530/miniml/",
  "GSE254530_family.xml.tgz"
)

destination <- paste0(
  "D:/Bac-sRNA/04_metadata/manual_download/",
  "GSE254530_family.xml.tgz"
)

handle <- new_handle(
  useragent = paste0(
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) ",
    "AppleWebKit/537.36 Chrome/150 Safari/537.36"
  )
)

curl_download(
  url = url,
  destfile = destination,
  handle = handle,
  quiet = FALSE
)

file.info(destination)[, c("size", "mtime")]