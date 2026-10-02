# ============================================================
# BacRegRNA Project
# Script 04C: Download and parse GEO MINiML sample metadata
# ============================================================

# ------------------------------------------------------------
# 1. Project folders
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"
output_folder <- file.path(project_folder, "04_metadata")
miniml_folder <- file.path(output_folder, "MINiML_files")

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(miniml_folder, recursive = TRUE, showWarnings = FALSE)

setwd(project_folder)

# Increase timeout for GEO downloads
options(timeout = 1200)

# Prefer a stable Windows download method
options(download.file.method = "libcurl")

# ------------------------------------------------------------
# 2. Install packages
# ------------------------------------------------------------

required_packages <- c(
  "xml2",
  "dplyr",
  "purrr",
  "stringr",
  "tidyr",
  "tibble",
  "readr"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# ------------------------------------------------------------
# 3. Load packages
# ------------------------------------------------------------

library(xml2)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(tibble)
library(readr)

# ------------------------------------------------------------
# 4. Datasets
# ------------------------------------------------------------

dataset_table <- tribble(
  ~dataset,    ~assay_expected,    ~study_role,
  "GSE254530", "RNA-seq",          "Discovery",
  "GSE254531", "Ribo-seq",         "Discovery",
  "GSE254532", "RNase III-CLASH",  "Discovery",
  "GSE158830", "RNase III-CLASH",  "Validation"
)

write_csv(
  dataset_table,
  file.path(output_folder, "BacRegRNA_dataset_list.csv")
)

# ------------------------------------------------------------
# 5. Build GEO FTP subfolder
# Example:
# GSE254530 -> GSE254nnn
# GSE158830 -> GSE158nnn
# ------------------------------------------------------------

geo_series_folder <- function(gse) {
  
  number <- as.integer(str_remove(gse, "^GSE"))
  prefix <- floor(number / 1000)
  
  paste0("GSE", prefix, "nnn")
}

# ------------------------------------------------------------
# 6. Safe XML helper functions
# ------------------------------------------------------------

safe_xml_text <- function(node, xpath) {
  
  result <- xml_find_first(node, xpath)
  
  if (inherits(result, "xml_missing")) {
    return(NA_character_)
  }
  
  value <- xml_text(result, trim = TRUE)
  
  if (length(value) == 0 || is.na(value) || value == "") {
    return(NA_character_)
  }
  
  value
}

safe_xml_multiple <- function(node, xpath) {
  
  results <- xml_find_all(node, xpath)
  
  if (length(results) == 0) {
    return(NA_character_)
  }
  
  values <- xml_text(results, trim = TRUE)
  values <- values[!is.na(values) & values != ""]
  
  if (length(values) == 0) {
    return(NA_character_)
  }
  
  paste(unique(values), collapse = " | ")
}

# ------------------------------------------------------------
# 7. Download one MINiML family file
# ------------------------------------------------------------

download_miniml <- function(gse) {
  
  series_folder <- geo_series_folder(gse)
  
  url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    series_folder,
    "/",
    gse,
    "/miniml/",
    gse,
    "_family.xml.tgz"
  )
  
  tgz_file <- file.path(
    miniml_folder,
    paste0(gse, "_family.xml.tgz")
  )
  
  extraction_folder <- file.path(
    miniml_folder,
    gse
  )
  
  dir.create(
    extraction_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  cat("\n============================================\n")
  cat("Dataset:", gse, "\n")
  cat("URL:", url, "\n")
  cat("============================================\n")
  
  if (!file.exists(tgz_file)) {
    
    status <- tryCatch(
      {
        download.file(
          url = url,
          destfile = tgz_file,
          mode = "wb",
          quiet = FALSE
        )
        
        TRUE
      },
      error = function(e) {
        
        message(
          "Download failed for ",
          gse,
          ": ",
          conditionMessage(e)
        )
        
        FALSE
      }
    )
    
    if (!status) {
      return(NA_character_)
    }
    
  } else {
    
    cat("Existing downloaded file detected; skipping download.\n")
  }
  
  if (
    !file.exists(tgz_file) ||
    file.info(tgz_file)$size == 0
  ) {
    message("Downloaded file is missing or empty: ", tgz_file)
    return(NA_character_)
  }
  
  extraction_status <- tryCatch(
    {
      untar(
        tgz_file,
        exdir = extraction_folder
      )
      
      TRUE
    },
    error = function(e) {
      
      message(
        "Extraction failed for ",
        gse,
        ": ",
        conditionMessage(e)
      )
      
      FALSE
    }
  )
  
  if (!extraction_status) {
    return(NA_character_)
  }
  
  xml_files <- list.files(
    extraction_folder,
    pattern = "\\.xml$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(xml_files) == 0) {
    message("No XML file found for ", gse)
    return(NA_character_)
  }
  
  xml_files[1]
}

# ------------------------------------------------------------
# 8. Parse one MINiML XML file
# ------------------------------------------------------------

parse_miniml_samples <- function(
    xml_file,
    dataset,
    assay_expected,
    study_role
) {
  
  if (
    is.na(xml_file) ||
    !file.exists(xml_file)
  ) {
    return(tibble())
  }
  
  cat("Parsing XML:", basename(xml_file), "\n")
  
  doc <- tryCatch(
    read_xml(xml_file),
    error = function(e) {
      
      message(
        "XML parsing failed for ",
        dataset,
        ": ",
        conditionMessage(e)
      )
      
      return(NULL)
    }
  )
  
  if (is.null(doc)) {
    return(tibble())
  }
  
  # local-name() avoids XML namespace problems
  sample_nodes <- xml_find_all(
    doc,
    ".//*[local-name()='Sample']"
  )
  
  if (length(sample_nodes) == 0) {
    message("No sample nodes found for ", dataset)
    return(tibble())
  }
  
  cat("Samples found:", length(sample_nodes), "\n")
  
  map_dfr(
    sample_nodes,
    function(sample_node) {
      
      sample_id <- xml_attr(sample_node, "iid")
      
      if (
        is.na(sample_id) ||
        sample_id == ""
      ) {
        
        sample_id <- safe_xml_text(
          sample_node,
          ".//*[local-name()='Accession']"
        )
      }
      
      title <- safe_xml_text(
        sample_node,
        ".//*[local-name()='Title']"
      )
      
      source <- safe_xml_text(
        sample_node,
        ".//*[local-name()='Source']"
      )
      
      organism <- safe_xml_text(
        sample_node,
        ".//*[local-name()='Organism']"
      )
      
      molecule <- safe_xml_text(
        sample_node,
        ".//*[local-name()='Molecule']"
      )
      
      description <- safe_xml_multiple(
        sample_node,
        ".//*[local-name()='Description']"
      )
      
      characteristics_nodes <- xml_find_all(
        sample_node,
        ".//*[local-name()='Characteristics']"
      )
      
      characteristics <- if (
        length(characteristics_nodes) == 0
      ) {
        
        NA_character_
        
      } else {
        
        map_chr(
          characteristics_nodes,
          function(x) {
            
            tag <- xml_attr(x, "tag")
            value <- xml_text(x, trim = TRUE)
            
            if (
              !is.na(tag) &&
              tag != ""
            ) {
              paste0(tag, ": ", value)
            } else {
              value
            }
          }
        ) |>
          paste(collapse = " | ")
      }
      
      extract_protocol <- safe_xml_multiple(
        sample_node,
        ".//*[local-name()='Extract-Protocol']"
      )
      
      treatment_protocol <- safe_xml_multiple(
        sample_node,
        ".//*[local-name()='Treatment-Protocol']"
      )
      
      growth_protocol <- safe_xml_multiple(
        sample_node,
        ".//*[local-name()='Growth-Protocol']"
      )
      
      library_strategy <- safe_xml_text(
        sample_node,
        ".//*[local-name()='Library-Strategy']"
      )
      
      library_source <- safe_xml_text(
        sample_node,
        ".//*[local-name()='Library-Source']"
      )
      
      library_selection <- safe_xml_text(
        sample_node,
        ".//*[local-name()='Library-Selection']"
      )
      
      platform_ref_node <- xml_find_first(
        sample_node,
        ".//*[local-name()='Platform-Ref']"
      )
      
      platform <- if (
        inherits(platform_ref_node, "xml_missing")
      ) {
        NA_character_
      } else {
        xml_attr(platform_ref_node, "ref")
      }
      
      supplementary_files <- safe_xml_multiple(
        sample_node,
        ".//*[local-name()='Supplementary-Data']"
      )
      
      relation_values <- safe_xml_multiple(
        sample_node,
        ".//*[local-name()='Relation']"
      )
      
      all_text <- str_to_lower(
        str_squish(
          paste(
            title,
            source,
            characteristics,
            description,
            treatment_protocol,
            growth_protocol,
            extract_protocol,
            library_strategy,
            library_source,
            library_selection,
            relation_values
          )
        )
      )
      
      treatment_group <- case_when(
        
        str_detect(
          all_text,
          "untreated|un-treated|without vancomycin|no vancomycin|control"
        ) ~ "Control",
        
        str_detect(
          all_text,
          "vancomycin|vanco"
        ) ~ "Vancomycin",
        
        TRUE ~ "Unclear"
      )
      
      strain_inferred <- case_when(
        
        str_detect(all_text, "jkd6008") ~ "JKD6008",
        str_detect(all_text, "jkd6009") ~ "JKD6009",
        str_detect(all_text, "usa300") ~ "USA300",
        str_detect(all_text, "newman") ~ "Newman",
        str_detect(all_text, "hg003") ~ "HG003",
        
        TRUE ~ str_match(
          all_text,
          "strain\\s*[:=]\\s*([^|;,]+)"
        )[, 2]
      )
      
      concentration_inferred <- str_extract(
        all_text,
        "[0-9.]+\\s*(µg/ml|μg/ml|ug/ml|mg/l)"
      )
      
      exposure_time_inferred <- str_extract(
        all_text,
        "[0-9.]+\\s*(min|mins|minute|minutes|hour|hours|hr|hrs)"
      )
      
      replicate_inferred <- str_extract(
        all_text,
        "(biological\\s*)?(replicate|rep)\\s*[:=_-]?\\s*[0-9]+"
      )
      
      sra_accession <- str_extract(
        all_text,
        "\\b(SRR|SRX|SRS|SRP)[0-9]+\\b"
      )
      
      tibble(
        dataset = dataset,
        study_role = study_role,
        assay_expected = assay_expected,
        sample_id = sample_id,
        sample_title = title,
        organism = organism,
        source = source,
        strain_inferred = strain_inferred,
        treatment_group = treatment_group,
        concentration_inferred = concentration_inferred,
        exposure_time_inferred = exposure_time_inferred,
        replicate_inferred = replicate_inferred,
        molecule = molecule,
        library_strategy = library_strategy,
        library_source = library_source,
        library_selection = library_selection,
        platform = platform,
        characteristics = characteristics,
        treatment_protocol = treatment_protocol,
        growth_protocol = growth_protocol,
        extract_protocol = extract_protocol,
        description = description,
        sra_accession = sra_accession,
        relations = relation_values,
        supplementary_files = supplementary_files,
        combined_metadata = all_text
      )
    }
  )
}

# ------------------------------------------------------------
# 9. Download and parse all datasets
# ------------------------------------------------------------

all_metadata <- pmap_dfr(
  dataset_table,
  function(dataset, assay_expected, study_role) {
    
    xml_file <- download_miniml(dataset)
    
    Sys.sleep(1)
    
    parse_miniml_samples(
      xml_file = xml_file,
      dataset = dataset,
      assay_expected = assay_expected,
      study_role = study_role
    )
  }
)

# ------------------------------------------------------------
# 10. Verify retrieval
# ------------------------------------------------------------

if (nrow(all_metadata) == 0) {
  
  stop(
    "No sample metadata were retrieved. ",
    "Check the download messages shown above in the Console."
  )
}

# ------------------------------------------------------------
# 11. Add review flags
# ------------------------------------------------------------

all_metadata <- all_metadata |>
  mutate(
    needs_treatment_review = if_else(
      treatment_group == "Unclear",
      "Yes",
      "No"
    ),
    
    needs_strain_review = if_else(
      is.na(strain_inferred) |
        strain_inferred == "",
      "Yes",
      "No"
    ),
    
    needs_replicate_review = if_else(
      is.na(replicate_inferred) |
        replicate_inferred == "",
      "Yes",
      "No"
    ),
    
    metadata_status = if_else(
      needs_treatment_review == "No" &
        needs_strain_review == "No" &
        needs_replicate_review == "No",
      "Complete",
      "Needs review"
    )
  ) |>
  arrange(
    dataset,
    treatment_group,
    sample_id
  )

# ------------------------------------------------------------
# 12. Dataset-level summary
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
    
    sra_accessions_detected = sum(
      !is.na(sra_accession)
    ),
    
    samples_needing_review = sum(
      metadata_status == "Needs review"
    ),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 13. Experimental-design summary
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 14. Manual-review table
# ------------------------------------------------------------

manual_review <- all_metadata |>
  filter(metadata_status == "Needs review")

# ------------------------------------------------------------
# 15. SRA accession table
# ------------------------------------------------------------

sra_table <- all_metadata |>
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
    treatment_group,
    sra_accession,
    relations
  )

# ------------------------------------------------------------
# 16. Save outputs
# ------------------------------------------------------------

write_csv(
  all_metadata,
  file.path(
    output_folder,
    "BacRegRNA_MINiML_all_sample_metadata.csv"
  )
)

write_csv(
  dataset_summary,
  file.path(
    output_folder,
    "BacRegRNA_MINiML_dataset_summary.csv"
  )
)

write_csv(
  design_summary,
  file.path(
    output_folder,
    "BacRegRNA_MINiML_design_summary.csv"
  )
)

write_csv(
  manual_review,
  file.path(
    output_folder,
    "BacRegRNA_MINiML_manual_review.csv"
  )
)

write_csv(
  sra_table,
  file.path(
    output_folder,
    "BacRegRNA_MINiML_SRA_accessions.csv"
  )
)

# ------------------------------------------------------------
# 17. Final report
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 04C COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat("Total samples retrieved:", nrow(all_metadata), "\n")

cat(
  "Datasets retrieved:",
  paste(unique(all_metadata$dataset), collapse = ", "),
  "\n"
)

cat(
  "Samples needing manual review:",
  nrow(manual_review),
  "\n\n"
)

print(dataset_summary)

cat("\nExperimental design:\n")
print(design_summary)

cat(
  "\nFiles saved in:\n",
  output_folder,
  "\n"
)