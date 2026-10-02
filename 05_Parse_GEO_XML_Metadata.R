# ============================================================
# BacRegRNA Project
# Script 05: Parse local GEO MINiML XML files
# ============================================================

# ------------------------------------------------------------
# 1. Project folders
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

xml_folder <- file.path(
  project_folder,
  "04_metadata",
  "extracted"
)

output_folder <- file.path(
  project_folder,
  "05_parsed_metadata"
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
# 4. Define study information
# ------------------------------------------------------------

study_table <- tribble(
  ~dataset,    ~assay_expected,    ~study_role,
  "GSE254530", "RNA-seq",          "Discovery",
  "GSE254531", "Ribo-seq",         "Discovery",
  "GSE254532", "RNase III-CLASH",  "Discovery",
  "GSE158830", "RNase III-CLASH",  "Validation"
)

# ------------------------------------------------------------
# 5. Locate all XML files
# ------------------------------------------------------------

xml_files <- list.files(
  xml_folder,
  pattern = "_family\\.xml$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(xml_files) == 0) {
  stop(
    "No XML files were found inside: ",
    xml_folder
  )
}

xml_file_table <- tibble(
  xml_file = xml_files,
  dataset = str_extract(
    basename(xml_files),
    "GSE[0-9]+"
  )
) |>
  left_join(
    study_table,
    by = "dataset"
  )

cat("\nXML files detected:\n")
print(xml_file_table)

# ------------------------------------------------------------
# 6. Helper functions
# ------------------------------------------------------------

safe_text <- function(node, xpath) {
  
  result <- xml_find_first(node, xpath)
  
  if (inherits(result, "xml_missing")) {
    return(NA_character_)
  }
  
  value <- xml_text(result, trim = TRUE)
  
  if (
    length(value) == 0 ||
    is.na(value) ||
    value == ""
  ) {
    return(NA_character_)
  }
  
  str_squish(value)
}

safe_multiple_text <- function(node, xpath) {
  
  results <- xml_find_all(node, xpath)
  
  if (length(results) == 0) {
    return(NA_character_)
  }
  
  values <- xml_text(results, trim = TRUE)
  
  values <- str_squish(values)
  
  values <- values[
    !is.na(values) &
      values != ""
  ]
  
  if (length(values) == 0) {
    return(NA_character_)
  }
  
  paste(
    unique(values),
    collapse = " | "
  )
}

extract_characteristics <- function(sample_node) {
  
  nodes <- xml_find_all(
    sample_node,
    ".//*[local-name()='Characteristics']"
  )
  
  if (length(nodes) == 0) {
    return(NA_character_)
  }
  
  values <- map_chr(
    nodes,
    function(node) {
      
      tag <- xml_attr(node, "tag")
      value <- xml_text(node, trim = TRUE)
      
      value <- str_squish(value)
      
      if (
        !is.na(tag) &&
        tag != ""
      ) {
        paste0(tag, ": ", value)
      } else {
        value
      }
    }
  )
  
  paste(
    unique(values),
    collapse = " | "
  )
}

extract_relations <- function(sample_node) {
  
  nodes <- xml_find_all(
    sample_node,
    ".//*[local-name()='Relation']"
  )
  
  if (length(nodes) == 0) {
    return(NA_character_)
  }
  
  values <- map_chr(
    nodes,
    function(node) {
      
      relation_type <- xml_attr(node, "type")
      relation_target <- xml_attr(node, "target")
      relation_text <- xml_text(node, trim = TRUE)
      
      pieces <- c(
        relation_type,
        relation_target,
        relation_text
      )
      
      pieces <- pieces[
        !is.na(pieces) &
          pieces != ""
      ]
      
      paste(
        pieces,
        collapse = ": "
      )
    }
  )
  
  paste(
    unique(values),
    collapse = " | "
  )
}

extract_supplementary <- function(sample_node) {
  
  nodes <- xml_find_all(
    sample_node,
    ".//*[local-name()='Supplementary-Data']"
  )
  
  if (length(nodes) == 0) {
    return(NA_character_)
  }
  
  values <- map_chr(
    nodes,
    function(node) {
      
      supplementary_type <- xml_attr(node, "type")
      supplementary_text <- xml_text(node, trim = TRUE)
      
      pieces <- c(
        supplementary_type,
        supplementary_text
      )
      
      pieces <- pieces[
        !is.na(pieces) &
          pieces != ""
      ]
      
      paste(
        pieces,
        collapse = ": "
      )
    }
  )
  
  paste(
    unique(values),
    collapse = " | "
  )
}

# ------------------------------------------------------------
# 7. Parse one GEO XML file
# ------------------------------------------------------------

parse_one_xml <- function(
    xml_file,
    dataset,
    assay_expected,
    study_role
) {
  
  cat("\n============================================\n")
  cat("Parsing dataset:", dataset, "\n")
  cat("File:", xml_file, "\n")
  cat("============================================\n")
  
  document <- tryCatch(
    read_xml(xml_file),
    error = function(e) {
      message(
        "Could not read ",
        dataset,
        ": ",
        conditionMessage(e)
      )
      
      return(NULL)
    }
  )
  
  if (is.null(document)) {
    return(tibble())
  }
  
  sample_nodes <- xml_find_all(
    document,
    ".//*[local-name()='Sample']"
  )
  
  cat(
    "Number of Sample nodes:",
    length(sample_nodes),
    "\n"
  )
  
  if (length(sample_nodes) == 0) {
    return(tibble())
  }
  
  map_dfr(
    sample_nodes,
    function(sample_node) {
      
      sample_id <- xml_attr(
        sample_node,
        "iid"
      )
      
      if (
        is.na(sample_id) ||
        sample_id == ""
      ) {
        sample_id <- safe_text(
          sample_node,
          ".//*[local-name()='Accession']"
        )
      }
      
      sample_title <- safe_text(
        sample_node,
        "./*[local-name()='Title']"
      )
      
      accession <- safe_text(
        sample_node,
        "./*[local-name()='Accession']"
      )
      
      source_name <- safe_text(
        sample_node,
        ".//*[local-name()='Source']"
      )
      
      organism <- safe_text(
        sample_node,
        ".//*[local-name()='Organism']"
      )
      
      molecule <- safe_text(
        sample_node,
        ".//*[local-name()='Molecule']"
      )
      
      characteristics <- extract_characteristics(
        sample_node
      )
      
      treatment_protocol <- safe_multiple_text(
        sample_node,
        ".//*[local-name()='Treatment-Protocol']"
      )
      
      growth_protocol <- safe_multiple_text(
        sample_node,
        ".//*[local-name()='Growth-Protocol']"
      )
      
      extract_protocol <- safe_multiple_text(
        sample_node,
        ".//*[local-name()='Extract-Protocol']"
      )
      
      labeling_protocol <- safe_multiple_text(
        sample_node,
        ".//*[local-name()='Label-Protocol']"
      )
      
      hybridization_protocol <- safe_multiple_text(
        sample_node,
        ".//*[local-name()='Hybridization-Protocol']"
      )
      
      description <- safe_multiple_text(
        sample_node,
        ".//*[local-name()='Description']"
      )
      
      library_strategy <- safe_text(
        sample_node,
        ".//*[local-name()='Library-Strategy']"
      )
      
      library_source <- safe_text(
        sample_node,
        ".//*[local-name()='Library-Source']"
      )
      
      library_selection <- safe_text(
        sample_node,
        ".//*[local-name()='Library-Selection']"
      )
      
      platform_node <- xml_find_first(
        sample_node,
        ".//*[local-name()='Platform-Ref']"
      )
      
      platform <- if (
        inherits(platform_node, "xml_missing")
      ) {
        NA_character_
      } else {
        xml_attr(platform_node, "ref")
      }
      
      relations <- extract_relations(
        sample_node
      )
      
      supplementary_files <- extract_supplementary(
        sample_node
      )
      
      combined_metadata <- str_to_lower(
        str_squish(
          paste(
            sample_title,
            source_name,
            characteristics,
            treatment_protocol,
            growth_protocol,
            extract_protocol,
            description,
            library_strategy,
            library_source,
            library_selection,
            relations,
            supplementary_files
          )
        )
      )
      
      treatment_group <- case_when(
        
        str_detect(
          combined_metadata,
          paste0(
            "untreated|",
            "un-treated|",
            "without vancomycin|",
            "no vancomycin|",
            "non-treated|",
            "\\bcontrol\\b"
          )
        ) ~ "Control",
        
        str_detect(
          combined_metadata,
          "vancomycin|vanco"
        ) ~ "Vancomycin",
        
        TRUE ~ "Unclear"
      )
      
      strain_inferred <- case_when(
        
        str_detect(
          combined_metadata,
          "jkd6008"
        ) ~ "JKD6008",
        
        str_detect(
          combined_metadata,
          "jkd6009"
        ) ~ "JKD6009",
        
        str_detect(
          combined_metadata,
          "usa300"
        ) ~ "USA300",
        
        str_detect(
          combined_metadata,
          "newman"
        ) ~ "Newman",
        
        str_detect(
          combined_metadata,
          "hg003"
        ) ~ "HG003",
        
        TRUE ~ str_match(
          combined_metadata,
          "strain\\s*[:=]\\s*([^|;,]+)"
        )[, 2]
      )
      
      concentration_inferred <- str_extract(
        combined_metadata,
        paste0(
          "[0-9.]+\\s*",
          "(µg/ml|μg/ml|ug/ml|mg/l|mg/ml)"
        )
      )
      
      exposure_time_inferred <- str_extract(
        combined_metadata,
        paste0(
          "[0-9.]+\\s*",
          "(min|mins|minute|minutes|",
          "hour|hours|hr|hrs)"
        )
      )
      
      replicate_inferred <- str_extract(
        combined_metadata,
        paste0(
          "(biological\\s*)?",
          "(replicate|rep)\\s*",
          "[:=_-]?\\s*[0-9]+"
        )
      )
      
      sra_accession <- str_extract(
        combined_metadata,
        "\\b(SRR|SRX|SRS|SRP)[0-9]+\\b"
      )
      
      bioproject_accession <- str_extract(
        combined_metadata,
        "\\b(PRJNA|PRJEB)[0-9]+\\b"
      )
      
      tibble(
        dataset = dataset,
        study_role = study_role,
        assay_expected = assay_expected,
        sample_id = sample_id,
        gsm_accession = accession,
        sample_title = sample_title,
        organism = organism,
        source_name = source_name,
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
        labeling_protocol = labeling_protocol,
        hybridization_protocol = hybridization_protocol,
        description = description,
        sra_accession = sra_accession,
        bioproject_accession = bioproject_accession,
        relations = relations,
        supplementary_files = supplementary_files,
        combined_metadata = combined_metadata
      )
    }
  )
}

# ------------------------------------------------------------
# 8. Parse all four studies
# ------------------------------------------------------------

all_metadata <- pmap_dfr(
  xml_file_table,
  function(
    xml_file,
    dataset,
    assay_expected,
    study_role
  ) {
    
    parse_one_xml(
      xml_file = xml_file,
      dataset = dataset,
      assay_expected = assay_expected,
      study_role = study_role
    )
  }
)

# ------------------------------------------------------------
# 9. Check parsing result
# ------------------------------------------------------------

if (nrow(all_metadata) == 0) {
  stop(
    "The XML files were read, but no samples were extracted."
  )
}

# ------------------------------------------------------------
# 10. Add metadata-review flags
# ------------------------------------------------------------

all_metadata <- all_metadata |>
  mutate(
    
    treatment_review = if_else(
      treatment_group == "Unclear",
      "Needs review",
      "OK"
    ),
    
    strain_review = if_else(
      is.na(strain_inferred) |
        strain_inferred == "",
      "Needs review",
      "OK"
    ),
    
    replicate_review = if_else(
      is.na(replicate_inferred) |
        replicate_inferred == "",
      "Needs review",
      "OK"
    ),
    
    sra_review = if_else(
      is.na(sra_accession) |
        sra_accession == "",
      "Needs review",
      "OK"
    ),
    
    metadata_status = case_when(
      
      treatment_review == "OK" &
        strain_review == "OK" &
        replicate_review == "OK" ~ "Complete",
      
      TRUE ~ "Needs review"
    )
  ) |>
  arrange(
    dataset,
    treatment_group,
    sample_id
  )

# ------------------------------------------------------------
# 11. Build summary tables
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
      sort(
        unique(
          na.omit(strain_inferred)
        )
      ),
      collapse = "; "
    ),
    
    samples_with_SRA = sum(
      !is.na(sra_accession),
      na.rm = TRUE
    ),
    
    samples_needing_review = sum(
      metadata_status == "Needs review",
      na.rm = TRUE
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
  filter(
    metadata_status == "Needs review"
  )

sra_table <- all_metadata |>
  select(
    dataset,
    sample_id,
    gsm_accession,
    sample_title,
    assay_expected,
    treatment_group,
    strain_inferred,
    sra_accession,
    bioproject_accession,
    relations
  ) |>
  distinct()

# ------------------------------------------------------------
# 12. Save files
# ------------------------------------------------------------

write_csv(
  all_metadata,
  file.path(
    output_folder,
    "BacRegRNA_all_sample_metadata.csv"
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
    "BacRegRNA_design_summary.csv"
  )
)

write_csv(
  manual_review,
  file.path(
    output_folder,
    "BacRegRNA_metadata_manual_review.csv"
  )
)

write_csv(
  sra_table,
  file.path(
    output_folder,
    "BacRegRNA_sample_SRA_table.csv"
  )
)

# ------------------------------------------------------------
# 13. Final console report
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 05 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat(
  "Total samples extracted:",
  nrow(all_metadata),
  "\n"
)

cat(
  "Datasets extracted:",
  paste(
    unique(all_metadata$dataset),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Samples requiring review:",
  nrow(manual_review),
  "\n\n"
)

print(dataset_summary)

cat("\nExperimental design summary:\n")
print(design_summary)

cat(
  "\nFiles saved inside:\n",
  output_folder,
  "\n"
)