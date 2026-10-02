# ============================================================
# BacRegRNA Project
# Script 07: Discover and download processed GEO files
#
# Inputs:
#   Local GEO MINiML XML files from Step 04E
#
# Outputs:
#   07_processed_data/
#     ├── processed_file_catalog.csv
#     ├── processed_download_report.csv
#     ├── selected_processed_files.csv
#     └── downloaded_files/
# ============================================================

# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

xml_folder <- file.path(
  project_folder,
  "04_metadata",
  "extracted"
)

output_folder <- file.path(
  project_folder,
  "07_processed_data"
)

download_folder <- file.path(
  output_folder,
  "downloaded_files"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  download_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

setwd(project_folder)

options(timeout = 1800)

# ------------------------------------------------------------
# 2. Install and load required packages
# ------------------------------------------------------------

required_packages <- c(
  "xml2",
  "curl",
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

library(xml2)
library(curl)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(tibble)
library(readr)

# ------------------------------------------------------------
# 3. Define datasets
# ------------------------------------------------------------

study_table <- tribble(
  ~dataset,    ~assay_expected,    ~study_role,
  "GSE254530", "RNA-seq",          "Discovery",
  "GSE254531", "Ribo-seq",         "Discovery support",
  "GSE254532", "RNase III-CLASH",  "Discovery",
  "GSE158830", "Mixed assays",      "Validation"
)

# ------------------------------------------------------------
# 4. Locate XML files
# ------------------------------------------------------------

xml_files <- list.files(
  xml_folder,
  pattern = "_family\\.xml$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(xml_files) == 0) {
  stop(
    "No GEO XML files were found inside:\n",
    xml_folder
  )
}

xml_table <- tibble(
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
print(xml_table)

# ------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------

clean_url <- function(x) {
  
  if (is.na(x) || x == "") {
    return(NA_character_)
  }
  
  x <- str_squish(x)
  
  # Convert FTP links to HTTPS
  x <- str_replace(
    x,
    "^ftp://ftp\\.ncbi\\.nlm\\.nih\\.gov",
    "https://ftp.ncbi.nlm.nih.gov"
  )
  
  x <- str_replace(
    x,
    "^ftp://",
    "https://"
  )
  
  x
}

safe_filename <- function(url, fallback_name) {
  
  if (is.na(url) || url == "") {
    return(fallback_name)
  }
  
  name <- basename(
    sub("\\?.*$", "", url)
  )
  
  name <- URLdecode(name)
  
  name <- str_replace_all(
    name,
    "[<>:\"/\\\\|?*]",
    "_"
  )
  
  if (is.na(name) || name == "") {
    name <- fallback_name
  }
  
  name
}

classify_processed_file <- function(
    file_name,
    supplementary_type,
    url
) {
  
  text <- str_to_lower(
    paste(
      file_name,
      supplementary_type,
      url
    )
  )
  
  case_when(
    
    str_detect(
      text,
      "fastq|\\.fq(\\.gz)?$|\\.sra$"
    ) ~ "Raw sequencing data",
    
    str_detect(
      text,
      "count|counts|readcount|read_count|featurecounts"
    ) ~ "Count matrix or count table",
    
    str_detect(
      text,
      "deseq|differential|deg|fold.change|log2fc|results"
    ) ~ "Differential-expression results",
    
    str_detect(
      text,
      "clash|hybrid|interaction|chimera|duplex"
    ) ~ "CLASH interaction data",
    
    str_detect(
      text,
      "ribo|ribosome|footprint|rpf"
    ) ~ "Ribo-seq processed data",
    
    str_detect(
      text,
      "bed|bedgraph|bigwig|\\.bw$|\\.wig"
    ) ~ "Genome-track or interval data",
    
    str_detect(
      text,
      "annotation|gff|gtf|gene.?info|feature"
    ) ~ "Annotation file",
    
    str_detect(
      text,
      "\\.csv|\\.tsv|\\.txt|\\.xlsx|\\.xls|\\.gz|\\.zip"
    ) ~ "Other processed table or archive",
    
    TRUE ~ "Unclassified"
  )
}

# ------------------------------------------------------------
# 6. Parse Series-level and Sample-level supplementary links
# ------------------------------------------------------------

extract_supplementary_from_xml <- function(
    xml_file,
    dataset,
    assay_expected,
    study_role
) {
  
  cat("\nParsing supplementary links for:", dataset, "\n")
  
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
  
  # ----------------------------------------------------------
  # Series-level supplementary files
  # ----------------------------------------------------------
  
  series_nodes <- xml_find_all(
    document,
    ".//*[local-name()='Series']/*[local-name()='Supplementary-Data']"
  )
  
  series_table <- if (length(series_nodes) == 0) {
    
    tibble()
    
  } else {
    
    map_dfr(
      seq_along(series_nodes),
      function(i) {
        
        node <- series_nodes[[i]]
        
        url <- xml_text(
          node,
          trim = TRUE
        )
        
        supplementary_type <- xml_attr(
          node,
          "type"
        )
        
        url <- clean_url(url)
        
        file_name <- safe_filename(
          url,
          paste0(
            dataset,
            "_series_file_",
            i
          )
        )
        
        tibble(
          dataset = dataset,
          study_role = study_role,
          assay_expected = assay_expected,
          source_level = "Series",
          sample_id = NA_character_,
          gsm_accession = NA_character_,
          supplementary_type = supplementary_type,
          url = url,
          file_name = file_name
        )
      }
    )
  }
  
  # ----------------------------------------------------------
  # Sample-level supplementary files
  # ----------------------------------------------------------
  
  sample_nodes <- xml_find_all(
    document,
    ".//*[local-name()='Sample']"
  )
  
  sample_table <- if (length(sample_nodes) == 0) {
    
    tibble()
    
  } else {
    
    map_dfr(
      sample_nodes,
      function(sample_node) {
        
        sample_id <- xml_attr(
          sample_node,
          "iid"
        )
        
        gsm_accession <- xml_text(
          xml_find_first(
            sample_node,
            "./*[local-name()='Accession']"
          ),
          trim = TRUE
        )
        
        supplementary_nodes <- xml_find_all(
          sample_node,
          ".//*[local-name()='Supplementary-Data']"
        )
        
        if (length(supplementary_nodes) == 0) {
          return(tibble())
        }
        
        map_dfr(
          seq_along(supplementary_nodes),
          function(i) {
            
            node <- supplementary_nodes[[i]]
            
            url <- xml_text(
              node,
              trim = TRUE
            )
            
            supplementary_type <- xml_attr(
              node,
              "type"
            )
            
            url <- clean_url(url)
            
            file_name <- safe_filename(
              url,
              paste0(
                dataset,
                "_",
                gsm_accession,
                "_sample_file_",
                i
              )
            )
            
            tibble(
              dataset = dataset,
              study_role = study_role,
              assay_expected = assay_expected,
              source_level = "Sample",
              sample_id = sample_id,
              gsm_accession = gsm_accession,
              supplementary_type = supplementary_type,
              url = url,
              file_name = file_name
            )
          }
        )
      }
    )
  }
  
  bind_rows(
    series_table,
    sample_table
  ) |>
    filter(
      !is.na(url),
      url != ""
    ) |>
    mutate(
      file_category = classify_processed_file(
        file_name,
        supplementary_type,
        url
      )
    )
}

# ------------------------------------------------------------
# 7. Build complete processed-file catalog
# ------------------------------------------------------------

processed_catalog <- pmap_dfr(
  xml_table,
  function(
    xml_file,
    dataset,
    assay_expected,
    study_role
  ) {
    
    extract_supplementary_from_xml(
      xml_file = xml_file,
      dataset = dataset,
      assay_expected = assay_expected,
      study_role = study_role
    )
  }
) |>
  distinct(
    dataset,
    url,
    .keep_all = TRUE
  ) |>
  arrange(
    dataset,
    source_level,
    file_category,
    file_name
  )

if (nrow(processed_catalog) == 0) {
  stop(
    "No supplementary-file links were found in the XML files."
  )
}

write_csv(
  processed_catalog,
  file.path(
    output_folder,
    "processed_file_catalog.csv"
  )
)

cat(
  "\nSupplementary links detected:",
  nrow(processed_catalog),
  "\n"
)

print(
  processed_catalog |>
    count(
      dataset,
      file_category,
      name = "number_of_files"
    ),
  n = Inf
)

# ------------------------------------------------------------
# 8. Select files relevant to the current project
# ------------------------------------------------------------

selected_files <- processed_catalog |>
  filter(
    file_category != "Raw sequencing data"
  ) |>
  mutate(
    priority = case_when(
      
      dataset == "GSE254530" &
        file_category %in% c(
          "Count matrix or count table",
          "Differential-expression results",
          "Other processed table or archive"
        ) ~ "Priority 1",
      
      dataset == "GSE254532" &
        file_category %in% c(
          "CLASH interaction data",
          "Count matrix or count table",
          "Other processed table or archive"
        ) ~ "Priority 1",
      
      dataset == "GSE158830" &
        file_category %in% c(
          "Count matrix or count table",
          "Differential-expression results",
          "CLASH interaction data",
          "Annotation file",
          "Other processed table or archive"
        ) ~ "Priority 1",
      
      dataset == "GSE254531" ~ "Priority 2",
      
      TRUE ~ "Priority 3"
    )
  ) |>
  arrange(
    priority,
    dataset,
    file_category,
    file_name
  )

write_csv(
  selected_files,
  file.path(
    output_folder,
    "selected_processed_files.csv"
  )
)

# ------------------------------------------------------------
# 9. Download function using curl
# ------------------------------------------------------------

download_processed_file <- function(
    dataset,
    url,
    file_name,
    file_category,
    priority
) {
  
  dataset_folder <- file.path(
    download_folder,
    dataset
  )
  
  dir.create(
    dataset_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  destination <- file.path(
    dataset_folder,
    file_name
  )
  
  cat("\n--------------------------------------------\n")
  cat("Dataset:", dataset, "\n")
  cat("Category:", file_category, "\n")
  cat("Priority:", priority, "\n")
  cat("File:", file_name, "\n")
  cat("--------------------------------------------\n")
  
  if (
    file.exists(destination) &&
    !is.na(file.info(destination)$size) &&
    file.info(destination)$size > 0
  ) {
    
    cat("Existing file detected; download skipped.\n")
    
    return(
      tibble(
        dataset = dataset,
        file_name = file_name,
        file_category = file_category,
        priority = priority,
        original_url = url,
        local_file = destination,
        status = "Already exists",
        size_bytes = file.info(destination)$size
      )
    )
  }
  
  handle <- new_handle()
  
  handle_setopt(
    handle,
    useragent = paste0(
      "Mozilla/5.0 (Windows NT 10.0; Win64; x64) ",
      "AppleWebKit/537.36 Chrome/150 Safari/537.36"
    ),
    followlocation = TRUE,
    timeout = 1200,
    connecttimeout = 120
  )
  
  handle_setheaders(
    handle,
    Accept = "*/*",
    `Accept-Language` = "en-US,en;q=0.9"
  )
  
  status <- tryCatch(
    {
      
      curl_download(
        url = url,
        destfile = destination,
        quiet = FALSE,
        mode = "wb",
        handle = handle
      )
      
      "Downloaded"
    },
    error = function(e) {
      
      message(
        "Download failed: ",
        conditionMessage(e)
      )
      
      "Failed"
    }
  )
  
  size_bytes <- if (
    file.exists(destination)
  ) {
    file.info(destination)$size
  } else {
    0
  }
  
  if (
    is.na(size_bytes) ||
    size_bytes == 0
  ) {
    
    status <- "Failed or empty"
    
    if (file.exists(destination)) {
      file.remove(destination)
    }
    
    size_bytes <- 0
  }
  
  tibble(
    dataset = dataset,
    file_name = file_name,
    file_category = file_category,
    priority = priority,
    original_url = url,
    local_file = destination,
    status = status,
    size_bytes = size_bytes
  )
}

# ------------------------------------------------------------
# 10. Download selected processed files
# ------------------------------------------------------------

download_report <- pmap_dfr(
  selected_files |>
    select(
      dataset,
      url,
      file_name,
      file_category,
      priority
    ),
  download_processed_file
)

write_csv(
  download_report,
  file.path(
    output_folder,
    "processed_download_report.csv"
  )
)

# ------------------------------------------------------------
# 11. Inspect downloaded file signatures
# ------------------------------------------------------------

inspect_file <- function(file_path) {
  
  if (
    is.na(file_path) ||
    !file.exists(file_path) ||
    file.info(file_path)$size == 0
  ) {
    return(
      tibble(
        detected_format = "Missing",
        readable_as_text = FALSE
      )
    )
  }
  
  extension <- str_to_lower(
    tools::file_ext(file_path)
  )
  
  detected_format <- case_when(
    extension == "gz" ~ "gzip",
    extension == "zip" ~ "zip",
    extension %in% c("csv", "tsv", "txt") ~ extension,
    extension %in% c("xlsx", "xls") ~ "Excel",
    extension %in% c("bed", "bedgraph", "wig", "bw") ~ extension,
    TRUE ~ paste0("Extension: ", extension)
  )
  
  readable_as_text <- extension %in% c(
    "csv",
    "tsv",
    "txt",
    "bed",
    "bedgraph",
    "wig"
  )
  
  tibble(
    detected_format = detected_format,
    readable_as_text = readable_as_text
  )
}

inspection_table <- download_report |>
  filter(
    status %in% c(
      "Downloaded",
      "Already exists"
    ),
    size_bytes > 0
  ) |>
  mutate(
    inspection = map(
      local_file,
      inspect_file
    )
  ) |>
  unnest(
    inspection
  ) |>
  arrange(
    dataset,
    file_category,
    file_name
  )

write_csv(
  inspection_table,
  file.path(
    output_folder,
    "downloaded_file_inspection.csv"
  )
)

# ------------------------------------------------------------
# 12. Create dataset-level download summary
# ------------------------------------------------------------

download_summary <- download_report |>
  group_by(
    dataset,
    file_category,
    priority
  ) |>
  summarise(
    files_attempted = n(),
    downloaded_or_existing = sum(
      status %in% c(
        "Downloaded",
        "Already exists"
      )
    ),
    failed_files = sum(
      !status %in% c(
        "Downloaded",
        "Already exists"
      )
    ),
    total_size_MB = round(
      sum(
        size_bytes,
        na.rm = TRUE
      ) / 1024^2,
      3
    ),
    .groups = "drop"
  ) |>
  arrange(
    priority,
    dataset,
    file_category
  )

write_csv(
  download_summary,
  file.path(
    output_folder,
    "processed_download_summary.csv"
  )
)

# ------------------------------------------------------------
# 13. Final report
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 07 COMPLETED\n")
cat("============================================\n")

cat(
  "Supplementary links found:",
  nrow(processed_catalog),
  "\n"
)

cat(
  "Processed files selected:",
  nrow(selected_files),
  "\n"
)

cat(
  "Successfully downloaded or already present:",
  sum(
    download_report$status %in% c(
      "Downloaded",
      "Already exists"
    )
  ),
  "\n"
)

cat(
  "Failed:",
  sum(
    !download_report$status %in% c(
      "Downloaded",
      "Already exists"
    )
  ),
  "\n\n"
)

print(
  download_summary,
  n = Inf
)

cat(
  "\nFiles and reports saved inside:\n",
  output_folder,
  "\n"
)