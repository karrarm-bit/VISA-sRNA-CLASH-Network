# ============================================================
# BacRegRNA Project
# Script 02B: Repair missing NCBI metadata
# ============================================================

project_folder <- "D:/Bac-sRNA"
input_file <- file.path(
  project_folder,
  "02_screening",
  "Staphylococcus_aureus_probable_sRNA_records.csv"
)

output_file <- file.path(
  project_folder,
  "02_screening",
  "Staphylococcus_aureus_metadata_repaired.csv"
)

required_packages <- c(
  "rentrez", "dplyr", "purrr",
  "stringr", "readr", "tibble"
)

new_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

library(rentrez)
library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(tibble)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

records <- read_csv(input_file, show_col_types = FALSE)

# Safely extract any field from an NCBI summary
safe_field <- function(x, possible_names) {
  
  for (nm in possible_names) {
    if (!is.null(x[[nm]]) && length(x[[nm]]) > 0) {
      value <- paste(as.character(x[[nm]]), collapse = " | ")
      
      if (!is.na(value) && nzchar(value)) {
        return(value)
      }
    }
  }
  
  NA_character_
}

# Retrieve one record at a time
fetch_one_record <- function(database, record_id) {
  
  Sys.sleep(0.45)
  
  result <- tryCatch(
    entrez_summary(
      db = database,
      id = as.character(record_id),
      always_return_list = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(result)) {
    return(
      tibble(
        database = database,
        record_id = as.character(record_id),
        title_new = NA_character_,
        accession_new = NA_character_,
        organism_new = NA_character_,
        date_new = NA_character_,
        metadata_text = NA_character_,
        retrieval_status = "Failed"
      )
    )
  }
  
  metadata_text <- paste(
    unlist(result, recursive = TRUE, use.names = TRUE),
    collapse = " | "
  )
  
  # Search the complete metadata text for common accessions
  accession_from_text <- str_extract(
    metadata_text,
    paste0(
      "\\b(",
      "GSE[0-9]+|",
      "GSM[0-9]+|",
      "SRP[0-9]+|",
      "SRR[0-9]+|",
      "SRS[0-9]+|",
      "SRX[0-9]+|",
      "PRJNA[0-9]+|",
      "PRJEB[0-9]+",
      ")\\b"
    )
  )
  
  title_value <- safe_field(
    result,
    c(
      "title",
      "datasettitle",
      "studytitle",
      "caption"
    )
  )
  
  accession_value <- safe_field(
    result,
    c(
      "accession",
      "acc",
      "gse",
      "caption"
    )
  )
  
  if (is.na(accession_value)) {
    accession_value <- accession_from_text
  }
  
  organism_value <- safe_field(
    result,
    c(
      "organism",
      "organismname",
      "taxon",
      "taxname"
    )
  )
  
  date_value <- safe_field(
    result,
    c(
      "pubdate",
      "pdat",
      "createdate",
      "updatedate",
      "submissiondate"
    )
  )
  
  tibble(
    database = database,
    record_id = as.character(record_id),
    title_new = title_value,
    accession_new = accession_value,
    organism_new = organism_value,
    date_new = date_value,
    metadata_text = metadata_text,
    retrieval_status = "Retrieved"
  )
}

cat("Retrieving metadata for", nrow(records), "records...\n")

metadata_fixed <- map2_dfr(
  records$database,
  records$record_id,
  function(db, id) {
    
    cat("Retrieving:", db, id, "\n")
    
    fetch_one_record(
      database = db,
      record_id = id
    )
  }
)

final_records <- records |>
  mutate(record_id = as.character(record_id)) |>
  select(
    -any_of(c(
      "title",
      "accession",
      "organism",
      "publication_date",
      "summary_text"
    ))
  ) |>
  left_join(
    metadata_fixed,
    by = c("database", "record_id")
  ) |>
  mutate(
    combined_text = str_to_lower(
      paste(
        title_new,
        metadata_text,
        srna_terms,
        resistance_terms
      )
    ),
    
    srna_relevance = case_when(
      str_detect(
        combined_text,
        paste0(
          "small regulatory rna|",
          "small non.?coding rna|",
          "srna-seq|",
          "small rna sequencing|",
          "regulatory srna"
        )
      ) ~ "High",
      
      str_detect(
        combined_text,
        "\\bsrna\\b|small rna"
      ) ~ "Possible",
      
      TRUE ~ "Low"
    ),
    
    amr_relevance = case_when(
      str_detect(
        combined_text,
        paste0(
          "antimicrobial resistance|",
          "antibiotic resistance|",
          "multidrug resistant|",
          "drug resistance|",
          "methicillin resistant|",
          "\\bmrsa\\b"
        )
      ) ~ "Yes",
      
      TRUE ~ "Unclear or no"
    )
  ) |>
  select(
    database,
    record_id,
    accession_new,
    title_new,
    organism_new,
    date_new,
    srna_relevance,
    amr_relevance,
    retrieval_status,
    times_retrieved,
    resistance_category,
    srna_terms,
    resistance_terms,
    metadata_text
  ) |>
  arrange(
    desc(srna_relevance),
    desc(amr_relevance),
    desc(times_retrieved)
  )

write_csv(final_records, output_file)

cat("\n====================================\n")
cat("Metadata repair completed\n")
cat("====================================\n")
cat("Total records:", nrow(final_records), "\n")
cat(
  "Successfully retrieved:",
  sum(final_records$retrieval_status == "Retrieved"),
  "\n"
)
cat(
  "Failed:",
  sum(final_records$retrieval_status == "Failed"),
  "\n"
)
cat("Saved to:\n", output_file, "\n")

print(
  final_records |>
    count(
      retrieval_status,
      srna_relevance,
      amr_relevance
    )
)