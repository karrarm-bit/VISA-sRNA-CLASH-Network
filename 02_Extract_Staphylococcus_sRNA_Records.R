# ============================================================
# BacRegRNA Project
# Script 02: Extract unique Staphylococcus aureus records
# ============================================================

# 1. Set project folder
project_folder <- "D:/Bac-sRNA"
setwd(project_folder)

# 2. Install required packages
required_packages <- c(
  "rentrez",
  "dplyr",
  "purrr",
  "stringr",
  "tibble",
  "readr",
  "tidyr"
)

new_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

# 3. Load packages
library(rentrez)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(readr)
library(tidyr)

# 4. Read the complete search-results file
input_file <- file.path(
  project_folder,
  "BacRegRNA_initial_search_results.csv"
)

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    input_file,
    "\nMake sure BacRegRNA_initial_search_results.csv is inside D:/Bac-sRNA"
  )
}

search_results <- read_csv(
  input_file,
  show_col_types = FALSE
)

# 5. Filter Staphylococcus aureus results
staph_results <- search_results |>
  filter(
    bacterium == "Staphylococcus aureus",
    !is.na(record_ids),
    record_ids != "",
    number_found > 0
  )

cat(
  "\nNumber of search rows for Staphylococcus aureus:",
  nrow(staph_results),
  "\n"
)

# 6. Split record IDs into individual rows
staph_record_ids <- staph_results |>
  select(
    database,
    srna_term,
    resistance_term,
    query,
    record_ids
  ) |>
  separate_rows(
    record_ids,
    sep = ";"
  ) |>
  mutate(
    record_id = str_trim(record_ids)
  ) |>
  filter(
    !is.na(record_id),
    record_id != ""
  ) |>
  select(
    database,
    record_id,
    srna_term,
    resistance_term,
    query
  )

# 7. Summarize how many times each record appeared
record_query_summary <- staph_record_ids |>
  group_by(
    database,
    record_id
  ) |>
  summarise(
    times_retrieved = n(),
    srna_terms = paste(
      sort(unique(srna_term)),
      collapse = "; "
    ),
    resistance_terms = paste(
      sort(unique(resistance_term)),
      collapse = "; "
    ),
    .groups = "drop"
  ) |>
  arrange(
    database,
    desc(times_retrieved)
  )

cat(
  "Unique records before metadata retrieval:",
  nrow(record_query_summary),
  "\n"
)

# 8. Function to retrieve metadata in batches
fetch_summaries_batch <- function(
    ids,
    database,
    batch_size = 100
) {
  
  if (length(ids) == 0) {
    return(tibble())
  }
  
  id_batches <- split(
    ids,
    ceiling(seq_along(ids) / batch_size)
  )
  
  map_dfr(
    seq_along(id_batches),
    function(i) {
      
      current_ids <- id_batches[[i]]
      
      cat(
        "Fetching",
        database,
        "batch",
        i,
        "of",
        length(id_batches),
        "with",
        length(current_ids),
        "records...\n"
      )
      
      Sys.sleep(0.5)
      
      summaries <- tryCatch(
        entrez_summary(
          db = database,
          id = current_ids,
          always_return_list = TRUE
        ),
        error = function(e) {
          message(
            "Batch failed: ",
            conditionMessage(e)
          )
          return(NULL)
        }
      )
      
      if (is.null(summaries)) {
        return(
          tibble(
            database = database,
            record_id = current_ids,
            title = NA_character_,
            accession = NA_character_,
            organism = NA_character_,
            publication_date = NA_character_,
            summary_text = NA_character_
          )
        )
      }
      
      map_dfr(
        seq_along(summaries),
        function(j) {
          
          x <- summaries[[j]]
          
          tibble(
            database = database,
            record_id = as.character(current_ids[j]),
            
            title = dplyr::coalesce(
              as.character(x$title),
              NA_character_
            ),
            
            accession = dplyr::coalesce(
              as.character(x$accession),
              as.character(x$acc),
              as.character(x$gse),
              NA_character_
            ),
            
            organism = dplyr::coalesce(
              as.character(x$organism),
              as.character(x$organismname),
              as.character(x$taxon),
              NA_character_
            ),
            
            publication_date = dplyr::coalesce(
              as.character(x$pubdate),
              as.character(x$pdat),
              as.character(x$updatedate),
              NA_character_
            ),
            
            summary_text = paste(
              unlist(x),
              collapse = " | "
            )
          )
        }
      )
    }
  )
}

# 9. Retrieve SRA metadata
sra_ids <- record_query_summary |>
  filter(database == "sra") |>
  pull(record_id)

sra_metadata <- fetch_summaries_batch(
  ids = sra_ids,
  database = "sra",
  batch_size = 100
)

# 10. Retrieve GEO metadata
geo_ids <- record_query_summary |>
  filter(database == "gds") |>
  pull(record_id)

geo_metadata <- fetch_summaries_batch(
  ids = geo_ids,
  database = "gds",
  batch_size = 100
)

# 11. Combine metadata
all_metadata <- bind_rows(
  sra_metadata,
  geo_metadata
)

# 12. Join metadata with search-term information
staph_unique_records <- record_query_summary |>
  left_join(
    all_metadata,
    by = c(
      "database",
      "record_id"
    )
  ) |>
  mutate(
    combined_text = str_to_lower(
      paste(
        title,
        summary_text,
        srna_terms,
        resistance_terms
      )
    ),
    
    probable_srna = case_when(
      str_detect(
        combined_text,
        "small regulatory rna|srna-seq|small rna sequencing|small non-coding rna|small noncoding rna|regulatory rna"
      ) ~ "Yes",
      
      str_detect(
        combined_text,
        "\\bsrna\\b|small rna"
      ) ~ "Possible",
      
      TRUE ~ "No"
    ),
    
    resistance_category = case_when(
      str_detect(
        combined_text,
        "antibiotic resistance|antimicrobial resistance|multidrug resistance|drug resistance"
      ) ~ "Antimicrobial resistance",
      
      str_detect(
        combined_text,
        "biofilm"
      ) ~ "Biofilm",
      
      str_detect(
        combined_text,
        "virulence|pathogenicity"
      ) ~ "Virulence",
      
      str_detect(
        combined_text,
        "persister|tolerance"
      ) ~ "Persistence or tolerance",
      
      str_detect(
        combined_text,
        "stress response|oxidative stress|acid stress|heat shock"
      ) ~ "Stress response",
      
      TRUE ~ "Other or unclear"
    )
  ) |>
  select(
    database,
    record_id,
    accession,
    title,
    organism,
    publication_date,
    probable_srna,
    resistance_category,
    times_retrieved,
    srna_terms,
    resistance_terms,
    summary_text
  ) |>
  arrange(
    desc(probable_srna),
    desc(times_retrieved)
  )

# 13. Create output folder
output_folder <- file.path(
  project_folder,
  "02_screening"
)

if (!dir.exists(output_folder)) {
  dir.create(
    output_folder,
    recursive = TRUE
  )
}

# 14. Save all unique records
write_csv(
  staph_unique_records,
  file.path(
    output_folder,
    "Staphylococcus_aureus_unique_records.csv"
  )
)

# 15. Save probable sRNA records only
probable_srna_records <- staph_unique_records |>
  filter(
    probable_srna %in% c(
      "Yes",
      "Possible"
    )
  )

write_csv(
  probable_srna_records,
  file.path(
    output_folder,
    "Staphylococcus_aureus_probable_sRNA_records.csv"
  )
)

# 16. Save category summary
category_summary <- staph_unique_records |>
  count(
    database,
    probable_srna,
    resistance_category,
    name = "number_of_records"
  ) |>
  arrange(
    database,
    desc(number_of_records)
  )

write_csv(
  category_summary,
  file.path(
    output_folder,
    "Staphylococcus_aureus_record_summary.csv"
  )
)

# 17. Print final summary
cat("\n============================================\n")
cat("Step 2 completed successfully\n")
cat("============================================\n")

cat(
  "Total unique records:",
  nrow(staph_unique_records),
  "\n"
)

cat(
  "Probable or possible sRNA records:",
  nrow(probable_srna_records),
  "\n"
)

cat(
  "Files saved inside:\n",
  output_folder,
  "\n\n"
)

print(category_summary)