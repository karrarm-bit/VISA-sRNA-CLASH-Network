# ============================================================
# BacRegRNA Project
# Script 01: Search GEO and SRA for bacterial sRNA datasets
# ============================================================

# 1. Install required packages
required_packages <- c(
  "rentrez",
  "dplyr",
  "purrr",
  "stringr",
  "tibble",
  "readr"
)

new_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

# 2. Load packages
library(rentrez)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(readr)

# 3. Create project folders
folders <- c(
  "01_search_results",
  "02_screening",
  "03_selected_datasets",
  "04_raw_data",
  "05_results",
  "06_figures"
)

walk(folders, ~ dir.create(.x, showWarnings = FALSE))

# 4. Define bacterial species of interest
bacteria <- c(
  "Staphylococcus aureus",
  "Pseudomonas aeruginosa",
  "Escherichia coli",
  "Klebsiella pneumoniae",
  "Acinetobacter baumannii",
  "Enterococcus faecium",
  "Enterococcus faecalis"
)

# 5. Define small-RNA terms
srna_terms <- c(
  "\"small RNA\"",
  "\"small regulatory RNA\"",
  "\"regulatory RNA\"",
  "sRNA",
  "\"non-coding RNA\"",
  "\"noncoding RNA\"",
  "\"small RNA sequencing\"",
  "sRNA-seq"
)

# 6. Define antimicrobial-resistance terms
resistance_terms <- c(
  "\"antibiotic resistance\"",
  "\"antimicrobial resistance\"",
  "\"multidrug resistance\"",
  "antibiotic",
  "biofilm",
  "virulence",
  "persister",
  "tolerance"
)

# 7. Function to build an NCBI search query
build_query <- function(bacterium, srna_term, resistance_term) {
  
  paste0(
    "\"", bacterium, "\"[Organism] AND ",
    "(", srna_term, ") AND ",
    "(", resistance_term, ")"
  )
}

# 8. Construct all search combinations
search_grid <- tidyr::crossing(
  bacterium = bacteria,
  srna_term = srna_terms,
  resistance_term = resistance_terms
) |>
  mutate(
    query = pmap_chr(
      list(bacterium, srna_term, resistance_term),
      build_query
    )
  )

# 9. Function to search one NCBI database
search_ncbi <- function(query, database = "sra", max_results = 500) {
  
  tryCatch({
    
    result <- entrez_search(
      db = database,
      term = query,
      retmax = max_results,
      use_history = TRUE
    )
    
    tibble(
      database = database,
      query = query,
      number_found = result$count,
      record_ids = paste(result$ids, collapse = ";")
    )
    
  }, error = function(e) {
    
    tibble(
      database = database,
      query = query,
      number_found = NA_integer_,
      record_ids = NA_character_
    )
  })
}

# 10. Search SRA
message("Searching SRA...")

sra_search_results <- map_dfr(
  search_grid$query,
  search_ncbi,
  database = "sra"
)

# 11. Search GEO DataSets
message("Searching GEO DataSets...")

geo_search_results <- map_dfr(
  search_grid$query,
  search_ncbi,
  database = "gds"
)

# 12. Combine search results
all_search_results <- bind_rows(
  sra_search_results,
  geo_search_results
) |>
  left_join(search_grid, by = "query") |>
  select(
    database,
    bacterium,
    srna_term,
    resistance_term,
    number_found,
    query,
    record_ids
  ) |>
  arrange(
    desc(number_found),
    database,
    bacterium
  )

# 13. Save full results
write_csv(
  all_search_results,
  "01_search_results/BacRegRNA_initial_search_results.csv"
)

# 14. Create a summary by bacterium and database
search_summary <- all_search_results |>
  group_by(database, bacterium) |>
  summarise(
    total_query_hits = sum(number_found, na.rm = TRUE),
    queries_with_results = sum(number_found > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(database, desc(total_query_hits))

write_csv(
  search_summary,
  "01_search_results/BacRegRNA_search_summary.csv"
)

# 15. Print summary
print(search_summary)

message(
  "\nSearch completed.\n",
  "Files saved in the folder: 01_search_results"
)
# ============================================================
# Save all outputs directly inside D:/Bac-sRNA
# ============================================================

output_folder <- "D:/Bac-sRNA"

# Create the folder if it does not exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Save the complete search results
readr::write_csv(
  all_search_results,
  file.path(output_folder, "BacRegRNA_initial_search_results.csv")
)

# Save the summary table
readr::write_csv(
  search_summary,
  file.path(output_folder, "BacRegRNA_search_summary.csv")
)

# Save an Excel-compatible tab-separated copy as backup
readr::write_tsv(
  all_search_results,
  file.path(output_folder, "BacRegRNA_initial_search_results.tsv")
)

# Confirm that the files were created
saved_files <- c(
  file.path(output_folder, "BacRegRNA_initial_search_results.csv"),
  file.path(output_folder, "BacRegRNA_search_summary.csv"),
  file.path(output_folder, "BacRegRNA_initial_search_results.tsv")
)

cat("\nFiles saved successfully in:\n", output_folder, "\n\n")

print(
  data.frame(
    file = basename(saved_files),
    exists = file.exists(saved_files),
    size_KB = round(file.info(saved_files)$size / 1024, 2)
  )
)