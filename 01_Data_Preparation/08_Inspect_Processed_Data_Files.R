# ============================================================
# BacRegRNA Project
# Script 08: Inspect downloaded processed-data files
# ============================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(tidyr)

project_folder <- "D:/Bac-sRNA"

input_folder <- file.path(
  project_folder,
  "07_processed_data",
  "downloaded_files"
)

output_folder <- file.path(
  project_folder,
  "08_file_inspection"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Find downloaded files
# ------------------------------------------------------------

files <- list.files(
  input_folder,
  recursive = TRUE,
  full.names = TRUE
)

files <- files[
  !dir.exists(files)
]

if (length(files) == 0) {
  stop("No downloaded files were found.")
}

file_table <- tibble(
  file_path = files,
  dataset = str_extract(
    files,
    "GSE[0-9]+"
  ),
  file_name = basename(files),
  extension = str_to_lower(
    tools::file_ext(files)
  ),
  size_bytes = file.info(files)$size
)

# ------------------------------------------------------------
# 2. Read first lines safely
# ------------------------------------------------------------

read_preview <- function(path, n_lines = 8) {
  
  tryCatch(
    {
      if (str_detect(path, "\\.gz$")) {
        connection <- gzfile(path, open = "rt")
      } else {
        connection <- file(path, open = "rt")
      }
      
      on.exit(close(connection), add = TRUE)
      
      lines <- readLines(
        connection,
        n = n_lines,
        warn = FALSE
      )
      
      paste(lines, collapse = "\n")
    },
    error = function(e) {
      paste0(
        "PREVIEW_FAILED: ",
        conditionMessage(e)
      )
    }
  )
}

preview_table <- file_table |>
  mutate(
    preview = map_chr(
      file_path,
      read_preview
    )
  )

write_csv(
  preview_table,
  file.path(
    output_folder,
    "all_file_previews.csv"
  )
)

# ------------------------------------------------------------
# 3. Inspect selected key files in detail
# ------------------------------------------------------------

key_patterns <- c(
  "GSM8045612",
  "GSM8045615",
  "GSE254530_RNAseq_VAN30_DEseq2",
  "GSE254532_Final_RATT_VISACLASH",
  "GSE158830_MasterTable_GEO",
  "delta_275_vs_WT",
  "CRISPRi_275_vs_pSD",
  "MUTvsWT"
)

key_files <- file_table |>
  filter(
    str_detect(
      file_name,
      paste(key_patterns, collapse = "|")
    )
  )

# ------------------------------------------------------------
# 4. Guess delimiters and columns
# ------------------------------------------------------------

inspect_tabular_file <- function(path) {
  
  result <- tryCatch(
    {
      
      if (str_detect(path, "\\.xlsx$")) {
        return(
          tibble(
            readable = "Excel file",
            rows_detected = NA_integer_,
            columns_detected = NA_integer_,
            column_names = NA_character_
          )
        )
      }
      
      if (str_detect(path, "\\.xls\\.gz$")) {
        return(
          tibble(
            readable = "Compressed Excel-like file",
            rows_detected = NA_integer_,
            columns_detected = NA_integer_,
            column_names = NA_character_
          )
        )
      }
      
      data <- suppressMessages(
        read_delim(
          path,
          delim = NULL,
          n_max = 50,
          show_col_types = FALSE,
          progress = FALSE
        )
      )
      
      tibble(
        readable = "Yes",
        rows_detected = nrow(data),
        columns_detected = ncol(data),
        column_names = paste(
          names(data),
          collapse = " | "
        )
      )
    },
    error = function(e) {
      
      tibble(
        readable = paste0(
          "Failed: ",
          conditionMessage(e)
        ),
        rows_detected = NA_integer_,
        columns_detected = NA_integer_,
        column_names = NA_character_
      )
    }
  )
  
  result
}

key_file_inspection <- key_files |>
  mutate(
    inspection = map(
      file_path,
      inspect_tabular_file
    )
  ) |>
  unnest(inspection)

write_csv(
  key_file_inspection,
  file.path(
    output_folder,
    "key_file_structure.csv"
  )
)

# ------------------------------------------------------------
# 5. Save human-readable preview text
# ------------------------------------------------------------

preview_output <- file.path(
  output_folder,
  "key_file_previews.txt"
)

sink(preview_output)

cat("BACREGRNA KEY FILE PREVIEWS\n")
cat("============================\n\n")

for (i in seq_len(nrow(key_files))) {
  
  cat("\n--------------------------------------------\n")
  cat("Dataset:", key_files$dataset[i], "\n")
  cat("File:", key_files$file_name[i], "\n")
  cat("--------------------------------------------\n")
  
  cat(
    read_preview(
      key_files$file_path[i],
      n_lines = 12
    )
  )
  
  cat("\n")
}

sink()

# ------------------------------------------------------------
# 6. Final summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 08 COMPLETED\n")
cat("============================================\n")

cat("Total files inspected:", nrow(file_table), "\n")
cat("Key files inspected:", nrow(key_files), "\n")

cat(
  "\nOutputs saved inside:\n",
  output_folder,
  "\n\n"
)

print(
  key_file_inspection |>
    select(
      dataset,
      file_name,
      readable,
      rows_detected,
      columns_detected,
      column_names
    ),
  n = Inf
)