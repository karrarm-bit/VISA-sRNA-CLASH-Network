# ============================================================
# BacRegRNA Project
# Script 04E: Download GEO MINiML files using curl
# ============================================================

# ------------------------------------------------------------
# 1. Project folders
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

output_folder <- file.path(
  project_folder,
  "04_metadata",
  "manual_download"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

setwd(project_folder)

# ------------------------------------------------------------
# 2. Install curl only if needed
# ------------------------------------------------------------

if (!requireNamespace("curl", quietly = TRUE)) {
  install.packages("curl")
}

library(curl)

# ------------------------------------------------------------
# 3. Datasets to download
# ------------------------------------------------------------

datasets <- c(
  "GSE254530",
  "GSE254531",
  "GSE254532",
  "GSE158830"
)

# ------------------------------------------------------------
# 4. Function to create GEO folder name
# Example:
# GSE254530 -> GSE254nnn
# ------------------------------------------------------------

geo_parent_folder <- function(gse) {
  
  gse_number <- as.integer(
    sub("^GSE", "", gse)
  )
  
  prefix <- floor(gse_number / 1000)
  
  paste0(
    "GSE",
    prefix,
    "nnn"
  )
}

# ------------------------------------------------------------
# 5. Download one GEO MINiML file
# ------------------------------------------------------------

download_one_geo <- function(gse) {
  
  parent_folder <- geo_parent_folder(gse)
  
  url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    parent_folder,
    "/",
    gse,
    "/miniml/",
    gse,
    "_family.xml.tgz"
  )
  
  destination <- file.path(
    output_folder,
    paste0(
      gse,
      "_family.xml.tgz"
    )
  )
  
  cat("\n============================================\n")
  cat("Dataset:", gse, "\n")
  cat("URL:", url, "\n")
  cat("Destination:", destination, "\n")
  cat("============================================\n")
  
  # Skip if file already exists and is not empty
  if (
    file.exists(destination) &&
    file.info(destination)$size > 1000
  ) {
    
    cat(
      "File already exists. Download skipped.\n"
    )
    
    return(
      data.frame(
        dataset = gse,
        status = "Already exists",
        file = destination,
        size_bytes = file.info(destination)$size,
        stringsAsFactors = FALSE
      )
    )
  }
  
  handle <- new_handle()
  
  handle_setopt(
    handle,
    useragent = paste0(
      "Mozilla/5.0 (Windows NT 10.0; Win64; x64) ",
      "AppleWebKit/537.36 ",
      "Chrome/150.0 Safari/537.36"
    ),
    followlocation = TRUE,
    timeout = 600,
    connecttimeout = 60
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
        "Download failed for ",
        gse,
        ": ",
        conditionMessage(e)
      )
      
      "Failed"
    }
  )
  
  file_size <- if (
    file.exists(destination)
  ) {
    file.info(destination)$size
  } else {
    0
  }
  
  # Remove failed or empty files
  if (
    status == "Failed" ||
    is.na(file_size) ||
    file_size < 1000
  ) {
    
    if (file.exists(destination)) {
      file.remove(destination)
    }
    
    file_size <- 0
  }
  
  data.frame(
    dataset = gse,
    status = status,
    file = destination,
    size_bytes = file_size,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 6. Download all datasets
# ------------------------------------------------------------

download_results <- do.call(
  rbind,
  lapply(
    datasets,
    download_one_geo
  )
)

# ------------------------------------------------------------
# 7. Save download report
# ------------------------------------------------------------

report_file <- file.path(
  output_folder,
  "GEO_download_report.csv"
)

write.csv(
  download_results,
  report_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Print final results
# ------------------------------------------------------------

cat("\n============================================\n")
cat("DOWNLOAD PROCESS FINISHED\n")
cat("============================================\n\n")

print(download_results)

cat(
  "\nFiles should be inside:\n",
  output_folder,
  "\n"
)

cat(
  "\nDownload report saved as:\n",
  report_file,
  "\n"
)

# ------------------------------------------------------------
# 9. Check successful files
# ------------------------------------------------------------

successful_files <- download_results[
  download_results$size_bytes > 1000,
]

cat(
  "\nSuccessfully downloaded files:",
  nrow(successful_files),
  "out of",
  nrow(download_results),
  "\n"
)

if (nrow(successful_files) == 0) {
  
  cat(
    "\nNo GEO files were downloaded.\n",
    "This means NCBI is still blocking automated downloads ",
    "from your connection or RStudio session.\n",
    "We will then use manual browser download.\n"
  )
  
} else {
  
  cat(
    "\nSuccessful files:\n"
  )
  
  print(
    successful_files[
      ,
      c(
        "dataset",
        "status",
        "size_bytes",
        "file"
      )
    ]
  )
}
files <- list.files(
  "D:/Bac-sRNA/04_metadata/manual_download",
  pattern = "\\.tgz$",
  full.names = TRUE
)

check_files <- data.frame(
  file = basename(files),
  size_bytes = file.info(files)$size,
  first_bytes = sapply(
    files,
    function(f) {
      con <- file(f, "rb")
      on.exit(close(con), add = TRUE)
      rawToChar(readBin(con, "raw", n = 20))
    }
  )
)

print(check_files)
files <- list.files(
  "D:/Bac-sRNA/04_metadata/manual_download",
  pattern = "\\.tgz$",
  full.names = TRUE
)

extract_folder <- "D:/Bac-sRNA/04_metadata/extracted"

dir.create(
  extract_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

for(f in files){
  
  cat("Extracting:", basename(f), "\n")
  
  untar(
    tarfile = f,
    exdir = file.path(
      extract_folder,
      tools::file_path_sans_ext(
        tools::file_path_sans_ext(
          basename(f)
        )
      )
    )
  )
  
}

cat("\nExtraction finished.\n")
list.files(
  "D:/Bac-sRNA/04_metadata/extracted",
  recursive = TRUE,
  full.names = TRUE
)