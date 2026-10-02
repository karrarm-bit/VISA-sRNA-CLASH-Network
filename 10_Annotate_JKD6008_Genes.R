# ============================================================
# BacRegRNA Project
# Script 10: Annotate Staphylococcus aureus JKD6008 genes
#
# Reference strain:
#   Staphylococcus aureus subsp. aureus JKD6008
#
# Reference accessions:
#   NC_017341.1  = RefSeq chromosome
#   CP002120.1   = Original GenBank chromosome
#
# Main purposes:
#   1. Download strain-specific GenBank annotation.
#   2. Parse locus tags, old locus tags, gene symbols and products.
#   3. Map SAA6008 identifiers from GSE254530.
#   4. Annotate all DESeq2 results.
#   5. Produce annotation tables for publication figures.
# ============================================================


# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

deseq_results_file <- file.path(
  project_folder,
  "09_RNAseq_analysis",
  "tables",
  "GSE254530_DESeq2_results_all_genes.csv"
)

significant_results_file <- file.path(
  project_folder,
  "09_RNAseq_analysis",
  "tables",
  "GSE254530_DESeq2_significant_genes.csv"
)

output_folder <- file.path(
  project_folder,
  "10_annotation"
)

reference_folder <- file.path(
  output_folder,
  "reference_files"
)

table_folder <- file.path(
  output_folder,
  "tables"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  reference_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

setwd(project_folder)

options(timeout = 1800)


# ------------------------------------------------------------
# 2. Install and load packages
# ------------------------------------------------------------

required_packages <- c(
  "curl",
  "readr",
  "dplyr",
  "stringr",
  "purrr",
  "tibble",
  "tidyr"
)

missing_packages <- required_packages[
  !required_packages %in%
    rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(curl)
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(tidyr)


# ------------------------------------------------------------
# 3. Check DESeq2 result files
# ------------------------------------------------------------

if (!file.exists(deseq_results_file)) {
  stop(
    "DESeq2 results file was not found:\n",
    deseq_results_file
  )
}

if (!file.exists(significant_results_file)) {
  stop(
    "Significant-gene file was not found:\n",
    significant_results_file
  )
}


# ------------------------------------------------------------
# 4. Define JKD6008 reference accessions
# ------------------------------------------------------------

reference_table <- tibble::tribble(
  ~accession,      ~reference_type,       ~priority,
  "NC_017341.1",   "Current RefSeq",       1,
  "CP002120.1",    "Original GenBank",     2
)


# ------------------------------------------------------------
# 5. Download one GenBank record from NCBI
# ------------------------------------------------------------

download_genbank_record <- function(
    accession,
    reference_type
) {
  
  destination <- file.path(
    reference_folder,
    paste0(
      accession,
      "_JKD6008.gb"
    )
  )
  
  cat("\n============================================\n")
  cat("Reference:", accession, "\n")
  cat("Type:", reference_type, "\n")
  cat("============================================\n")
  
  if (
    file.exists(destination) &&
    !is.na(file.info(destination)$size) &&
    file.info(destination)$size > 100000
  ) {
    
    cat(
      "Existing valid GenBank file detected; ",
      "download skipped.\n"
    )
    
    return(
      tibble(
        accession = accession,
        reference_type = reference_type,
        local_file = destination,
        status = "Already exists",
        size_bytes = file.info(destination)$size
      )
    )
  }
  
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
    "efetch.fcgi?",
    "db=nuccore",
    "&id=", accession,
    "&rettype=gbwithparts",
    "&retmode=text"
  )
  
  handle <- curl::new_handle()
  
  curl::handle_setopt(
    handle,
    useragent = paste0(
      "BacRegRNA academic research ",
      "contact: local-R-analysis"
    ),
    followlocation = TRUE,
    timeout = 1200,
    connecttimeout = 120
  )
  
  curl::handle_setheaders(
    handle,
    Accept = "text/plain",
    `Accept-Language` = "en-US,en;q=0.9"
  )
  
  status <- tryCatch(
    {
      curl::curl_download(
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
        accession,
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
  
  # A complete bacterial GenBank file should be much larger
  # than a small HTML error page.
  if (
    is.na(file_size) ||
    file_size < 100000
  ) {
    
    status <- "Failed or invalid file"
    
    if (file.exists(destination)) {
      file.remove(destination)
    }
    
    file_size <- 0
  }
  
  tibble(
    accession = accession,
    reference_type = reference_type,
    local_file = destination,
    status = status,
    size_bytes = file_size
  )
}


# ------------------------------------------------------------
# 6. Download both reference records
# ------------------------------------------------------------

download_report <- purrr::pmap_dfr(
  reference_table |>
    select(
      accession,
      reference_type
    ),
  download_genbank_record
)

readr::write_csv(
  download_report,
  file.path(
    table_folder,
    "JKD6008_reference_download_report.csv"
  )
)

successful_references <- download_report |>
  filter(
    status %in% c(
      "Downloaded",
      "Already exists"
    ),
    size_bytes > 100000
  )

if (nrow(successful_references) == 0) {
  
  stop(
    "Neither JKD6008 GenBank reference could be downloaded.\n",
    "Check the NCBI connection or download NC_017341.1 ",
    "manually as a GenBank file."
  )
}


# ------------------------------------------------------------
# 7. Helper: extract one GenBank qualifier
#
# Examples:
#   /locus_tag="SAA6008_00001"
#   /gene="dnaA"
#   /product="chromosomal replication initiator protein"
# ------------------------------------------------------------

extract_qualifier <- function(
    feature_text,
    qualifier
) {
  
  quoted_pattern <- paste0(
    "/",
    qualifier,
    "=\"([^\"]*)\""
  )
  
  quoted_match <- stringr::str_match(
    feature_text,
    stringr::regex(
      quoted_pattern,
      dotall = TRUE
    )
  )[, 2]
  
  if (
    !is.na(quoted_match) &&
    quoted_match != ""
  ) {
    
    return(
      stringr::str_squish(
        quoted_match
      )
    )
  }
  
  unquoted_pattern <- paste0(
    "/",
    qualifier,
    "=([^ /]+)"
  )
  
  unquoted_match <- stringr::str_match(
    feature_text,
    unquoted_pattern
  )[, 2]
  
  if (
    is.na(unquoted_match) ||
    unquoted_match == ""
  ) {
    return(NA_character_)
  }
  
  stringr::str_squish(
    unquoted_match
  )
}


# ------------------------------------------------------------
# 8. Helper: extract feature type and location
# ------------------------------------------------------------

extract_feature_header <- function(
    first_line
) {
  
  feature_type <- stringr::str_match(
    first_line,
    "^\\s{5}(\\S+)"
  )[, 2]
  
  feature_location <- stringr::str_match(
    first_line,
    "^\\s{5}\\S+\\s+(.+)$"
  )[, 2]
  
  tibble(
    feature_type = feature_type,
    feature_location = stringr::str_squish(
      feature_location
    )
  )
}


# ------------------------------------------------------------
# 9. Parse one GenBank file
# ------------------------------------------------------------

parse_genbank_features <- function(
    genbank_file,
    accession,
    reference_type
) {
  
  cat("\nParsing GenBank reference:", accession, "\n")
  
  lines <- readLines(
    genbank_file,
    warn = FALSE,
    encoding = "UTF-8"
  )
  
  features_start <- which(
    stringr::str_detect(
      lines,
      "^FEATURES\\s+Location/Qualifiers"
    )
  )[1]
  
  origin_start <- which(
    stringr::str_detect(
      lines,
      "^ORIGIN"
    )
  )[1]
  
  if (
    is.na(features_start) ||
    is.na(origin_start) ||
    origin_start <= features_start
  ) {
    
    stop(
      "The GenBank FEATURES section could not be detected in:\n",
      genbank_file
    )
  }
  
  feature_lines <- lines[
    (features_start + 1):
      (origin_start - 1)
  ]
  
  feature_start_indices <- which(
    stringr::str_detect(
      feature_lines,
      "^\\s{5}\\S"
    )
  )
  
  if (length(feature_start_indices) == 0) {
    
    stop(
      "No GenBank features were detected in:\n",
      genbank_file
    )
  }
  
  feature_end_indices <- c(
    feature_start_indices[-1] - 1,
    length(feature_lines)
  )
  
  feature_table <- purrr::map2_dfr(
    feature_start_indices,
    feature_end_indices,
    function(start_index, end_index) {
      
      block_lines <- feature_lines[
        start_index:end_index
      ]
      
      header <- extract_feature_header(
        block_lines[1]
      )
      
      block_text <- paste(
        block_lines,
        collapse = " "
      ) |>
        stringr::str_squish()
      
      tibble(
        reference_accession = accession,
        reference_type = reference_type,
        
        feature_type =
          header$feature_type,
        
        feature_location =
          header$feature_location,
        
        locus_tag = extract_qualifier(
          block_text,
          "locus_tag"
        ),
        
        old_locus_tag = extract_qualifier(
          block_text,
          "old_locus_tag"
        ),
        
        gene_symbol = extract_qualifier(
          block_text,
          "gene"
        ),
        
        product = extract_qualifier(
          block_text,
          "product"
        ),
        
        protein_id = extract_qualifier(
          block_text,
          "protein_id"
        ),
        
        gene_synonym = extract_qualifier(
          block_text,
          "gene_synonym"
        ),
        
        ncRNA_class = extract_qualifier(
          block_text,
          "ncRNA_class"
        ),
        
        note = extract_qualifier(
          block_text,
          "note"
        ),
        
        db_xref = extract_qualifier(
          block_text,
          "db_xref"
        )
      )
    }
  )
  
  feature_table |>
    filter(
      feature_type %in% c(
        "CDS",
        "gene",
        "rRNA",
        "tRNA",
        "tmRNA",
        "ncRNA",
        "misc_RNA"
      )
    )
}


# ------------------------------------------------------------
# 10. Parse all successfully downloaded references
# ------------------------------------------------------------

all_features <- successful_references |>
  left_join(
    reference_table,
    by = c(
      "accession",
      "reference_type"
    )
  ) |>
  arrange(priority) |>
  transmute(
    accession,
    reference_type,
    priority,
    local_file
  ) |>
  purrr::pmap_dfr(
    function(
    accession,
    reference_type,
    priority,
    local_file
    ) {
      
      parsed <- parse_genbank_features(
        genbank_file = local_file,
        accession = accession,
        reference_type = reference_type
      )
      
      parsed |>
        mutate(
          reference_priority = priority
        )
    }
  )

readr::write_csv(
  all_features,
  file.path(
    table_folder,
    "JKD6008_all_parsed_GenBank_features.csv"
  )
)


# ------------------------------------------------------------
# 11. Assign the SAA6008 gene identifier
#
# New RefSeq records may store SAA6008 as old_locus_tag.
# Original GenBank records may store it as locus_tag.
# ------------------------------------------------------------

features_with_gene_id <- all_features |>
  mutate(
    
    gene_id = case_when(
      
      stringr::str_detect(
        locus_tag,
        "^SAA6008_[0-9]+$"
      ) ~ locus_tag,
      
      stringr::str_detect(
        old_locus_tag,
        "^SAA6008_[0-9]+$"
      ) ~ old_locus_tag,
      
      TRUE ~ NA_character_
    ),
    
    current_locus_tag = case_when(
      
      !is.na(locus_tag) &
        locus_tag != gene_id ~ locus_tag,
      
      TRUE ~ NA_character_
    )
  ) |>
  filter(
    !is.na(gene_id)
  )


# ------------------------------------------------------------
# 12. Prioritize biologically informative feature types
#
# CDS is preferred for protein-coding genes.
# RNA-specific feature types are retained for RNA genes.
# ------------------------------------------------------------

features_with_gene_id <- features_with_gene_id |>
  mutate(
    
    feature_priority = case_when(
      feature_type == "CDS" ~ 1,
      feature_type == "ncRNA" ~ 2,
      feature_type == "rRNA" ~ 2,
      feature_type == "tRNA" ~ 2,
      feature_type == "tmRNA" ~ 2,
      feature_type == "misc_RNA" ~ 3,
      feature_type == "gene" ~ 4,
      TRUE ~ 5
    ),
    
    information_score =
      as.integer(
        !is.na(gene_symbol) &
          gene_symbol != ""
      ) * 4 +
      
      as.integer(
        !is.na(product) &
          product != ""
      ) * 3 +
      
      as.integer(
        !is.na(protein_id) &
          protein_id != ""
      ) * 2 +
      
      as.integer(
        !is.na(current_locus_tag) &
          current_locus_tag != ""
      )
  )


# ------------------------------------------------------------
# 13. Select the best annotation record per SAA6008 gene
# ------------------------------------------------------------

best_annotation <- features_with_gene_id |>
  arrange(
    gene_id,
    reference_priority,
    feature_priority,
    desc(information_score)
  ) |>
  group_by(gene_id) |>
  slice_head(n = 1) |>
  ungroup()


# ------------------------------------------------------------
# 14. Recover missing information from duplicate feature rows
# ------------------------------------------------------------

annotation_collapsed <- features_with_gene_id |>
  arrange(
    gene_id,
    reference_priority,
    feature_priority,
    desc(information_score)
  ) |>
  group_by(gene_id) |>
  summarise(
    
    current_locus_tag = {
      values <- na.omit(current_locus_tag)
      if (length(values) == 0) {
        NA_character_
      } else {
        values[1]
      }
    },
    
    gene_symbol = {
      values <- na.omit(gene_symbol)
      values <- values[
        values != ""
      ]
      
      if (length(values) == 0) {
        NA_character_
      } else {
        values[1]
      }
    },
    
    gene_synonym = {
      values <- na.omit(gene_synonym)
      values <- values[
        values != ""
      ]
      
      if (length(values) == 0) {
        NA_character_
      } else {
        values[1]
      }
    },
    
    product = {
      values <- na.omit(product)
      values <- values[
        values != ""
      ]
      
      if (length(values) == 0) {
        NA_character_
      } else {
        values[1]
      }
    },
    
    protein_id = {
      values <- na.omit(protein_id)
      values <- values[
        values != ""
      ]
      
      if (length(values) == 0) {
        NA_character_
      } else {
        values[1]
      }
    },
    
    feature_type = {
      values <- na.omit(feature_type)
      if (length(values) == 0) {
        NA_character_
      } else {
        values[1]
      }
    },
    
    feature_location = {
      values <- na.omit(feature_location)
      if (length(values) == 0) {
        NA_character_
      } else {
        values[1]
      }
    },
    
    ncRNA_class = {
      values <- na.omit(ncRNA_class)
      values <- values[
        values != ""
      ]
      
      if (length(values) == 0) {
        NA_character_
      } else {
        values[1]
      }
    },
    
    annotation_source = paste(
      unique(reference_accession),
      collapse = "; "
    ),
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 15. Clean gene symbols and product descriptions
# ------------------------------------------------------------

annotation_table <- annotation_collapsed |>
  mutate(
    
    gene_symbol = na_if(
      str_squish(gene_symbol),
      ""
    ),
    
    product = na_if(
      str_squish(product),
      ""
    ),
    
    # Some records may use the locus tag as a pseudo-symbol.
    valid_gene_symbol = case_when(
      
      is.na(gene_symbol) ~ NA_character_,
      
      gene_symbol == gene_id ~ NA_character_,
      
      !is.na(current_locus_tag) &
        gene_symbol == current_locus_tag ~
        NA_character_,
      
      str_detect(
        gene_symbol,
        "^SAA6008_[0-9]+$"
      ) ~ NA_character_,
      
      TRUE ~ gene_symbol
    ),
    
    product_short = product |>
      str_remove(
        regex(
          "^putative\\s+",
          ignore_case = TRUE
        )
      ) |>
      str_remove(
        regex(
          "^probable\\s+",
          ignore_case = TRUE
        )
      ) |>
      str_replace(
        regex(
          " family protein$",
          ignore_case = TRUE
        ),
        " protein"
      ) |>
      str_squish(),
    
    # Figure-label hierarchy:
    # gene symbol -> short product -> original locus tag
    display_label = case_when(
      
      !is.na(valid_gene_symbol) &
        valid_gene_symbol != "" ~
        valid_gene_symbol,
      
      !is.na(product_short) &
        product_short != "" &
        str_length(product_short) <= 38 ~
        product_short,
      
      !is.na(product_short) &
        product_short != "" ~
        paste0(
          str_sub(
            product_short,
            1,
            35
          ),
          "..."
        ),
      
      TRUE ~ gene_id
    ),
    
    annotation_status = case_when(
      
      !is.na(valid_gene_symbol) ~
        "Gene symbol available",
      
      !is.na(product) ~
        "Product annotation only",
      
      TRUE ~
        "Locus tag only"
    )
  ) |>
  select(
    gene_id,
    current_locus_tag,
    valid_gene_symbol,
    gene_synonym,
    product,
    product_short,
    protein_id,
    feature_type,
    ncRNA_class,
    feature_location,
    display_label,
    annotation_status,
    annotation_source
  ) |>
  rename(
    gene_symbol = valid_gene_symbol
  ) |>
  arrange(gene_id)


# ------------------------------------------------------------
# 16. Save whole-genome annotation database
# ------------------------------------------------------------

readr::write_csv(
  annotation_table,
  file.path(
    table_folder,
    "JKD6008_genome_annotation.csv"
  )
)


# ------------------------------------------------------------
# 17. Read DESeq2 result tables
# ------------------------------------------------------------

all_deseq_results <- readr::read_csv(
  deseq_results_file,
  show_col_types = FALSE
)

significant_deseq_results <- readr::read_csv(
  significant_results_file,
  show_col_types = FALSE
)

if (!"gene_id" %in% names(all_deseq_results)) {
  stop(
    "The all-gene DESeq2 table does not contain gene_id."
  )
}

if (!"gene_id" %in% names(significant_deseq_results)) {
  stop(
    "The significant-gene table does not contain gene_id."
  )
}


# ------------------------------------------------------------
# 18. Annotate all DESeq2 genes
# ------------------------------------------------------------

annotated_all_results <- all_deseq_results |>
  left_join(
    annotation_table,
    by = "gene_id"
  ) |>
  mutate(
    
    display_label = coalesce(
      display_label,
      gene_id
    ),
    
    annotation_status = coalesce(
      annotation_status,
      "Unmapped"
    ),
    
    annotation_source = coalesce(
      annotation_source,
      "No matching GenBank annotation"
    )
  ) |>
  relocate(
    gene_id,
    current_locus_tag,
    gene_symbol,
    product,
    product_short,
    protein_id,
    display_label,
    annotation_status,
    annotation_source
  )

readr::write_csv(
  annotated_all_results,
  file.path(
    table_folder,
    "GSE254530_DESeq2_annotated_all_genes.csv"
  )
)


# ------------------------------------------------------------
# 19. Annotate significant DEGs
# ------------------------------------------------------------

annotated_significant_results <-
  significant_deseq_results |>
  left_join(
    annotation_table,
    by = "gene_id"
  ) |>
  mutate(
    
    display_label = coalesce(
      display_label,
      gene_id
    ),
    
    annotation_status = coalesce(
      annotation_status,
      "Unmapped"
    ),
    
    annotation_source = coalesce(
      annotation_source,
      "No matching GenBank annotation"
    )
  ) |>
  relocate(
    gene_id,
    current_locus_tag,
    gene_symbol,
    product,
    product_short,
    protein_id,
    display_label,
    annotation_status,
    annotation_source
  )

readr::write_csv(
  annotated_significant_results,
  file.path(
    table_folder,
    "GSE254530_DESeq2_annotated_significant_genes.csv"
  )
)


# ------------------------------------------------------------
# 20. Create a unique figure label
#
# Duplicate symbols/products are supplemented with the locus tag
# to prevent duplicated heatmap row names.
# ------------------------------------------------------------

annotated_significant_results <-
  annotated_significant_results |>
  group_by(display_label) |>
  mutate(
    
    figure_label = case_when(
      
      n() == 1 ~ display_label,
      
      TRUE ~ paste0(
        display_label,
        " [",
        gene_id,
        "]"
      )
    )
  ) |>
  ungroup()

readr::write_csv(
  annotated_significant_results,
  file.path(
    table_folder,
    "GSE254530_DESeq2_annotated_significant_genes.csv"
  )
)


# ------------------------------------------------------------
# 21. Prepare top genes for the final volcano plot
#
# Select highly significant genes from both directions.
# ------------------------------------------------------------

top_volcano_labels <- bind_rows(
  
  annotated_significant_results |>
    filter(
      significance == "Upregulated",
      !is.na(adjusted_p_value)
    ) |>
    arrange(
      adjusted_p_value,
      desc(abs(log2_fold_change))
    ) |>
    slice_head(n = 6),
  
  annotated_significant_results |>
    filter(
      significance == "Downregulated",
      !is.na(adjusted_p_value)
    ) |>
    arrange(
      adjusted_p_value,
      desc(abs(log2_fold_change))
    ) |>
    slice_head(n = 6)
) |>
  distinct(
    gene_id,
    .keep_all = TRUE
  )

readr::write_csv(
  top_volcano_labels,
  file.path(
    table_folder,
    "GSE254530_top_volcano_labels.csv"
  )
)


# ------------------------------------------------------------
# 22. Prepare top 30 genes for final heatmap
#
# Balanced selection:
#   15 upregulated
#   15 downregulated
# ------------------------------------------------------------

top_heatmap_genes <- bind_rows(
  
  annotated_significant_results |>
    filter(
      significance == "Upregulated",
      !is.na(adjusted_p_value)
    ) |>
    arrange(
      adjusted_p_value,
      desc(abs(log2_fold_change))
    ) |>
    slice_head(n = 15),
  
  annotated_significant_results |>
    filter(
      significance == "Downregulated",
      !is.na(adjusted_p_value)
    ) |>
    arrange(
      adjusted_p_value,
      desc(abs(log2_fold_change))
    ) |>
    slice_head(n = 15)
) |>
  distinct(
    gene_id,
    .keep_all = TRUE
  )

readr::write_csv(
  top_heatmap_genes,
  file.path(
    table_folder,
    "GSE254530_top30_annotated_heatmap_genes.csv"
  )
)


# ------------------------------------------------------------
# 23. Save unmapped genes for manual review
# ------------------------------------------------------------

unmapped_genes <- annotated_all_results |>
  filter(
    annotation_status == "Unmapped"
  ) |>
  select(
    gene_id,
    base_mean,
    log2_fold_change,
    adjusted_p_value,
    significance
  ) |>
  arrange(
    adjusted_p_value
  )

readr::write_csv(
  unmapped_genes,
  file.path(
    table_folder,
    "GSE254530_unmapped_gene_ids.csv"
  )
)


# ------------------------------------------------------------
# 24. Annotation coverage summary
# ------------------------------------------------------------

annotation_summary <- annotated_all_results |>
  group_by(annotation_status) |>
  summarise(
    number_of_genes = dplyr::n(),
    .groups = "drop"
  ) |>
  mutate(
    percentage = round(
      100 *
        number_of_genes /
        sum(number_of_genes),
      2
    )
  ) |>
  arrange(
    desc(number_of_genes)
  )

significant_annotation_summary <-
  annotated_significant_results |>
  group_by(annotation_status) |>
  summarise(
    number_of_genes = dplyr::n(),
    .groups = "drop"
  ) |>
  mutate(
    percentage = round(
      100 *
        number_of_genes /
        sum(number_of_genes),
      2
    )
  ) |>
  arrange(
    desc(number_of_genes)
  )

readr::write_csv(
  annotation_summary,
  file.path(
    table_folder,
    "GSE254530_annotation_coverage_all_genes.csv"
  )
)

readr::write_csv(
  significant_annotation_summary,
  file.path(
    table_folder,
    "GSE254530_annotation_coverage_significant_genes.csv"
  )
)


# ------------------------------------------------------------
# 25. Save an audit log for reproducibility
# ------------------------------------------------------------

annotation_audit <- tibble(
  item = c(
    "Organism",
    "Strain",
    "Genome assembly",
    "Current RefSeq chromosome",
    "Original GenBank chromosome",
    "DESeq2 gene identifier pattern",
    "Annotation priority"
  ),
  
  value = c(
    "Staphylococcus aureus subsp. aureus",
    "JKD6008",
    "GCA_000145595.1",
    "NC_017341.1",
    "CP002120.1",
    "SAA6008_[0-9]+",
    paste(
      "Gene symbol;",
      "short product;",
      "original locus tag"
    )
  )
)

readr::write_csv(
  annotation_audit,
  file.path(
    table_folder,
    "JKD6008_annotation_audit_log.csv"
  )
)


# ------------------------------------------------------------
# 26. Save session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    output_folder,
    "sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 27. Final console report
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 10 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat(
  "Genome annotation records:",
  nrow(annotation_table),
  "\n"
)

cat(
  "DESeq2 genes:",
  nrow(annotated_all_results),
  "\n"
)

cat(
  "Significant DEGs:",
  nrow(annotated_significant_results),
  "\n"
)

cat(
  "Genes with gene symbols:",
  sum(
    annotated_all_results$
      annotation_status ==
      "Gene symbol available"
  ),
  "\n"
)

cat(
  "Genes with product annotation only:",
  sum(
    annotated_all_results$
      annotation_status ==
      "Product annotation only"
  ),
  "\n"
)

cat(
  "Unmapped genes:",
  nrow(unmapped_genes),
  "\n\n"
)

cat("Annotation coverage for all genes:\n")

print(
  annotation_summary,
  n = Inf
)

cat(
  "\nTop volcano labels:\n"
)

print(
  top_volcano_labels |>
    select(
      gene_id,
      gene_symbol,
      product,
      figure_label,
      log2_fold_change,
      adjusted_p_value
    ),
  n = Inf
)

cat(
  "\nOutputs saved inside:\n",
  output_folder,
  "\n"
)
# ============================================================
# STEP 10 COMPLETION FIX
# Continue JKD6008 gene annotation without re-downloading files
# ============================================================

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)

# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

parsed_features_file <- file.path(
  project_folder,
  "10_annotation",
  "tables",
  "JKD6008_all_parsed_GenBank_features.csv"
)

deseq_results_file <- file.path(
  project_folder,
  "09_RNAseq_analysis",
  "tables",
  "GSE254530_DESeq2_results_all_genes.csv"
)

significant_results_file <- file.path(
  project_folder,
  "09_RNAseq_analysis",
  "tables",
  "GSE254530_DESeq2_significant_genes.csv"
)

output_folder <- file.path(
  project_folder,
  "10_annotation"
)

table_folder <- file.path(
  output_folder,
  "tables"
)

dir.create(
  table_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Check required files
# ------------------------------------------------------------

required_files <- c(
  parsed_features_file,
  deseq_results_file,
  significant_results_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Required files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ------------------------------------------------------------
# 3. Read parsed GenBank features
# ------------------------------------------------------------

all_features <- readr::read_csv(
  parsed_features_file,
  show_col_types = FALSE
)

cat(
  "Parsed feature records:",
  nrow(all_features),
  "\n"
)

# ------------------------------------------------------------
# 4. Identify the original SAA6008 identifier
# ------------------------------------------------------------

features_with_gene_id <- all_features |>
  dplyr::mutate(
    
    gene_id = dplyr::case_when(
      
      !is.na(locus_tag) &
        stringr::str_detect(
          locus_tag,
          "^SAA6008_[0-9]+$"
        ) ~ locus_tag,
      
      !is.na(old_locus_tag) &
        stringr::str_detect(
          old_locus_tag,
          "^SAA6008_[0-9]+$"
        ) ~ old_locus_tag,
      
      TRUE ~ NA_character_
    ),
    
    current_locus_tag = dplyr::case_when(
      
      !is.na(locus_tag) &
        locus_tag != gene_id ~ locus_tag,
      
      TRUE ~ NA_character_
    ),
    
    feature_priority = dplyr::case_when(
      feature_type == "CDS" ~ 1L,
      feature_type %in% c(
        "ncRNA",
        "rRNA",
        "tRNA",
        "tmRNA"
      ) ~ 2L,
      feature_type == "misc_RNA" ~ 3L,
      feature_type == "gene" ~ 4L,
      TRUE ~ 5L
    ),
    
    information_score =
      as.integer(
        !is.na(gene_symbol) &
          gene_symbol != ""
      ) * 4L +
      as.integer(
        !is.na(product) &
          product != ""
      ) * 3L +
      as.integer(
        !is.na(protein_id) &
          protein_id != ""
      ) * 2L +
      as.integer(
        !is.na(current_locus_tag) &
          current_locus_tag != ""
      )
  ) |>
  dplyr::filter(
    !is.na(gene_id)
  ) |>
  dplyr::arrange(
    gene_id,
    reference_priority,
    feature_priority,
    dplyr::desc(information_score)
  )

cat(
  "Feature rows linked to SAA6008 IDs:",
  nrow(features_with_gene_id),
  "\n"
)

# ------------------------------------------------------------
# 5. Helper to select the first useful value
# ------------------------------------------------------------

first_valid_character <- function(x) {
  
  x <- as.character(x)
  
  x <- x[
    !is.na(x) &
      stringr::str_squish(x) != ""
  ]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  stringr::str_squish(x[1])
}

# ------------------------------------------------------------
# 6. Collapse duplicate feature rows per gene
# ------------------------------------------------------------

annotation_collapsed <- features_with_gene_id |>
  dplyr::group_by(gene_id) |>
  dplyr::summarise(
    
    current_locus_tag =
      first_valid_character(current_locus_tag),
    
    original_gene_symbol =
      first_valid_character(gene_symbol),
    
    gene_synonym =
      first_valid_character(gene_synonym),
    
    product =
      first_valid_character(product),
    
    protein_id =
      first_valid_character(protein_id),
    
    feature_type =
      first_valid_character(feature_type),
    
    ncRNA_class =
      first_valid_character(ncRNA_class),
    
    feature_location =
      first_valid_character(feature_location),
    
    annotation_source = paste(
      unique(
        stats::na.omit(reference_accession)
      ),
      collapse = "; "
    ),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 7. Clean symbols and build publication labels
#
# This explicitly creates gene_symbol and avoids the previous
# valid_gene_symbol / rename error.
# ------------------------------------------------------------

annotation_table <- annotation_collapsed |>
  dplyr::mutate(
    
    original_gene_symbol = dplyr::na_if(
      stringr::str_squish(
        original_gene_symbol
      ),
      ""
    ),
    
    product = dplyr::na_if(
      stringr::str_squish(product),
      ""
    ),
    
    gene_symbol = dplyr::case_when(
      
      is.na(original_gene_symbol) ~
        NA_character_,
      
      original_gene_symbol == gene_id ~
        NA_character_,
      
      !is.na(current_locus_tag) &
        original_gene_symbol ==
        current_locus_tag ~
        NA_character_,
      
      stringr::str_detect(
        original_gene_symbol,
        "^SAA6008_[0-9]+$"
      ) ~ NA_character_,
      
      TRUE ~ original_gene_symbol
    ),
    
    product_short = product |>
      stringr::str_remove(
        stringr::regex(
          "^putative\\s+",
          ignore_case = TRUE
        )
      ) |>
      stringr::str_remove(
        stringr::regex(
          "^probable\\s+",
          ignore_case = TRUE
        )
      ) |>
      stringr::str_replace(
        stringr::regex(
          "\\s+family protein$",
          ignore_case = TRUE
        ),
        " protein"
      ) |>
      stringr::str_squish(),
    
    display_label = dplyr::case_when(
      
      !is.na(gene_symbol) &
        gene_symbol != "" ~
        gene_symbol,
      
      !is.na(product_short) &
        product_short != "" &
        stringr::str_length(
          product_short
        ) <= 35 ~
        product_short,
      
      !is.na(product_short) &
        product_short != "" ~
        paste0(
          stringr::str_sub(
            product_short,
            1,
            32
          ),
          "..."
        ),
      
      TRUE ~ gene_id
    ),
    
    annotation_status = dplyr::case_when(
      
      !is.na(gene_symbol) ~
        "Gene symbol available",
      
      !is.na(product) ~
        "Product annotation only",
      
      TRUE ~
        "Locus tag only"
    )
  ) |>
  dplyr::transmute(
    gene_id = gene_id,
    current_locus_tag = current_locus_tag,
    gene_symbol = gene_symbol,
    gene_synonym = gene_synonym,
    product = product,
    product_short = product_short,
    protein_id = protein_id,
    feature_type = feature_type,
    ncRNA_class = ncRNA_class,
    feature_location = feature_location,
    display_label = display_label,
    annotation_status = annotation_status,
    annotation_source = annotation_source
  ) |>
  dplyr::arrange(gene_id)

# ------------------------------------------------------------
# 8. Save whole-genome annotation
# ------------------------------------------------------------

readr::write_csv(
  annotation_table,
  file.path(
    table_folder,
    "JKD6008_genome_annotation.csv"
  )
)

# ------------------------------------------------------------
# 9. Read DESeq2 tables
# ------------------------------------------------------------

all_deseq_results <- readr::read_csv(
  deseq_results_file,
  show_col_types = FALSE
)

significant_deseq_results <-
  readr::read_csv(
    significant_results_file,
    show_col_types = FALSE
  )

# ------------------------------------------------------------
# 10. Annotate all analysed genes
# ------------------------------------------------------------

annotated_all_results <- all_deseq_results |>
  dplyr::left_join(
    annotation_table,
    by = "gene_id"
  ) |>
  dplyr::mutate(
    
    display_label = dplyr::coalesce(
      display_label,
      gene_id
    ),
    
    annotation_status =
      dplyr::coalesce(
        annotation_status,
        "Unmapped"
      ),
    
    annotation_source =
      dplyr::coalesce(
        annotation_source,
        "No matching GenBank annotation"
      )
  )

readr::write_csv(
  annotated_all_results,
  file.path(
    table_folder,
    "GSE254530_DESeq2_annotated_all_genes.csv"
  )
)

# ------------------------------------------------------------
# 11. Annotate significant DEGs
# ------------------------------------------------------------

annotated_significant_results <-
  significant_deseq_results |>
  dplyr::left_join(
    annotation_table,
    by = "gene_id"
  ) |>
  dplyr::mutate(
    
    display_label = dplyr::coalesce(
      display_label,
      gene_id
    ),
    
    annotation_status =
      dplyr::coalesce(
        annotation_status,
        "Unmapped"
      ),
    
    annotation_source =
      dplyr::coalesce(
        annotation_source,
        "No matching GenBank annotation"
      )
  ) |>
  dplyr::group_by(display_label) |>
  dplyr::mutate(
    
    figure_label = dplyr::case_when(
      
      dplyr::n() == 1L ~
        display_label,
      
      TRUE ~ paste0(
        display_label,
        " [",
        gene_id,
        "]"
      )
    )
  ) |>
  dplyr::ungroup()

readr::write_csv(
  annotated_significant_results,
  file.path(
    table_folder,
    "GSE254530_DESeq2_annotated_significant_genes.csv"
  )
)

# ------------------------------------------------------------
# 12. Select volcano labels
# Six upregulated and six downregulated
# ------------------------------------------------------------

top_volcano_labels <- dplyr::bind_rows(
  
  annotated_significant_results |>
    dplyr::filter(
      significance == "Upregulated",
      !is.na(adjusted_p_value)
    ) |>
    dplyr::arrange(
      adjusted_p_value,
      dplyr::desc(
        abs(log2_fold_change)
      )
    ) |>
    dplyr::slice_head(n = 6),
  
  annotated_significant_results |>
    dplyr::filter(
      significance == "Downregulated",
      !is.na(adjusted_p_value)
    ) |>
    dplyr::arrange(
      adjusted_p_value,
      dplyr::desc(
        abs(log2_fold_change)
      )
    ) |>
    dplyr::slice_head(n = 6)
) |>
  dplyr::distinct(
    gene_id,
    .keep_all = TRUE
  )

readr::write_csv(
  top_volcano_labels,
  file.path(
    table_folder,
    "GSE254530_top_volcano_labels.csv"
  )
)

# ------------------------------------------------------------
# 13. Select top 30 heatmap genes
# 15 upregulated + 15 downregulated
# ------------------------------------------------------------

top_heatmap_genes <- dplyr::bind_rows(
  
  annotated_significant_results |>
    dplyr::filter(
      significance == "Upregulated",
      !is.na(adjusted_p_value)
    ) |>
    dplyr::arrange(
      adjusted_p_value,
      dplyr::desc(
        abs(log2_fold_change)
      )
    ) |>
    dplyr::slice_head(n = 15),
  
  annotated_significant_results |>
    dplyr::filter(
      significance == "Downregulated",
      !is.na(adjusted_p_value)
    ) |>
    dplyr::arrange(
      adjusted_p_value,
      dplyr::desc(
        abs(log2_fold_change)
      )
    ) |>
    dplyr::slice_head(n = 15)
) |>
  dplyr::distinct(
    gene_id,
    .keep_all = TRUE
  )

readr::write_csv(
  top_heatmap_genes,
  file.path(
    table_folder,
    "GSE254530_top30_annotated_heatmap_genes.csv"
  )
)

# ------------------------------------------------------------
# 14. Unmapped genes
# ------------------------------------------------------------

unmapped_genes <- annotated_all_results |>
  dplyr::filter(
    annotation_status == "Unmapped"
  ) |>
  dplyr::arrange(
    adjusted_p_value
  )

readr::write_csv(
  unmapped_genes,
  file.path(
    table_folder,
    "GSE254530_unmapped_gene_ids.csv"
  )
)

# ------------------------------------------------------------
# 15. Annotation coverage summaries
# ------------------------------------------------------------

annotation_summary <- annotated_all_results |>
  dplyr::group_by(
    annotation_status
  ) |>
  dplyr::summarise(
    number_of_genes = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    percentage = round(
      100 *
        number_of_genes /
        sum(number_of_genes),
      2
    )
  )

significant_annotation_summary <-
  annotated_significant_results |>
  dplyr::group_by(
    annotation_status
  ) |>
  dplyr::summarise(
    number_of_genes = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    percentage = round(
      100 *
        number_of_genes /
        sum(number_of_genes),
      2
    )
  )

readr::write_csv(
  annotation_summary,
  file.path(
    table_folder,
    "GSE254530_annotation_coverage_all_genes.csv"
  )
)

readr::write_csv(
  significant_annotation_summary,
  file.path(
    table_folder,
    "GSE254530_annotation_coverage_significant_genes.csv"
  )
)

# ------------------------------------------------------------
# 16. Final report
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 10 ANNOTATION COMPLETED\n")
cat("============================================\n")

cat(
  "Unique SAA6008 annotation records:",
  nrow(annotation_table),
  "\n"
)

cat(
  "All analysed genes:",
  nrow(annotated_all_results),
  "\n"
)

cat(
  "Significant DEGs:",
  nrow(annotated_significant_results),
  "\n"
)

cat(
  "Gene symbols available:",
  sum(
    annotated_all_results$
      annotation_status ==
      "Gene symbol available"
  ),
  "\n"
)

cat(
  "Product-only annotations:",
  sum(
    annotated_all_results$
      annotation_status ==
      "Product annotation only"
  ),
  "\n"
)

cat(
  "Unmapped genes:",
  nrow(unmapped_genes),
  "\n\n"
)

print(
  annotation_summary,
  n = Inf
)

cat(
  "\nFiles saved inside:\n",
  table_folder,
  "\n"
)