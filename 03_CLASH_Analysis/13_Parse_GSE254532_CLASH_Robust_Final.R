# ============================================================
# BacRegRNA Project
# Script 13 — ROBUST FINAL VERSION
#
# Parse and integrate GSE254532 RNase III-CLASH interactions
# with the JKD6008 master knowledgebase.
#
# Central objective:
# Reconstructing the sRNA regulatory network underlying the
# vancomycin response in VISA Staphylococcus aureus.
#
# Main outputs:
# 1. Standardized CLASH interaction table
# 2. Collapsed sRNA-mRNA edge table
# 3. Vancomycin-responsive CLASH targets
# 4. Regulatory RNA summary
# 5. Network-node table
# 6. CLASH-updated master knowledgebase
# 7. Mapping and reproducibility audit
# ============================================================


# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

clash_file <- file.path(
  project_folder,
  "07_processed_data",
  "downloaded_files",
  "GSE254532",
  "GSE254532_Final_RATT_VISACLASH.txt.gz"
)

master_kb_file <- file.path(
  project_folder,
  "12_master_knowledgebase",
  "tables",
  "JKD6008_master_gene_knowledgebase.csv"
)

output_folder <- file.path(
  project_folder,
  "13_CLASH_network"
)

table_folder <- file.path(
  output_folder,
  "tables"
)

audit_folder <- file.path(
  output_folder,
  "audit"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  audit_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

setwd(project_folder)


# ------------------------------------------------------------
# 2. Install packages
# ------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "tibble",
  "tidyr",
  "purrr"
)

installed_package_names <- rownames(
  installed.packages()
)

missing_packages <- setdiff(
  required_packages,
  installed_package_names
)

if (length(missing_packages) > 0) {
  
  install.packages(
    missing_packages
  )
}


# ------------------------------------------------------------
# 3. Load packages
# ------------------------------------------------------------

library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(purrr)


# ------------------------------------------------------------
# 4. Verify input files
# ------------------------------------------------------------

required_input_files <- c(
  clash_file,
  master_kb_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  
  stop(
    "Required input files were not found:\n",
    paste(
      missing_input_files,
      collapse = "\n"
    )
  )
}

cat(
  "\nInput files detected successfully.\n"
)

cat(
  "\nCLASH file:\n",
  clash_file,
  "\n"
)

cat(
  "\nMaster knowledgebase:\n",
  master_kb_file,
  "\n"
)


# ============================================================
# HELPER FUNCTIONS
# ============================================================


# ------------------------------------------------------------
# 5. Clean column names
# ------------------------------------------------------------

clean_column_names <- function(x) {
  
  cleaned <- x |>
    as.character() |>
    stringr::str_to_lower() |>
    stringr::str_replace_all(
      "[^a-z0-9]+",
      "_"
    ) |>
    stringr::str_replace_all(
      "^_+|_+$",
      ""
    )
  
  make.unique(
    cleaned,
    sep = "_"
  )
}


# ------------------------------------------------------------
# 6. Normalize identifiers
# ------------------------------------------------------------

normalize_identifier <- function(x) {
  
  x <- as.character(x)
  
  x <- stringr::str_squish(x)
  
  x <- stringr::str_remove_all(
    x,
    "^['\"]|['\"]$"
  )
  
  x <- stringr::str_replace_all(
    x,
    "[\\.\\-]",
    "_"
  )
  
  x <- stringr::str_replace_all(
    x,
    "\\s+",
    ""
  )
  
  x <- stringr::str_to_lower(x)
  
  x[
    is.na(x) |
      x %in% c(
        "",
        "na",
        "n_a",
        "nan",
        "null"
      )
  ] <- NA_character_
  
  x
}


# ------------------------------------------------------------
# 7. Safe numeric conversion
# ------------------------------------------------------------

safe_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


# ------------------------------------------------------------
# 8. Safe logical conversion
# ------------------------------------------------------------

safe_logical <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  x_character <- stringr::str_to_lower(
    stringr::str_squish(
      as.character(x)
    )
  )
  
  dplyr::case_when(
    
    x_character %in% c(
      "true",
      "t",
      "1",
      "yes",
      "y"
    ) ~ TRUE,
    
    x_character %in% c(
      "false",
      "f",
      "0",
      "no",
      "n"
    ) ~ FALSE,
    
    TRUE ~ FALSE
  )
}


# ------------------------------------------------------------
# 9. First non-missing character
# ------------------------------------------------------------

first_character_or_na <- function(x) {
  
  x <- as.character(x)
  
  keep <- !is.na(x) &
    stringr::str_squish(x) != ""
  
  x <- x[keep]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  stringr::str_squish(
    x[1]
  )
}


# ------------------------------------------------------------
# 10. First non-missing numeric
# ------------------------------------------------------------

first_numeric_or_na <- function(x) {
  
  x <- safe_numeric(x)
  
  x <- x[
    !is.na(x)
  ]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  x[1]
}


# ------------------------------------------------------------
# 11. Minimum numeric value
# ------------------------------------------------------------

minimum_or_na <- function(x) {
  
  x <- safe_numeric(x)
  
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  
  min(
    x,
    na.rm = TRUE
  )
}


# ------------------------------------------------------------
# 12. Maximum numeric value
# ------------------------------------------------------------

maximum_or_na <- function(x) {
  
  x <- safe_numeric(x)
  
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  
  max(
    x,
    na.rm = TRUE
  )
}


# ------------------------------------------------------------
# 13. Sum numeric values
# ------------------------------------------------------------

sum_or_na <- function(x) {
  
  x <- safe_numeric(x)
  
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  
  sum(
    x,
    na.rm = TRUE
  )
}


# ------------------------------------------------------------
# 14. Collapse unique values
# ------------------------------------------------------------

collapse_unique <- function(x) {
  
  x <- as.character(x)
  
  x <- stringr::str_squish(x)
  
  x <- unique(
    x[
      !is.na(x) &
        x != ""
    ]
  )
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(
    x,
    collapse = "; "
  )
}


# ------------------------------------------------------------
# 15. Add a missing column
# ------------------------------------------------------------

add_missing_column <- function(
    data,
    column_name,
    default_value
) {
  
  if (!column_name %in% names(data)) {
    
    data[[column_name]] <- rep(
      default_value,
      nrow(data)
    )
  }
  
  data
}


# ------------------------------------------------------------
# 16. Convert accidental list column to atomic character
# ------------------------------------------------------------

flatten_character_column <- function(x) {
  
  if (is.list(x)) {
    
    return(
      vapply(
        x,
        function(value) {
          
          if (
            is.null(value) ||
            length(value) == 0 ||
            all(is.na(value))
          ) {
            return(NA_character_)
          }
          
          as.character(
            value[1]
          )
        },
        FUN.VALUE = character(1)
      )
    )
  }
  
  as.character(x)
}


# ------------------------------------------------------------
# 17. Create empty data frames safely
# ------------------------------------------------------------

empty_srna_edge_table <- function() {
  
  tibble::tibble(
    
    source_node = character(),
    target_node = character(),
    source_label = character(),
    target_label = character(),
    source_raw_names = character(),
    target_raw_names = character(),
    source_type = character(),
    target_type = character(),
    interaction_class = character(),
    number_of_supporting_rows = integer(),
    total_hybrid_count = numeric(),
    best_raw_p_value = numeric(),
    best_adjusted_p_value = numeric(),
    best_connection_score = numeric(),
    maximum_number_of_experiments = numeric(),
    experimental_dataset = character(),
    evidence_type = character(),
    target_vancomycin_response = character(),
    target_log2_fold_change = numeric(),
    target_deseq_adjusted_p = numeric(),
    target_is_significant_deg = logical(),
    target_strong_response = logical(),
    evidence_tier = character()
  )
}


# ============================================================
# READ CLASH TABLE
# ============================================================


# ------------------------------------------------------------
# 18. Read compressed CLASH file
# ------------------------------------------------------------

raw_clash <- readr::read_tsv(
  
  file = clash_file,
  
  col_types = readr::cols(
    .default = readr::col_character()
  ),
  
  trim_ws = TRUE,
  
  na = c(
    "",
    "NA",
    "N/A",
    "NaN",
    "null",
    "NULL"
  ),
  
  show_col_types = FALSE,
  
  progress = FALSE
)

if (!is.data.frame(raw_clash)) {
  
  stop(
    "The CLASH object is not a data frame. Class: ",
    paste(
      class(raw_clash),
      collapse = ", "
    )
  )
}

if (nrow(raw_clash) == 0) {
  
  stop(
    "The CLASH file was read successfully but contains zero rows."
  )
}

original_clash_column_names <- names(
  raw_clash
)

names(raw_clash) <- clean_column_names(
  original_clash_column_names
)

cat(
  "\nCLASH dimensions:",
  nrow(raw_clash),
  "rows x",
  ncol(raw_clash),
  "columns\n"
)

cat(
  "\nStandardized CLASH columns:\n"
)

print(
  names(raw_clash)
)


# ------------------------------------------------------------
# 19. Essential CLASH columns
# ------------------------------------------------------------

essential_clash_columns <- c(
  "left_name",
  "locus_tag",
  "common_name",
  "rna_class",
  "name_1",
  "locus_tag_1",
  "common_name_1",
  "rna_class_1"
)

missing_essential_clash_columns <- setdiff(
  essential_clash_columns,
  names(raw_clash)
)

if (
  length(
    missing_essential_clash_columns
  ) > 0
) {
  
  stop(
    "Essential CLASH columns were not found:\n",
    paste(
      missing_essential_clash_columns,
      collapse = "\n"
    ),
    "\n\nAvailable columns:\n",
    paste(
      names(raw_clash),
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------
# 20. Add optional CLASH columns
# ------------------------------------------------------------

optional_clash_columns <- list(
  
  hyb_count = NA_character_,
  hyb_count_1 = NA_character_,
  
  total_hybrids = NA_character_,
  
  p_value = NA_character_,
  bh_adj_p_value = NA_character_,
  
  connection_score = NA_character_,
  number_expts = NA_character_,
  
  chromo = NA_character_,
  start = NA_character_,
  end = NA_character_,
  strand = NA_character_,
  
  chromo_1 = NA_character_,
  start_1 = NA_character_,
  end_1 = NA_character_,
  strand_1 = NA_character_
)

for (
  column_name in
  names(optional_clash_columns)
) {
  
  default_value <-
    optional_clash_columns[[column_name]]
  
  raw_clash <- add_missing_column(
    raw_clash,
    column_name,
    default_value
  )
}


# ------------------------------------------------------------
# 21. Save CLASH column audit
# ------------------------------------------------------------

clash_column_dictionary <- tibble::tibble(
  
  original_column =
    original_clash_column_names,
  
  standardized_column =
    names(raw_clash)[
      seq_along(
        original_clash_column_names
      )
    ]
)

readr::write_csv(
  clash_column_dictionary,
  file.path(
    audit_folder,
    "GSE254532_column_dictionary.csv"
  )
)

readr::write_csv(
  utils::head(
    raw_clash,
    100
  ),
  file.path(
    audit_folder,
    "GSE254532_first_100_rows.csv"
  )
)


# ============================================================
# READ MASTER KNOWLEDGEBASE
# ============================================================


# ------------------------------------------------------------
# 22. Read master table
# ------------------------------------------------------------

master_kb <- readr::read_csv(
  master_kb_file,
  show_col_types = FALSE
)

if (!is.data.frame(master_kb)) {
  
  stop(
    "master_kb is not a data frame. Class: ",
    paste(
      class(master_kb),
      collapse = ", "
    )
  )
}

if (nrow(master_kb) == 0) {
  
  stop(
    "The master knowledgebase contains zero rows."
  )
}


# ------------------------------------------------------------
# 23. Required master columns and defaults
# ------------------------------------------------------------

required_master_columns <- list(
  
  gene_id = NA_character_,
  current_locus_tag = NA_character_,
  gene_symbol = NA_character_,
  gene_synonym = NA_character_,
  master_label = NA_character_,
  product = NA_character_,
  
  is_rna_feature = FALSE,
  
  vancomycin_response = NA_character_,
  log2_fold_change = NA_real_,
  adjusted_p_value = NA_real_,
  
  is_significant_deg = FALSE,
  strong_response = FALSE,
  
  clash_target_status = FALSE,
  clash_srna_status = FALSE,
  
  number_of_clash_srna_regulators = 0L,
  number_of_clash_mrna_targets = 0L,
  
  clash_srna_regulators = NA_character_,
  clash_mrna_targets = NA_character_,
  
  clash_total_hybrid_count = NA_real_,
  clash_best_adjusted_p = NA_real_,
  clash_best_connection_score = NA_real_,
  
  clash_evidence_level = NA_character_,
  
  srna_network_role = NA_character_,
  interaction_evidence = NA_character_,
  evidence_tier = NA_character_
)


# ------------------------------------------------------------
# 24. Add missing master columns
#
# Important corrected syntax:
# required_master_columns[[column_name]]
# ------------------------------------------------------------

for (
  column_name in
  names(required_master_columns)
) {
  
  default_value <-
    required_master_columns[[column_name]]
  
  master_kb <- add_missing_column(
    master_kb,
    column_name,
    default_value
  )
}


# ------------------------------------------------------------
# 25. Flatten and standardize master columns
# ------------------------------------------------------------

master_kb$gene_id <- flatten_character_column(
  master_kb$gene_id
)

master_kb$current_locus_tag <-
  flatten_character_column(
    master_kb$current_locus_tag
  )

master_kb$gene_symbol <-
  flatten_character_column(
    master_kb$gene_symbol
  )

master_kb$gene_synonym <-
  flatten_character_column(
    master_kb$gene_synonym
  )

master_kb$master_label <-
  flatten_character_column(
    master_kb$master_label
  )

master_kb$product <-
  flatten_character_column(
    master_kb$product
  )

master_kb$vancomycin_response <-
  flatten_character_column(
    master_kb$vancomycin_response
  )

master_kb$is_rna_feature <-
  safe_logical(
    master_kb$is_rna_feature
  )

master_kb$is_significant_deg <-
  safe_logical(
    master_kb$is_significant_deg
  )

master_kb$strong_response <-
  safe_logical(
    master_kb$strong_response
  )

master_kb$clash_target_status <-
  safe_logical(
    master_kb$clash_target_status
  )

master_kb$clash_srna_status <-
  safe_logical(
    master_kb$clash_srna_status
  )

master_kb$log2_fold_change <-
  safe_numeric(
    master_kb$log2_fold_change
  )

master_kb$adjusted_p_value <-
  safe_numeric(
    master_kb$adjusted_p_value
  )

master_kb$clash_total_hybrid_count <-
  safe_numeric(
    master_kb$clash_total_hybrid_count
  )

master_kb$clash_best_adjusted_p <-
  safe_numeric(
    master_kb$clash_best_adjusted_p
  )

master_kb$clash_best_connection_score <-
  safe_numeric(
    master_kb$clash_best_connection_score
  )

master_kb$number_of_clash_srna_regulators <-
  as.integer(
    safe_numeric(
      master_kb$
        number_of_clash_srna_regulators
    )
  )

master_kb$number_of_clash_mrna_targets <-
  as.integer(
    safe_numeric(
      master_kb$
        number_of_clash_mrna_targets
    )
  )

master_kb$number_of_clash_srna_regulators[
  is.na(
    master_kb$
      number_of_clash_srna_regulators
  )
] <- 0L

master_kb$number_of_clash_mrna_targets[
  is.na(
    master_kb$
      number_of_clash_mrna_targets
  )
] <- 0L


# ------------------------------------------------------------
# 26. Validate gene_id
# ------------------------------------------------------------

if (!"gene_id" %in% names(master_kb)) {
  
  stop(
    "The master knowledgebase does not contain gene_id."
  )
}

if (
  length(master_kb$gene_id) !=
  nrow(master_kb)
) {
  
  stop(
    "gene_id length does not match the number of master rows."
  )
}


# ------------------------------------------------------------
# 27. Duplicate gene-ID check without count()
#
# Uses base R table() to avoid package conflicts.
# ------------------------------------------------------------

valid_master_gene_ids <- master_kb$gene_id[
  !is.na(master_kb$gene_id) &
    stringr::str_squish(
      master_kb$gene_id
    ) != ""
]

gene_id_frequency <- table(
  valid_master_gene_ids,
  useNA = "no"
)

duplicate_gene_id_frequency <-
  gene_id_frequency[
    gene_id_frequency > 1
  ]

duplicate_master_ids <- tibble::tibble(
  
  gene_id =
    names(
      duplicate_gene_id_frequency
    ),
  
  number_of_rows =
    as.integer(
      duplicate_gene_id_frequency
    )
) |>
  dplyr::arrange(
    dplyr::desc(
      .data$number_of_rows
    )
  )

if (nrow(duplicate_master_ids) > 0) {
  
  readr::write_csv(
    duplicate_master_ids,
    file.path(
      audit_folder,
      "duplicate_master_gene_ids.csv"
    )
  )
  
  stop(
    "Duplicated gene_id values were found. See:\n",
    file.path(
      audit_folder,
      "duplicate_master_gene_ids.csv"
    )
  )
}

cat(
  "\nMaster gene-ID check passed.\n"
)

cat(
  "Unique master gene IDs:",
  length(
    unique(
      valid_master_gene_ids
    )
  ),
  "\n"
)


# ============================================================
# IDENTIFIER LOOKUP
# ============================================================


# ------------------------------------------------------------
# 28. Build long identifier lookup
# ------------------------------------------------------------

identifier_lookup <- dplyr::bind_rows(
  
  master_kb |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      identifier_type = "gene_id",
      identifier_value = .data$gene_id
    ),
  
  master_kb |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      identifier_type =
        "current_locus_tag",
      identifier_value =
        .data$current_locus_tag
    ),
  
  master_kb |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      identifier_type = "gene_symbol",
      identifier_value =
        .data$gene_symbol
    ),
  
  master_kb |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      identifier_type = "gene_synonym",
      identifier_value =
        .data$gene_synonym
    ),
  
  master_kb |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      identifier_type = "master_label",
      identifier_value =
        .data$master_label
    )
) |>
  dplyr::filter(
    !is.na(.data$identifier_value),
    stringr::str_squish(
      .data$identifier_value
    ) != ""
  ) |>
  dplyr::mutate(
    identifier_normalized =
      normalize_identifier(
        .data$identifier_value
      )
  ) |>
  dplyr::filter(
    !is.na(
      .data$identifier_normalized
    ),
    .data$identifier_normalized != ""
  ) |>
  dplyr::distinct(
    .data$identifier_normalized,
    .data$gene_id,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# 29. Remove ambiguous identifiers
#
# A symbol mapping to more than one gene is excluded.
# ------------------------------------------------------------

identifier_frequency <- identifier_lookup |>
  dplyr::group_by(
    .data$identifier_normalized
  ) |>
  dplyr::summarise(
    number_of_gene_ids =
      dplyr::n_distinct(
        .data$gene_id
      ),
    .groups = "drop"
  )

ambiguous_identifiers <-
  identifier_frequency |>
  dplyr::filter(
    .data$number_of_gene_ids > 1
  )

if (nrow(ambiguous_identifiers) > 0) {
  
  readr::write_csv(
    ambiguous_identifiers,
    file.path(
      audit_folder,
      "ambiguous_master_identifiers.csv"
    )
  )
}

identifier_lookup <- identifier_lookup |>
  dplyr::anti_join(
    ambiguous_identifiers,
    by = "identifier_normalized"
  ) |>
  dplyr::distinct(
    .data$identifier_normalized,
    .keep_all = TRUE
  )


# ============================================================
# STANDARDIZE CLASH INTERACTIONS
# ============================================================


# ------------------------------------------------------------
# 30. Build standardized endpoint table
# ------------------------------------------------------------

clash_standardized <- raw_clash |>
  dplyr::mutate(
    
    interaction_row_id =
      dplyr::row_number(),
    
    endpoint1_feature_name =
      dplyr::na_if(
        stringr::str_squish(
          as.character(
            .data$left_name
          )
        ),
        ""
      ),
    
    endpoint2_feature_name =
      dplyr::na_if(
        stringr::str_squish(
          as.character(
            .data$name_1
          )
        ),
        ""
      ),
    
    endpoint1_locus_raw =
      dplyr::na_if(
        stringr::str_squish(
          as.character(
            .data$locus_tag
          )
        ),
        ""
      ),
    
    endpoint2_locus_raw =
      dplyr::na_if(
        stringr::str_squish(
          as.character(
            .data$locus_tag_1
          )
        ),
        ""
      ),
    
    endpoint1_common_name =
      dplyr::na_if(
        stringr::str_squish(
          as.character(
            .data$common_name
          )
        ),
        ""
      ),
    
    endpoint2_common_name =
      dplyr::na_if(
        stringr::str_squish(
          as.character(
            .data$common_name_1
          )
        ),
        ""
      ),
    
    endpoint1_class_raw =
      dplyr::na_if(
        stringr::str_squish(
          as.character(
            .data$rna_class
          )
        ),
        ""
      ),
    
    endpoint2_class_raw =
      dplyr::na_if(
        stringr::str_squish(
          as.character(
            .data$rna_class_1
          )
        ),
        ""
      ),
    
    endpoint1_name =
      dplyr::coalesce(
        .data$endpoint1_feature_name,
        .data$endpoint1_common_name,
        .data$endpoint1_locus_raw
      ),
    
    endpoint2_name =
      dplyr::coalesce(
        .data$endpoint2_feature_name,
        .data$endpoint2_common_name,
        .data$endpoint2_locus_raw
      ),
    
    endpoint1_hybrid_count =
      safe_numeric(
        .data$hyb_count
      ),
    
    endpoint2_hybrid_count =
      safe_numeric(
        .data$hyb_count_1
      ),
    
    total_hybrid_count =
      safe_numeric(
        .data$total_hybrids
      ),
    
    raw_p_value =
      safe_numeric(
        .data$p_value
      ),
    
    clash_adjusted_p_value =
      safe_numeric(
        .data$bh_adj_p_value
      ),
    
    connection_score_numeric =
      safe_numeric(
        .data$connection_score
      ),
    
    number_of_experiments =
      safe_numeric(
        .data$number_expts
      ),
    
    endpoint1_start =
      safe_numeric(
        .data$start
      ),
    
    endpoint1_end =
      safe_numeric(
        .data$end
      ),
    
    endpoint2_start =
      safe_numeric(
        .data$start_1
      ),
    
    endpoint2_end =
      safe_numeric(
        .data$end_1
      )
  ) |>
  dplyr::filter(
    !is.na(.data$endpoint1_name),
    !is.na(.data$endpoint2_name)
  )

if (nrow(clash_standardized) == 0) {
  
  stop(
    "No valid two-endpoint CLASH interactions remained."
  )
}


# ============================================================
# MAP ENDPOINTS
# ============================================================


# ------------------------------------------------------------
# 31. Map one endpoint
# ------------------------------------------------------------

map_endpoint <- function(
    endpoint_name,
    endpoint_locus,
    endpoint_common_name,
    identifier_lookup_table
) {
  
  endpoint_table <- tibble::tibble(
    
    endpoint_row =
      seq_along(
        endpoint_name
      ),
    
    endpoint_name =
      as.character(
        endpoint_name
      ),
    
    endpoint_locus =
      as.character(
        endpoint_locus
      ),
    
    endpoint_common_name =
      as.character(
        endpoint_common_name
      )
  ) |>
    dplyr::mutate(
      
      locus_normalized =
        normalize_identifier(
          .data$endpoint_locus
        ),
      
      common_normalized =
        normalize_identifier(
          .data$endpoint_common_name
        ),
      
      name_normalized =
        normalize_identifier(
          .data$endpoint_name
        )
    )
  
  locus_lookup <- identifier_lookup_table |>
    dplyr::transmute(
      
      locus_normalized =
        .data$identifier_normalized,
      
      locus_gene_id =
        .data$gene_id,
      
      locus_identifier_type =
        .data$identifier_type
    )
  
  common_lookup <- identifier_lookup_table |>
    dplyr::transmute(
      
      common_normalized =
        .data$identifier_normalized,
      
      common_gene_id =
        .data$gene_id,
      
      common_identifier_type =
        .data$identifier_type
    )
  
  name_lookup <- identifier_lookup_table |>
    dplyr::transmute(
      
      name_normalized =
        .data$identifier_normalized,
      
      name_gene_id =
        .data$gene_id,
      
      name_identifier_type =
        .data$identifier_type
    )
  
  endpoint_table |>
    dplyr::left_join(
      locus_lookup,
      by = "locus_normalized"
    ) |>
    dplyr::left_join(
      common_lookup,
      by = "common_normalized"
    ) |>
    dplyr::left_join(
      name_lookup,
      by = "name_normalized"
    ) |>
    dplyr::transmute(
      
      endpoint_row =
        .data$endpoint_row,
      
      mapped_gene_id =
        dplyr::coalesce(
          .data$locus_gene_id,
          .data$common_gene_id,
          .data$name_gene_id
        ),
      
      mapping_method =
        dplyr::case_when(
          
          !is.na(
            .data$locus_gene_id
          ) ~
            paste0(
              "Locus match: ",
              .data$locus_identifier_type
            ),
          
          !is.na(
            .data$common_gene_id
          ) ~
            paste0(
              "Common-name match: ",
              .data$common_identifier_type
            ),
          
          !is.na(
            .data$name_gene_id
          ) ~
            paste0(
              "Feature-name match: ",
              .data$name_identifier_type
            ),
          
          TRUE ~ "Unmapped"
        )
    )
}


# ------------------------------------------------------------
# 32. Map endpoint 1
# ------------------------------------------------------------

endpoint1_map <- map_endpoint(
  
  endpoint_name =
    clash_standardized$
    endpoint1_name,
  
  endpoint_locus =
    clash_standardized$
    endpoint1_locus_raw,
  
  endpoint_common_name =
    clash_standardized$
    endpoint1_common_name,
  
  identifier_lookup_table =
    identifier_lookup
) |>
  dplyr::rename(
    
    endpoint1_mapped_gene_id =
      .data$mapped_gene_id,
    
    endpoint1_mapping_method =
      .data$mapping_method
  )


# ------------------------------------------------------------
# 33. Map endpoint 2
# ------------------------------------------------------------

endpoint2_map <- map_endpoint(
  
  endpoint_name =
    clash_standardized$
    endpoint2_name,
  
  endpoint_locus =
    clash_standardized$
    endpoint2_locus_raw,
  
  endpoint_common_name =
    clash_standardized$
    endpoint2_common_name,
  
  identifier_lookup_table =
    identifier_lookup
) |>
  dplyr::rename(
    
    endpoint2_mapped_gene_id =
      .data$mapped_gene_id,
    
    endpoint2_mapping_method =
      .data$mapping_method
  )


# ------------------------------------------------------------
# 34. Join endpoint mappings
# ------------------------------------------------------------

clash_mapped <- clash_standardized |>
  dplyr::mutate(
    endpoint_row =
      dplyr::row_number()
  ) |>
  dplyr::left_join(
    endpoint1_map,
    by = "endpoint_row"
  ) |>
  dplyr::left_join(
    endpoint2_map,
    by = "endpoint_row"
  ) |>
  dplyr::select(
    -.data$endpoint_row
  )


# ============================================================
# ADD MASTER ANNOTATIONS
# ============================================================


# ------------------------------------------------------------
# 35. Endpoint 1 annotation
# ------------------------------------------------------------

endpoint1_annotation <- master_kb |>
  dplyr::transmute(
    
    endpoint1_mapped_gene_id =
      .data$gene_id,
    
    endpoint1_gene_symbol =
      .data$gene_symbol,
    
    endpoint1_master_label =
      .data$master_label,
    
    endpoint1_product =
      .data$product,
    
    endpoint1_is_rna_feature =
      .data$is_rna_feature,
    
    endpoint1_vancomycin_response =
      .data$vancomycin_response,
    
    endpoint1_log2_fold_change =
      .data$log2_fold_change,
    
    endpoint1_deseq_adjusted_p =
      .data$adjusted_p_value,
    
    endpoint1_is_significant_deg =
      .data$is_significant_deg,
    
    endpoint1_strong_response =
      .data$strong_response
  )


# ------------------------------------------------------------
# 36. Endpoint 2 annotation
# ------------------------------------------------------------

endpoint2_annotation <- master_kb |>
  dplyr::transmute(
    
    endpoint2_mapped_gene_id =
      .data$gene_id,
    
    endpoint2_gene_symbol =
      .data$gene_symbol,
    
    endpoint2_master_label =
      .data$master_label,
    
    endpoint2_product =
      .data$product,
    
    endpoint2_is_rna_feature =
      .data$is_rna_feature,
    
    endpoint2_vancomycin_response =
      .data$vancomycin_response,
    
    endpoint2_log2_fold_change =
      .data$log2_fold_change,
    
    endpoint2_deseq_adjusted_p =
      .data$adjusted_p_value,
    
    endpoint2_is_significant_deg =
      .data$is_significant_deg,
    
    endpoint2_strong_response =
      .data$strong_response
  )


# ------------------------------------------------------------
# 37. Attach annotations
# ------------------------------------------------------------

clash_mapped <- clash_mapped |>
  dplyr::left_join(
    endpoint1_annotation,
    by = "endpoint1_mapped_gene_id"
  ) |>
  dplyr::left_join(
    endpoint2_annotation,
    by = "endpoint2_mapped_gene_id"
  )


# ============================================================
# CLASSIFY ENDPOINTS
# ============================================================


# ------------------------------------------------------------
# 38. Endpoint classification function
# ------------------------------------------------------------

classify_endpoint <- function(
    raw_class,
    feature_name,
    common_name,
    mapped_gene_id,
    is_rna_feature
) {
  
  class_text <- stringr::str_to_lower(
    dplyr::coalesce(
      as.character(raw_class),
      ""
    )
  )
  
  name_text <- paste(
    as.character(feature_name),
    as.character(common_name)
  )
  
  dplyr::case_when(
    
    stringr::str_detect(
      class_text,
      "srna|small.?rna|ncrna|antisense|asrna"
    ) ~ "Regulatory RNA",
    
    stringr::str_detect(
      name_text,
      stringr::regex(
        paste0(
          "Rsa[A-Za-z0-9/]+|",
          "Spr[A-Za-z0-9/]+|",
          "RNAIII|",
          "sRNA[-_][0-9]+|",
          "vigR|",
          "3.?UTR|",
          "5.?UTR"
        ),
        ignore_case = TRUE
      )
    ) ~ "Regulatory RNA",
    
    stringr::str_detect(
      class_text,
      "rrna"
    ) ~ "rRNA",
    
    stringr::str_detect(
      class_text,
      "trna"
    ) ~ "tRNA",
    
    stringr::str_detect(
      class_text,
      "intergenic"
    ) ~ "Intergenic RNA",
    
    stringr::str_detect(
      class_text,
      "utr"
    ) ~ "Regulatory RNA",
    
    stringr::str_detect(
      class_text,
      "cds|mrna|coding"
    ) ~ "mRNA",
    
    isTRUE(is_rna_feature) ~
      "Regulatory RNA",
    
    !is.na(mapped_gene_id) ~
      "mRNA",
    
    TRUE ~
      "Unclassified RNA"
  )
}


# ------------------------------------------------------------
# 39. Apply endpoint classifications
# ------------------------------------------------------------

clash_mapped <- clash_mapped |>
  dplyr::rowwise() |>
  dplyr::mutate(
    
    endpoint1_type =
      classify_endpoint(
        
        raw_class =
          .data$endpoint1_class_raw,
        
        feature_name =
          .data$endpoint1_feature_name,
        
        common_name =
          .data$endpoint1_common_name,
        
        mapped_gene_id =
          .data$endpoint1_mapped_gene_id,
        
        is_rna_feature =
          .data$endpoint1_is_rna_feature
      ),
    
    endpoint2_type =
      classify_endpoint(
        
        raw_class =
          .data$endpoint2_class_raw,
        
        feature_name =
          .data$endpoint2_feature_name,
        
        common_name =
          .data$endpoint2_common_name,
        
        mapped_gene_id =
          .data$endpoint2_mapped_gene_id,
        
        is_rna_feature =
          .data$endpoint2_is_rna_feature
      )
  ) |>
  dplyr::ungroup()


# ------------------------------------------------------------
# 40. Interaction classification
# ------------------------------------------------------------

clash_mapped <- clash_mapped |>
  dplyr::mutate(
    
    interaction_class =
      dplyr::case_when(
        
        .data$endpoint1_type ==
          "Regulatory RNA" &
          .data$endpoint2_type ==
          "mRNA" ~
          "sRNA-mRNA",
        
        .data$endpoint1_type ==
          "mRNA" &
          .data$endpoint2_type ==
          "Regulatory RNA" ~
          "sRNA-mRNA",
        
        .data$endpoint1_type ==
          "Intergenic RNA" &
          .data$endpoint2_type ==
          "mRNA" ~
          "candidate_sRNA-mRNA",
        
        .data$endpoint1_type ==
          "mRNA" &
          .data$endpoint2_type ==
          "Intergenic RNA" ~
          "candidate_sRNA-mRNA",
        
        .data$endpoint1_type ==
          "mRNA" &
          .data$endpoint2_type ==
          "mRNA" ~
          "mRNA-mRNA",
        
        .data$endpoint1_type ==
          "Regulatory RNA" &
          .data$endpoint2_type ==
          "Regulatory RNA" ~
          "sRNA-sRNA",
        
        .data$endpoint1_type %in% c(
          "rRNA",
          "tRNA"
        ) |
          .data$endpoint2_type %in% c(
            "rRNA",
            "tRNA"
          ) ~
          "Structural RNA-associated",
        
        TRUE ~
          "Other or unclassified"
      ),
    
    both_endpoints_mapped =
      !is.na(
        .data$endpoint1_mapped_gene_id
      ) &
      !is.na(
        .data$endpoint2_mapped_gene_id
      ),
    
    at_least_one_endpoint_mapped =
      !is.na(
        .data$endpoint1_mapped_gene_id
      ) |
      !is.na(
        .data$endpoint2_mapped_gene_id
      )
  )


# ============================================================
# CONSTRUCT sRNA-mRNA EDGES
# ============================================================


# ------------------------------------------------------------
# 41. Create stable unmapped RNA node IDs
# ------------------------------------------------------------

make_unmapped_node_id <- function(
    endpoint_name,
    endpoint_class,
    chromosome,
    start_position,
    end_position,
    strand_value
) {
  
  node_name <- normalize_identifier(
    endpoint_name
  )
  
  node_class <- normalize_identifier(
    endpoint_class
  )
  
  chromosome <- normalize_identifier(
    chromosome
  )
  
  node_name <- dplyr::coalesce(
    node_name,
    "unnamed"
  )
  
  node_class <- dplyr::coalesce(
    node_class,
    "unclassified"
  )
  
  chromosome <- dplyr::coalesce(
    chromosome,
    "unknown_chromosome"
  )
  
  paste0(
    "RNA:",
    node_name,
    ":",
    node_class,
    ":",
    chromosome,
    ":",
    dplyr::coalesce(
      as.character(
        start_position
      ),
      "NA"
    ),
    "-",
    dplyr::coalesce(
      as.character(
        end_position
      ),
      "NA"
    ),
    ":",
    dplyr::coalesce(
      as.character(
        strand_value
      ),
      "NA"
    )
  )
}


# ------------------------------------------------------------
# 42. Select regulatory RNA-mRNA interactions
# ------------------------------------------------------------

srna_candidate_rows <- clash_mapped |>
  dplyr::filter(
    .data$interaction_class %in% c(
      "sRNA-mRNA",
      "candidate_sRNA-mRNA"
    )
  )


# ------------------------------------------------------------
# 43. Orient regulatory RNA toward mRNA target
# ------------------------------------------------------------

if (nrow(srna_candidate_rows) > 0) {
  
  srna_mrna_rows <- srna_candidate_rows |>
    dplyr::mutate(
      
      source_is_endpoint1 =
        .data$endpoint1_type %in% c(
          "Regulatory RNA",
          "Intergenic RNA"
        ),
      
      source_raw_name =
        dplyr::if_else(
          .data$source_is_endpoint1,
          .data$endpoint1_name,
          .data$endpoint2_name
        ),
      
      target_raw_name =
        dplyr::if_else(
          .data$source_is_endpoint1,
          .data$endpoint2_name,
          .data$endpoint1_name
        ),
      
      source_gene_id =
        dplyr::if_else(
          .data$source_is_endpoint1,
          .data$endpoint1_mapped_gene_id,
          .data$endpoint2_mapped_gene_id
        ),
      
      target_gene_id =
        dplyr::if_else(
          .data$source_is_endpoint1,
          .data$endpoint2_mapped_gene_id,
          .data$endpoint1_mapped_gene_id
        ),
      
      source_label =
        dplyr::if_else(
          
          .data$source_is_endpoint1,
          
          dplyr::coalesce(
            .data$endpoint1_gene_symbol,
            .data$endpoint1_master_label,
            .data$endpoint1_common_name,
            .data$endpoint1_feature_name,
            .data$endpoint1_name
          ),
          
          dplyr::coalesce(
            .data$endpoint2_gene_symbol,
            .data$endpoint2_master_label,
            .data$endpoint2_common_name,
            .data$endpoint2_feature_name,
            .data$endpoint2_name
          )
        ),
      
      target_label =
        dplyr::if_else(
          
          .data$source_is_endpoint1,
          
          dplyr::coalesce(
            .data$endpoint2_gene_symbol,
            .data$endpoint2_master_label,
            .data$endpoint2_common_name,
            .data$endpoint2_feature_name,
            .data$endpoint2_name
          ),
          
          dplyr::coalesce(
            .data$endpoint1_gene_symbol,
            .data$endpoint1_master_label,
            .data$endpoint1_common_name,
            .data$endpoint1_feature_name,
            .data$endpoint1_name
          )
        ),
      
      target_vancomycin_response =
        dplyr::if_else(
          .data$source_is_endpoint1,
          .data$endpoint2_vancomycin_response,
          .data$endpoint1_vancomycin_response
        ),
      
      target_log2_fold_change =
        dplyr::if_else(
          .data$source_is_endpoint1,
          .data$endpoint2_log2_fold_change,
          .data$endpoint1_log2_fold_change
        ),
      
      target_deseq_adjusted_p =
        dplyr::if_else(
          .data$source_is_endpoint1,
          .data$endpoint2_deseq_adjusted_p,
          .data$endpoint1_deseq_adjusted_p
        ),
      
      target_is_significant_deg =
        dplyr::if_else(
          .data$source_is_endpoint1,
          .data$endpoint2_is_significant_deg,
          .data$endpoint1_is_significant_deg,
          missing = FALSE
        ),
      
      target_strong_response =
        dplyr::if_else(
          .data$source_is_endpoint1,
          .data$endpoint2_strong_response,
          .data$endpoint1_strong_response,
          missing = FALSE
        ),
      
      source_node =
        dplyr::if_else(
          
          !is.na(
            .data$source_gene_id
          ),
          
          .data$source_gene_id,
          
          dplyr::if_else(
            
            .data$source_is_endpoint1,
            
            make_unmapped_node_id(
              .data$endpoint1_name,
              .data$endpoint1_class_raw,
              .data$chromo,
              .data$endpoint1_start,
              .data$endpoint1_end,
              .data$strand
            ),
            
            make_unmapped_node_id(
              .data$endpoint2_name,
              .data$endpoint2_class_raw,
              .data$chromo_1,
              .data$endpoint2_start,
              .data$endpoint2_end,
              .data$strand_1
            )
          )
        ),
      
      target_node =
        dplyr::if_else(
          
          !is.na(
            .data$target_gene_id
          ),
          
          .data$target_gene_id,
          
          dplyr::if_else(
            
            .data$source_is_endpoint1,
            
            make_unmapped_node_id(
              .data$endpoint2_name,
              .data$endpoint2_class_raw,
              .data$chromo_1,
              .data$endpoint2_start,
              .data$endpoint2_end,
              .data$strand_1
            ),
            
            make_unmapped_node_id(
              .data$endpoint1_name,
              .data$endpoint1_class_raw,
              .data$chromo,
              .data$endpoint1_start,
              .data$endpoint1_end,
              .data$strand
            )
          )
        )
    ) |>
    dplyr::transmute(
      
      interaction_row_id =
        .data$interaction_row_id,
      
      source_node =
        .data$source_node,
      
      target_node =
        .data$target_node,
      
      source_label =
        .data$source_label,
      
      target_label =
        .data$target_label,
      
      source_raw_name =
        .data$source_raw_name,
      
      target_raw_name =
        .data$target_raw_name,
      
      source_type =
        "Regulatory RNA",
      
      target_type =
        "mRNA",
      
      interaction_class =
        .data$interaction_class,
      
      total_hybrid_count =
        .data$total_hybrid_count,
      
      endpoint1_hybrid_count =
        .data$endpoint1_hybrid_count,
      
      endpoint2_hybrid_count =
        .data$endpoint2_hybrid_count,
      
      raw_p_value =
        .data$raw_p_value,
      
      clash_adjusted_p_value =
        .data$clash_adjusted_p_value,
      
      connection_score =
        .data$connection_score_numeric,
      
      number_of_experiments =
        .data$number_of_experiments,
      
      experimental_dataset =
        "GSE254532",
      
      evidence_type =
        "Experimental RNase III-CLASH",
      
      target_vancomycin_response =
        .data$target_vancomycin_response,
      
      target_log2_fold_change =
        .data$target_log2_fold_change,
      
      target_deseq_adjusted_p =
        .data$target_deseq_adjusted_p,
      
      target_is_significant_deg =
        .data$target_is_significant_deg,
      
      target_strong_response =
        .data$target_strong_response
    )
  
} else {
  
  srna_mrna_rows <- tibble::tibble(
    
    interaction_row_id = integer(),
    source_node = character(),
    target_node = character(),
    source_label = character(),
    target_label = character(),
    source_raw_name = character(),
    target_raw_name = character(),
    source_type = character(),
    target_type = character(),
    interaction_class = character(),
    total_hybrid_count = numeric(),
    endpoint1_hybrid_count = numeric(),
    endpoint2_hybrid_count = numeric(),
    raw_p_value = numeric(),
    clash_adjusted_p_value = numeric(),
    connection_score = numeric(),
    number_of_experiments = numeric(),
    experimental_dataset = character(),
    evidence_type = character(),
    target_vancomycin_response = character(),
    target_log2_fold_change = numeric(),
    target_deseq_adjusted_p = numeric(),
    target_is_significant_deg = logical(),
    target_strong_response = logical()
  )
}


# ------------------------------------------------------------
# 44. Collapse repeated sRNA-mRNA edges
# ------------------------------------------------------------

if (nrow(srna_mrna_rows) > 0) {
  
  srna_mrna_edges <- srna_mrna_rows |>
    dplyr::group_by(
      .data$source_node,
      .data$target_node
    ) |>
    dplyr::summarise(
      
      source_label =
        first_character_or_na(
          .data$source_label
        ),
      
      target_label =
        first_character_or_na(
          .data$target_label
        ),
      
      source_raw_names =
        collapse_unique(
          .data$source_raw_name
        ),
      
      target_raw_names =
        collapse_unique(
          .data$target_raw_name
        ),
      
      source_type =
        "Regulatory RNA",
      
      target_type =
        "mRNA",
      
      interaction_class =
        collapse_unique(
          .data$interaction_class
        ),
      
      number_of_supporting_rows =
        dplyr::n(),
      
      total_hybrid_count =
        sum_or_na(
          .data$total_hybrid_count
        ),
      
      best_raw_p_value =
        minimum_or_na(
          .data$raw_p_value
        ),
      
      best_adjusted_p_value =
        minimum_or_na(
          .data$clash_adjusted_p_value
        ),
      
      best_connection_score =
        maximum_or_na(
          .data$connection_score
        ),
      
      maximum_number_of_experiments =
        maximum_or_na(
          .data$number_of_experiments
        ),
      
      experimental_dataset =
        "GSE254532",
      
      evidence_type =
        "Experimental RNase III-CLASH",
      
      target_vancomycin_response =
        first_character_or_na(
          .data$target_vancomycin_response
        ),
      
      target_log2_fold_change =
        first_numeric_or_na(
          .data$target_log2_fold_change
        ),
      
      target_deseq_adjusted_p =
        minimum_or_na(
          .data$target_deseq_adjusted_p
        ),
      
      target_is_significant_deg =
        any(
          .data$target_is_significant_deg %in%
            TRUE
        ),
      
      target_strong_response =
        any(
          .data$target_strong_response %in%
            TRUE
        ),
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      
      evidence_tier =
        dplyr::case_when(
          
          .data$target_strong_response ~
            paste(
              "Tier 2-high:",
              "CLASH plus strong",
              "vancomycin response"
            ),
          
          .data$target_is_significant_deg ~
            paste(
              "Tier 2:",
              "CLASH plus significant",
              "vancomycin response"
            ),
          
          !stringr::str_detect(
            .data$target_node,
            "^RNA:"
          ) ~
            paste(
              "CLASH-supported:",
              "mapped target without",
              "significant response"
            ),
          
          TRUE ~
            "CLASH-supported; target unmapped"
        )
    )
  
} else {
  
  srna_mrna_edges <-
    empty_srna_edge_table()
}


# ============================================================
# OTHER INTERACTION TABLES
# ============================================================


# ------------------------------------------------------------
# 45. mRNA-mRNA edges
# ------------------------------------------------------------

mrna_mrna_source <- clash_mapped |>
  dplyr::filter(
    .data$interaction_class ==
      "mRNA-mRNA"
  )

if (nrow(mrna_mrna_source) > 0) {
  
  mrna_mrna_edges <- mrna_mrna_source |>
    dplyr::transmute(
      
      interaction_row_id =
        .data$interaction_row_id,
      
      source_node =
        dplyr::coalesce(
          .data$endpoint1_mapped_gene_id,
          paste0(
            "RNA:",
            normalize_identifier(
              .data$endpoint1_name
            )
          )
        ),
      
      target_node =
        dplyr::coalesce(
          .data$endpoint2_mapped_gene_id,
          paste0(
            "RNA:",
            normalize_identifier(
              .data$endpoint2_name
            )
          )
        ),
      
      source_label =
        dplyr::coalesce(
          .data$endpoint1_gene_symbol,
          .data$endpoint1_master_label,
          .data$endpoint1_common_name,
          .data$endpoint1_name
        ),
      
      target_label =
        dplyr::coalesce(
          .data$endpoint2_gene_symbol,
          .data$endpoint2_master_label,
          .data$endpoint2_common_name,
          .data$endpoint2_name
        ),
      
      interaction_type =
        "mRNA-mRNA RNase III-CLASH",
      
      total_hybrid_count =
        .data$total_hybrid_count,
      
      raw_p_value =
        .data$raw_p_value,
      
      adjusted_p_value =
        .data$clash_adjusted_p_value,
      
      connection_score =
        .data$connection_score_numeric,
      
      number_of_experiments =
        .data$number_of_experiments,
      
      experimental_dataset =
        "GSE254532"
    )
  
} else {
  
  mrna_mrna_edges <- tibble::tibble(
    
    interaction_row_id = integer(),
    source_node = character(),
    target_node = character(),
    source_label = character(),
    target_label = character(),
    interaction_type = character(),
    total_hybrid_count = numeric(),
    raw_p_value = numeric(),
    adjusted_p_value = numeric(),
    connection_score = numeric(),
    number_of_experiments = numeric(),
    experimental_dataset = character()
  )
}


# ============================================================
# FOCUSED NETWORK TABLES
# ============================================================


# ------------------------------------------------------------
# 46. Vancomycin-responsive CLASH targets
# ------------------------------------------------------------

vancomycin_responsive_targets <-
  srna_mrna_edges |>
  dplyr::filter(
    .data$target_is_significant_deg %in%
      TRUE
  ) |>
  dplyr::arrange(
    .data$target_deseq_adjusted_p,
    dplyr::desc(
      abs(
        .data$target_log2_fold_change
      )
    )
  )


# ------------------------------------------------------------
# 47. Regulatory RNA summary
# ------------------------------------------------------------

if (nrow(srna_mrna_edges) > 0) {
  
  srna_summary <- srna_mrna_edges |>
    dplyr::group_by(
      .data$source_node,
      .data$source_label
    ) |>
    dplyr::summarise(
      
      number_of_mrna_targets =
        dplyr::n_distinct(
          .data$target_node
        ),
      
      number_of_vancomycin_responsive_targets =
        sum(
          .data$target_is_significant_deg,
          na.rm = TRUE
        ),
      
      number_of_strong_response_targets =
        sum(
          .data$target_strong_response,
          na.rm = TRUE
        ),
      
      induced_targets =
        sum(
          .data$target_vancomycin_response ==
            "Induced",
          na.rm = TRUE
        ),
      
      repressed_targets =
        sum(
          .data$target_vancomycin_response ==
            "Repressed",
          na.rm = TRUE
        ),
      
      total_hybrid_count =
        sum_or_na(
          .data$total_hybrid_count
        ),
      
      best_adjusted_p_value =
        minimum_or_na(
          .data$best_adjusted_p_value
        ),
      
      best_connection_score =
        maximum_or_na(
          .data$best_connection_score
        ),
      
      targets =
        collapse_unique(
          .data$target_label
        ),
      
      target_gene_ids =
        collapse_unique(
          .data$target_node
        ),
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      
      network_priority =
        dplyr::case_when(
          
          .data$number_of_strong_response_targets >=
            2 ~
            "Priority 1",
          
          .data$number_of_vancomycin_responsive_targets >=
            2 ~
            "Priority 2",
          
          .data$number_of_vancomycin_responsive_targets ==
            1 ~
            "Priority 3",
          
          TRUE ~
            "Background CLASH RNA"
        ),
      
      priority_order =
        dplyr::case_when(
          
          .data$network_priority ==
            "Priority 1" ~ 1L,
          
          .data$network_priority ==
            "Priority 2" ~ 2L,
          
          .data$network_priority ==
            "Priority 3" ~ 3L,
          
          TRUE ~ 4L
        )
    ) |>
    dplyr::arrange(
      .data$priority_order,
      dplyr::desc(
        .data$
          number_of_vancomycin_responsive_targets
      ),
      dplyr::desc(
        .data$number_of_mrna_targets
      )
    ) |>
    dplyr::select(
      -.data$priority_order
    )
  
} else {
  
  srna_summary <- tibble::tibble(
    
    source_node = character(),
    source_label = character(),
    number_of_mrna_targets = integer(),
    number_of_vancomycin_responsive_targets =
      integer(),
    number_of_strong_response_targets =
      integer(),
    induced_targets = integer(),
    repressed_targets = integer(),
    total_hybrid_count = numeric(),
    best_adjusted_p_value = numeric(),
    best_connection_score = numeric(),
    targets = character(),
    target_gene_ids = character(),
    network_priority = character()
  )
}


# ------------------------------------------------------------
# 48. Build network nodes
# ------------------------------------------------------------

srna_nodes <- srna_mrna_edges |>
  dplyr::transmute(
    
    node_id =
      .data$source_node,
    
    node_label =
      .data$source_label,
    
    node_type =
      "Regulatory RNA",
    
    vancomycin_response =
      "Not directly measured",
    
    log2_fold_change =
      NA_real_,
    
    adjusted_p_value =
      NA_real_,
    
    is_significant_deg =
      FALSE,
    
    strong_response =
      FALSE,
    
    evidence_tier =
      "Experimental CLASH node"
  ) |>
  dplyr::distinct(
    .data$node_id,
    .keep_all = TRUE
  )

target_nodes <- srna_mrna_edges |>
  dplyr::transmute(
    
    node_id =
      .data$target_node,
    
    node_label =
      .data$target_label,
    
    node_type =
      "mRNA target",
    
    vancomycin_response =
      .data$target_vancomycin_response,
    
    log2_fold_change =
      .data$target_log2_fold_change,
    
    adjusted_p_value =
      .data$target_deseq_adjusted_p,
    
    is_significant_deg =
      .data$target_is_significant_deg,
    
    strong_response =
      .data$target_strong_response,
    
    evidence_tier =
      .data$evidence_tier
  ) |>
  dplyr::distinct(
    .data$node_id,
    .keep_all = TRUE
  )

network_nodes <- dplyr::bind_rows(
  srna_nodes,
  target_nodes
) |>
  dplyr::distinct(
    .data$node_id,
    .keep_all = TRUE
  )


# ============================================================
# UPDATE MASTER KNOWLEDGEBASE
# ============================================================


# ------------------------------------------------------------
# 49. Summarize mapped targets
# ------------------------------------------------------------

mapped_target_edges <- srna_mrna_edges |>
  dplyr::filter(
    !stringr::str_detect(
      .data$target_node,
      "^RNA:"
    )
  )

if (nrow(mapped_target_edges) > 0) {
  
  target_summary <- mapped_target_edges |>
    dplyr::group_by(
      .data$target_node
    ) |>
    dplyr::summarise(
      
      clash_target_status_new =
        TRUE,
      
      number_of_clash_srna_regulators_new =
        dplyr::n_distinct(
          .data$source_node
        ),
      
      clash_srna_regulators_new =
        collapse_unique(
          .data$source_label
        ),
      
      clash_total_hybrid_count_new =
        sum_or_na(
          .data$total_hybrid_count
        ),
      
      clash_best_adjusted_p_new =
        minimum_or_na(
          .data$best_adjusted_p_value
        ),
      
      clash_best_connection_score_new =
        maximum_or_na(
          .data$best_connection_score
        ),
      
      clash_evidence_level_new =
        first_character_or_na(
          .data$evidence_tier
        ),
      
      .groups = "drop"
    ) |>
    dplyr::rename(
      gene_id =
        .data$target_node
    )
  
} else {
  
  target_summary <- tibble::tibble(
    
    gene_id = character(),
    clash_target_status_new = logical(),
    number_of_clash_srna_regulators_new =
      integer(),
    clash_srna_regulators_new =
      character(),
    clash_total_hybrid_count_new =
      numeric(),
    clash_best_adjusted_p_new =
      numeric(),
    clash_best_connection_score_new =
      numeric(),
    clash_evidence_level_new =
      character()
  )
}


# ------------------------------------------------------------
# 50. Merge CLASH results into master table
# ------------------------------------------------------------

updated_master_kb <- master_kb |>
  dplyr::left_join(
    target_summary,
    by = "gene_id"
  ) |>
  dplyr::mutate(
    
    clash_target_status =
      dplyr::coalesce(
        .data$clash_target_status_new,
        .data$clash_target_status,
        FALSE
      ),
    
    number_of_clash_srna_regulators =
      dplyr::coalesce(
        .data$
          number_of_clash_srna_regulators_new,
        .data$
          number_of_clash_srna_regulators,
        0L
      ),
    
    clash_srna_regulators =
      dplyr::coalesce(
        .data$clash_srna_regulators_new,
        .data$clash_srna_regulators
      ),
    
    clash_total_hybrid_count =
      dplyr::coalesce(
        .data$clash_total_hybrid_count_new,
        .data$clash_total_hybrid_count
      ),
    
    clash_best_adjusted_p =
      dplyr::coalesce(
        .data$clash_best_adjusted_p_new,
        .data$clash_best_adjusted_p
      ),
    
    clash_best_connection_score =
      dplyr::coalesce(
        .data$
          clash_best_connection_score_new,
        .data$clash_best_connection_score
      ),
    
    clash_evidence_level =
      dplyr::coalesce(
        .data$clash_evidence_level_new,
        .data$clash_evidence_level
      ),
    
    srna_network_role =
      dplyr::case_when(
        
        .data$clash_target_status &
          .data$is_significant_deg ~
          paste(
            "Vancomycin-responsive",
            "CLASH mRNA target"
          ),
        
        .data$clash_target_status ~
          "CLASH mRNA target",
        
        TRUE ~
          .data$srna_network_role
      ),
    
    interaction_evidence =
      dplyr::case_when(
        
        .data$clash_target_status ~
          "GSE254532 RNase III-CLASH",
        
        TRUE ~
          .data$interaction_evidence
      ),
    
    evidence_tier =
      dplyr::case_when(
        
        .data$clash_target_status &
          .data$strong_response ~
          paste(
            "Tier 2-high:",
            "CLASH plus strong",
            "transcriptomic response"
          ),
        
        .data$clash_target_status &
          .data$is_significant_deg ~
          paste(
            "Tier 2:",
            "CLASH plus transcriptomic response"
          ),
        
        .data$clash_target_status ~
          "CLASH-supported interaction",
        
        TRUE ~
          .data$evidence_tier
      )
  ) |>
  dplyr::select(
    -dplyr::ends_with("_new")
  )


# ============================================================
# AUDIT AND SUMMARY
# ============================================================


# ------------------------------------------------------------
# 51. Endpoint mapping audit
# ------------------------------------------------------------

endpoint_mapping_audit <- dplyr::bind_rows(
  
  clash_mapped |>
    dplyr::transmute(
      
      endpoint =
        "Endpoint 1",
      
      raw_name =
        .data$endpoint1_name,
      
      raw_class =
        .data$endpoint1_class_raw,
      
      raw_locus =
        .data$endpoint1_locus_raw,
      
      mapped_gene_id =
        .data$endpoint1_mapped_gene_id,
      
      mapping_method =
        .data$endpoint1_mapping_method,
      
      assigned_type =
        .data$endpoint1_type
    ),
  
  clash_mapped |>
    dplyr::transmute(
      
      endpoint =
        "Endpoint 2",
      
      raw_name =
        .data$endpoint2_name,
      
      raw_class =
        .data$endpoint2_class_raw,
      
      raw_locus =
        .data$endpoint2_locus_raw,
      
      mapped_gene_id =
        .data$endpoint2_mapped_gene_id,
      
      mapping_method =
        .data$endpoint2_mapping_method,
      
      assigned_type =
        .data$endpoint2_type
    )
)


# ------------------------------------------------------------
# 52. Unmapped endpoint summary
#
# No dplyr::count() is used.
# ------------------------------------------------------------

unmapped_endpoint_rows <-
  endpoint_mapping_audit |>
  dplyr::filter(
    is.na(
      .data$mapped_gene_id
    )
  )

if (nrow(unmapped_endpoint_rows) > 0) {
  
  unmapped_endpoints <- unmapped_endpoint_rows |>
    dplyr::group_by(
      .data$raw_name,
      .data$raw_class,
      .data$assigned_type
    ) |>
    dplyr::summarise(
      number_of_occurrences =
        dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(
        .data$number_of_occurrences
      )
    )
  
} else {
  
  unmapped_endpoints <- tibble::tibble(
    
    raw_name = character(),
    raw_class = character(),
    assigned_type = character(),
    number_of_occurrences = integer()
  )
}


# ------------------------------------------------------------
# 53. Interaction-class summary
# ------------------------------------------------------------

interaction_class_summary <- clash_mapped |>
  dplyr::group_by(
    .data$interaction_class
  ) |>
  dplyr::summarise(
    number_of_interactions =
      dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(
      .data$number_of_interactions
    )
  )


# ------------------------------------------------------------
# 54. Main CLASH summary
# ------------------------------------------------------------

clash_summary <- tibble::tibble(
  
  metric = c(
    
    "Raw CLASH interactions",
    
    "Valid standardized interactions",
    
    "Interactions with both endpoints mapped",
    
    "Interactions with at least one endpoint mapped",
    
    "sRNA-mRNA interaction rows",
    
    "Unique collapsed sRNA-mRNA edges",
    
    "Unique regulatory RNA nodes",
    
    "Unique mRNA target nodes",
    
    "Vancomycin-responsive CLASH targets",
    
    "Strong vancomycin-responsive CLASH targets",
    
    "mRNA-mRNA interactions",
    
    "Candidate intergenic RNA-mRNA interactions",
    
    "Other or unclassified interactions"
  ),
  
  value = c(
    
    nrow(raw_clash),
    
    nrow(clash_mapped),
    
    sum(
      clash_mapped$
        both_endpoints_mapped,
      na.rm = TRUE
    ),
    
    sum(
      clash_mapped$
        at_least_one_endpoint_mapped,
      na.rm = TRUE
    ),
    
    nrow(srna_mrna_rows),
    
    nrow(srna_mrna_edges),
    
    dplyr::n_distinct(
      srna_mrna_edges$
        source_node
    ),
    
    dplyr::n_distinct(
      srna_mrna_edges$
        target_node
    ),
    
    sum(
      srna_mrna_edges$
        target_is_significant_deg,
      na.rm = TRUE
    ),
    
    sum(
      srna_mrna_edges$
        target_strong_response,
      na.rm = TRUE
    ),
    
    nrow(mrna_mrna_edges),
    
    sum(
      clash_mapped$
        interaction_class ==
        "candidate_sRNA-mRNA",
      na.rm = TRUE
    ),
    
    sum(
      clash_mapped$
        interaction_class ==
        "Other or unclassified",
      na.rm = TRUE
    )
  )
)


# ============================================================
# SAVE OUTPUTS
# ============================================================


# ------------------------------------------------------------
# 55. Save core tables
# ------------------------------------------------------------

readr::write_csv(
  clash_mapped,
  file.path(
    table_folder,
    "GSE254532_CLASH_standardized_all_interactions.csv"
  )
)

readr::write_csv(
  srna_mrna_rows,
  file.path(
    table_folder,
    "GSE254532_sRNA_mRNA_edges_all_rows.csv"
  )
)

readr::write_csv(
  srna_mrna_edges,
  file.path(
    table_folder,
    "GSE254532_sRNA_mRNA_edges_collapsed.csv"
  )
)

readr::write_csv(
  mrna_mrna_edges,
  file.path(
    table_folder,
    "GSE254532_mRNA_mRNA_edges.csv"
  )
)

readr::write_csv(
  vancomycin_responsive_targets,
  file.path(
    table_folder,
    "GSE254532_vancomycin_responsive_CLASH_targets.csv"
  )
)

readr::write_csv(
  srna_summary,
  file.path(
    table_folder,
    "GSE254532_regulatory_RNA_summary.csv"
  )
)

readr::write_csv(
  network_nodes,
  file.path(
    table_folder,
    "GSE254532_CLASH_network_nodes.csv"
  )
)

readr::write_csv(
  updated_master_kb,
  file.path(
    table_folder,
    "JKD6008_master_gene_knowledgebase_CLASH_updated.csv"
  )
)


# ------------------------------------------------------------
# 56. Save audit files
# ------------------------------------------------------------

readr::write_csv(
  endpoint_mapping_audit,
  file.path(
    audit_folder,
    "GSE254532_endpoint_mapping_audit.csv"
  )
)

readr::write_csv(
  unmapped_endpoints,
  file.path(
    audit_folder,
    "GSE254532_unmapped_endpoints.csv"
  )
)

readr::write_csv(
  interaction_class_summary,
  file.path(
    audit_folder,
    "GSE254532_interaction_class_summary.csv"
  )
)

readr::write_csv(
  clash_summary,
  file.path(
    audit_folder,
    "GSE254532_CLASH_summary.csv"
  )
)


# ------------------------------------------------------------
# 57. Output verification
# ------------------------------------------------------------

expected_outputs <- c(
  
  file.path(
    table_folder,
    "GSE254532_CLASH_standardized_all_interactions.csv"
  ),
  
  file.path(
    table_folder,
    "GSE254532_sRNA_mRNA_edges_all_rows.csv"
  ),
  
  file.path(
    table_folder,
    "GSE254532_sRNA_mRNA_edges_collapsed.csv"
  ),
  
  file.path(
    table_folder,
    "GSE254532_vancomycin_responsive_CLASH_targets.csv"
  ),
  
  file.path(
    table_folder,
    "GSE254532_regulatory_RNA_summary.csv"
  ),
  
  file.path(
    table_folder,
    "GSE254532_CLASH_network_nodes.csv"
  ),
  
  file.path(
    table_folder,
    "JKD6008_master_gene_knowledgebase_CLASH_updated.csv"
  )
)

output_verification <- tibble::tibble(
  
  file =
    expected_outputs,
  
  exists =
    file.exists(
      expected_outputs
    ),
  
  size_KB =
    round(
      file.info(
        expected_outputs
      )$size /
        1024,
      2
    )
)

readr::write_csv(
  output_verification,
  file.path(
    audit_folder,
    "Script13_output_verification.csv"
  )
)


# ------------------------------------------------------------
# 58. Save session information
# ------------------------------------------------------------

writeLines(
  
  capture.output(
    sessionInfo()
  ),
  
  con = file.path(
    audit_folder,
    "sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 59. Final console report
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 13 COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================\n\n"
)

cat(
  "CLASH integration summary:\n"
)

print(
  clash_summary,
  n = Inf
)

cat(
  "\nInteraction classes:\n"
)

print(
  interaction_class_summary,
  n = Inf
)

cat(
  "\nTop regulatory RNA candidates:\n"
)

if (nrow(srna_summary) > 0) {
  
  print(
    
    srna_summary |>
      dplyr::select(
        
        .data$source_label,
        
        .data$number_of_mrna_targets,
        
        .data$
          number_of_vancomycin_responsive_targets,
        
        .data$
          number_of_strong_response_targets,
        
        .data$induced_targets,
        
        .data$repressed_targets,
        
        .data$network_priority
      ) |>
      dplyr::slice_head(
        n = 20
      ),
    
    n = Inf
  )
  
} else {
  
  cat(
    "No sRNA-mRNA interactions were classified.\n"
  )
}

cat(
  "\nOutput verification:\n"
)

print(
  output_verification,
  n = Inf
)

cat(
  "\nMain edge file:\n",
  file.path(
    table_folder,
    "GSE254532_sRNA_mRNA_edges_collapsed.csv"
  ),
  "\n"
)

cat(
  "\nVancomycin-responsive CLASH targets:\n",
  file.path(
    table_folder,
    "GSE254532_vancomycin_responsive_CLASH_targets.csv"
  ),
  "\n"
)

cat(
  "\nRegulatory RNA summary:\n",
  file.path(
    table_folder,
    "GSE254532_regulatory_RNA_summary.csv"
  ),
  "\n"
)

cat(
  "\nUpdated master knowledgebase:\n",
  file.path(
    table_folder,
    "JKD6008_master_gene_knowledgebase_CLASH_updated.csv"
  ),
  "\n"
)

cat(
  "\nNext step:\n",
  "Curate the priority regulatory RNAs and reconstruct ",
  "the vancomycin-focused sRNA-mRNA subnetwork.\n"
)