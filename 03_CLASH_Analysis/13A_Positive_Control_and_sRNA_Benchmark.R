# ============================================================
# BacRegRNA Project
# Script 13A — Positive-Control and sRNA Benchmark
#
# Main objectives:
# 1. Benchmark the CLASH pipeline using known vancomycin-
#    responsive sRNAs.
# 2. Test recovery of RsaOI and SprA2/RsaJ.
# 3. Search for known/expected targets such as hpr and atl.
# 4. Determine exactly why an interaction survives or fails
#    the current vancomycin-responsive target filter.
# 5. Compare positive controls with novel nrdF candidates.
#
# IMPORTANT:
# - No thresholds are modified.
# - No Script 13 output is overwritten.
# - This is an audit/benchmark analysis.
# ============================================================


# ------------------------------------------------------------
# 1. PROJECT PATHS
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

input_file <- file.path(
  project_folder,
  "13_CLASH_network",
  "tables",
  "GSE254532_sRNA_mRNA_edges_collapsed.csv"
)

output_folder <- file.path(
  project_folder,
  "13A_positive_control_benchmark"
)

table_folder <- file.path(
  output_folder,
  "tables"
)

figure_folder <- file.path(
  output_folder,
  "figures"
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
  figure_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "tidyr",
  "tibble",
  "ggplot2"
)

installed_names <- rownames(
  installed.packages()
)

missing_packages <- setdiff(
  required_packages,
  installed_names
)

if (length(missing_packages) > 0) {
  
  install.packages(
    missing_packages
  )
}

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)


# ------------------------------------------------------------
# 3. VERIFY INPUT
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "\nRequired input file was not found:\n",
      input_file,
      "\n\nRun Script 13 first."
    )
  )
}

cat(
  "\n========================================\n",
  "SCRIPT 13A — POSITIVE CONTROL BENCHMARK\n",
  "========================================\n"
)

cat(
  "\nInput file:\n",
  input_file,
  "\n"
)


# ------------------------------------------------------------
# 4. READ COMPLETE CLASH EDGE TABLE
# ------------------------------------------------------------

edges <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)

if (nrow(edges) == 0) {
  
  stop(
    "The collapsed CLASH edge table contains zero rows."
  )
}

cat(
  "\nTotal collapsed CLASH edges:",
  nrow(edges),
  "\n"
)


# ------------------------------------------------------------
# 5. HELPER FUNCTIONS
# ------------------------------------------------------------

safe_character <- function(x) {
  
  x <- as.character(x)
  
  x[is.na(x)] <- ""
  
  x
}


safe_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


safe_logical <- function(x) {
  
  if (is.logical(x)) {
    
    x[is.na(x)] <- FALSE
    
    return(x)
  }
  
  value <- stringr::str_to_lower(
    stringr::str_squish(
      as.character(x)
    )
  )
  
  value %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}


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
# 6. ENSURE REQUIRED COLUMNS EXIST
# ------------------------------------------------------------

required_columns <- list(
  
  source_node = NA_character_,
  target_node = NA_character_,
  
  source_label = NA_character_,
  target_label = NA_character_,
  
  source_raw_names = NA_character_,
  target_raw_names = NA_character_,
  
  source_type = NA_character_,
  target_type = NA_character_,
  
  interaction_class = NA_character_,
  
  number_of_supporting_rows = NA_real_,
  total_hybrid_count = NA_real_,
  
  best_raw_p_value = NA_real_,
  best_adjusted_p_value = NA_real_,
  best_connection_score = NA_real_,
  
  maximum_number_of_experiments = NA_real_,
  
  target_vancomycin_response = NA_character_,
  target_log2_fold_change = NA_real_,
  target_deseq_adjusted_p = NA_real_,
  
  target_is_significant_deg = FALSE,
  target_strong_response = FALSE,
  
  evidence_tier = NA_character_
)

for (column_name in names(required_columns)) {
  
  edges <- add_missing_column(
    edges,
    column_name,
    required_columns[[column_name]]
  )
}


# ------------------------------------------------------------
# 7. STANDARDIZE COLUMN TYPES
# ------------------------------------------------------------

character_columns <- c(
  
  "source_node",
  "target_node",
  
  "source_label",
  "target_label",
  
  "source_raw_names",
  "target_raw_names",
  
  "source_type",
  "target_type",
  
  "interaction_class",
  
  "target_vancomycin_response",
  
  "evidence_tier"
)

for (column_name in character_columns) {
  
  edges[[column_name]] <-
    safe_character(
      edges[[column_name]]
    )
}


numeric_columns <- c(
  
  "number_of_supporting_rows",
  "total_hybrid_count",
  
  "best_raw_p_value",
  "best_adjusted_p_value",
  
  "best_connection_score",
  
  "maximum_number_of_experiments",
  
  "target_log2_fold_change",
  "target_deseq_adjusted_p"
)

for (column_name in numeric_columns) {
  
  edges[[column_name]] <-
    safe_numeric(
      edges[[column_name]]
    )
}


edges$target_is_significant_deg <-
  safe_logical(
    edges$target_is_significant_deg
  )

edges$target_strong_response <-
  safe_logical(
    edges$target_strong_response
  )


# ------------------------------------------------------------
# 8. BUILD SEARCHABLE TEXT
#
# Search across all available identifiers instead of relying
# only on source_label or target_label.
# ------------------------------------------------------------

edges <- edges |>
  dplyr::mutate(
    
    source_search_text =
      stringr::str_to_lower(
        paste(
          .data$source_node,
          .data$source_label,
          .data$source_raw_names,
          .data$source_type,
          .data$interaction_class,
          sep = " | "
        )
      ),
    
    target_search_text =
      stringr::str_to_lower(
        paste(
          .data$target_node,
          .data$target_label,
          .data$target_raw_names,
          .data$target_type,
          sep = " | "
        )
      )
  )


# ------------------------------------------------------------
# 9. IDENTIFY RsaOI
#
# Actual dataset identifier observed in the project:
# Sau-6477/RsaOI
# ------------------------------------------------------------

rsaoi_edges <- edges |>
  dplyr::filter(
    
    stringr::str_detect(
      .data$source_search_text,
      stringr::regex(
        "rsaoi|sau[-_ ]?6477",
        ignore_case = TRUE
      )
    )
  ) |>
  
  dplyr::mutate(
    benchmark_sRNA = "RsaOI"
  )


cat(
  "\n----------------------------------------\n",
  "RsaOI recovery\n",
  "----------------------------------------\n"
)

cat(
  "RsaOI CLASH edges recovered:",
  nrow(rsaoi_edges),
  "\n"
)


if (nrow(rsaoi_edges) > 0) {
  
  print(
    rsaoi_edges |>
      dplyr::select(
        source_label,
        source_raw_names,
        target_node,
        target_label,
        target_raw_names,
        total_hybrid_count,
        maximum_number_of_experiments,
        best_connection_score,
        target_log2_fold_change,
        target_deseq_adjusted_p,
        target_is_significant_deg,
        target_strong_response
      )
  )
}


# ------------------------------------------------------------
# 10. IDENTIFY SprA2 / RsaJ
# ------------------------------------------------------------

spra2_edges <- edges |>
  dplyr::filter(
    
    stringr::str_detect(
      .data$source_search_text,
      stringr::regex(
        "spra2|rsaj",
        ignore_case = TRUE
      )
    )
  ) |>
  
  dplyr::mutate(
    benchmark_sRNA = "SprA2/RsaJ"
  )


cat(
  "\n----------------------------------------\n",
  "SprA2/RsaJ recovery\n",
  "----------------------------------------\n"
)

cat(
  "SprA2/RsaJ CLASH edges recovered:",
  nrow(spra2_edges),
  "\n"
)


if (nrow(spra2_edges) > 0) {
  
  print(
    spra2_edges |>
      dplyr::select(
        source_label,
        source_raw_names,
        target_node,
        target_label,
        target_raw_names,
        total_hybrid_count,
        maximum_number_of_experiments,
        best_connection_score,
        target_log2_fold_change,
        target_deseq_adjusted_p,
        target_is_significant_deg
      )
  )
}


# ------------------------------------------------------------
# 11. SEARCH FOR KNOWN TARGET NAMES
#
# Diagnostic only.
# Do NOT interpret absence as proof that the interaction
# does not exist; identifier mapping may differ.
# ------------------------------------------------------------

known_target_patterns <- c(
  "hpr",
  "atl"
)

known_target_regex <- paste(
  known_target_patterns,
  collapse = "|"
)


known_target_edges <- edges |>
  dplyr::filter(
    
    stringr::str_detect(
      .data$target_search_text,
      stringr::regex(
        known_target_regex,
        ignore_case = TRUE
      )
    )
  )


cat(
  "\n----------------------------------------\n",
  "Known target-name search\n",
  "----------------------------------------\n"
)

cat(
  "Edges matching hpr/atl identifiers:",
  nrow(known_target_edges),
  "\n"
)


if (nrow(known_target_edges) > 0) {
  
  print(
    known_target_edges |>
      dplyr::select(
        source_label,
        source_raw_names,
        target_node,
        target_label,
        target_raw_names,
        total_hybrid_count,
        maximum_number_of_experiments,
        target_log2_fold_change,
        target_is_significant_deg
      )
  )
}


# ------------------------------------------------------------
# 12. SEARCH RsaOI × hpr/atl DIRECTLY
# ------------------------------------------------------------

rsaoi_known_targets <- rsaoi_edges |>
  dplyr::filter(
    
    stringr::str_detect(
      .data$target_search_text,
      stringr::regex(
        known_target_regex,
        ignore_case = TRUE
      )
    )
  )


cat(
  "\nRsaOI -> hpr/atl edges recovered by gene-name search:",
  nrow(rsaoi_known_targets),
  "\n"
)


# ------------------------------------------------------------
# 13. IDENTIFY nrdF CANDIDATES
# ------------------------------------------------------------

nrdf_edges <- edges |>
  dplyr::filter(
    
    stringr::str_detect(
      .data$target_search_text,
      stringr::regex(
        "\\bnrdf\\b",
        ignore_case = TRUE
      )
    )
  ) |>
  
  dplyr::mutate(
    benchmark_sRNA = "nrdF candidate"
  )


cat(
  "\n----------------------------------------\n",
  "nrdF candidate recovery\n",
  "----------------------------------------\n"
)

cat(
  "nrdF CLASH edges recovered:",
  nrow(nrdf_edges),
  "\n"
)


if (nrow(nrdf_edges) > 0) {
  
  print(
    nrdf_edges |>
      dplyr::select(
        source_label,
        source_raw_names,
        target_node,
        target_label,
        total_hybrid_count,
        maximum_number_of_experiments,
        best_connection_score,
        target_log2_fold_change,
        target_deseq_adjusted_p,
        target_is_significant_deg,
        target_strong_response
      )
  )
}


# ------------------------------------------------------------
# 14. DEFINE CURRENT SCRIPT-14 SURVIVAL LOGIC
#
# Current focused-network logic:
# target must be a significant DEG.
#
# This DOES NOT modify that rule.
# It only audits its consequence.
# ------------------------------------------------------------

benchmark_edges <- dplyr::bind_rows(
  
  rsaoi_edges,
  
  spra2_edges,
  
  nrdf_edges
  
) |>
  
  dplyr::distinct(
    .data$source_node,
    .data$target_node,
    .keep_all = TRUE
  ) |>
  
  dplyr::mutate(
    
    survives_target_DE_filter =
      .data$target_is_significant_deg,
    
    exclusion_reason =
      dplyr::case_when(
        
        .data$target_is_significant_deg ~
          "Retained: target is a significant DEG",
        
        is.na(
          .data$target_log2_fold_change
        ) ~
          paste0(
            "Excluded from DEG-focused network: ",
            "target has no mapped DE value"
          ),
        
        !is.na(
          .data$target_log2_fold_change
        ) &
          !.data$target_is_significant_deg ~
          paste0(
            "Excluded from DEG-focused network: ",
            "target does not pass DEG threshold"
          ),
        
        TRUE ~
          "Excluded from DEG-focused network"
      )
  )


# ------------------------------------------------------------
# 15. ASSIGN BENCHMARK CATEGORY
# ------------------------------------------------------------

benchmark_edges <- benchmark_edges |>
  dplyr::mutate(
    
    benchmark_category =
      dplyr::case_when(
        
        stringr::str_detect(
          .data$source_search_text,
          stringr::regex(
            "rsaoi|sau[-_ ]?6477",
            ignore_case = TRUE
          )
        ) ~
          "Positive-control sRNA: RsaOI",
        
        stringr::str_detect(
          .data$source_search_text,
          stringr::regex(
            "spra2|rsaj",
            ignore_case = TRUE
          )
        ) ~
          "Known sRNA: SprA2/RsaJ",
        
        stringr::str_detect(
          .data$target_search_text,
          stringr::regex(
            "\\bnrdf\\b",
            ignore_case = TRUE
          )
        ) ~
          "Novel nrdF candidate",
        
        TRUE ~
          "Other"
      )
  )


# ------------------------------------------------------------
# 16. CREATE COMPACT BENCHMARK TABLE
# ------------------------------------------------------------

benchmark_table <- benchmark_edges |>
  dplyr::transmute(
    
    benchmark_category,
    
    source_node,
    source_label,
    source_raw_names,
    
    target_node,
    target_label,
    target_raw_names,
    
    interaction_class,
    
    supporting_rows =
      number_of_supporting_rows,
    
    hybrid_count =
      total_hybrid_count,
    
    experiments =
      maximum_number_of_experiments,
    
    connection_score =
      best_connection_score,
    
    clash_p =
      best_raw_p_value,
    
    clash_adjusted_p =
      best_adjusted_p_value,
    
    target_log2FC =
      target_log2_fold_change,
    
    target_DE_padj =
      target_deseq_adjusted_p,
    
    target_significant_DEG =
      target_is_significant_deg,
    
    target_strong_response =
      target_strong_response,
    
    survives_target_DE_filter,
    
    exclusion_reason
  ) |>
  
  dplyr::arrange(
    benchmark_category,
    dplyr::desc(hybrid_count)
  )


# ------------------------------------------------------------
# 17. SAVE BENCHMARK TABLE
# ------------------------------------------------------------

readr::write_csv(
  
  benchmark_table,
  
  file.path(
    table_folder,
    "positive_control_benchmark_all.csv"
  )
)


# ------------------------------------------------------------
# 18. SAVE RsaOI-SPECIFIC TABLE
# ------------------------------------------------------------

readr::write_csv(
  
  rsaoi_edges,
  
  file.path(
    table_folder,
    "RsaOI_all_CLASH_edges.csv"
  )
)


# ------------------------------------------------------------
# 19. SAVE SprA2/RsaJ TABLE
# ------------------------------------------------------------

readr::write_csv(
  
  spra2_edges,
  
  file.path(
    table_folder,
    "SprA2_RsaJ_all_CLASH_edges.csv"
  )
)


# ------------------------------------------------------------
# 20. SAVE nrdF TABLE
# ------------------------------------------------------------

readr::write_csv(
  
  nrdf_edges,
  
  file.path(
    table_folder,
    "nrdF_candidate_CLASH_edges.csv"
  )
)


# ------------------------------------------------------------
# 21. SAVE KNOWN-TARGET SEARCH
# ------------------------------------------------------------

readr::write_csv(
  
  known_target_edges,
  
  file.path(
    table_folder,
    "known_target_hpr_atl_search.csv"
  )
)


readr::write_csv(
  
  rsaoi_known_targets,
  
  file.path(
    table_folder,
    "RsaOI_hpr_atl_direct_search.csv"
  )
)


# ------------------------------------------------------------
# 22. CREATE DIAGNOSTIC SUMMARY
# ------------------------------------------------------------

diagnostic_summary <- tibble::tibble(
  
  diagnostic = c(
    
    "Total collapsed CLASH edges",
    
    "RsaOI CLASH edges recovered",
    
    "SprA2/RsaJ CLASH edges recovered",
    
    "hpr/atl target-name matches",
    
    "RsaOI-hpr/atl direct matches",
    
    "nrdF CLASH edges",
    
    "RsaOI edges surviving current DEG filter",
    
    "SprA2/RsaJ edges surviving current DEG filter",
    
    "nrdF edges surviving current DEG filter"
  ),
  
  value = c(
    
    nrow(edges),
    
    nrow(rsaoi_edges),
    
    nrow(spra2_edges),
    
    nrow(known_target_edges),
    
    nrow(rsaoi_known_targets),
    
    nrow(nrdf_edges),
    
    sum(
      rsaoi_edges$target_is_significant_deg,
      na.rm = TRUE
    ),
    
    sum(
      spra2_edges$target_is_significant_deg,
      na.rm = TRUE
    ),
    
    sum(
      nrdf_edges$target_is_significant_deg,
      na.rm = TRUE
    )
  )
)


readr::write_csv(
  
  diagnostic_summary,
  
  file.path(
    table_folder,
    "positive_control_diagnostic_summary.csv"
  )
)


cat(
  "\n----------------------------------------\n",
  "DIAGNOSTIC SUMMARY\n",
  "----------------------------------------\n"
)

print(
  diagnostic_summary
)


# ------------------------------------------------------------
# 23. FIGURE — HYBRID SUPPORT
# ------------------------------------------------------------

plot_data <- benchmark_table |>
  
  dplyr::filter(
    !is.na(hybrid_count)
  ) |>
  
  dplyr::mutate(
    
    interaction_label =
      paste0(
        ifelse(
          source_label == "",
          source_raw_names,
          source_label
        ),
        " -> ",
        ifelse(
          target_label == "",
          target_node,
          target_label
        )
      )
  )


if (nrow(plot_data) > 0) {
  
  p1 <- ggplot(
    plot_data,
    aes(
      x = reorder(
        interaction_label,
        hybrid_count
      ),
      y = hybrid_count
    )
  ) +
    
    geom_col() +
    
    coord_flip() +
    
    facet_wrap(
      ~ benchmark_category,
      scales = "free_y"
    ) +
    
    labs(
      title =
        "Positive-control benchmark of CLASH-supported interactions",
      subtitle =
        "Known sRNAs and candidate nrdF interactions",
      x =
        "RNA-mRNA interaction",
      y =
        "Total CLASH hybrid count"
    ) +
    
    theme_bw(base_size = 12)
  
  
  ggsave(
    
    filename = file.path(
      figure_folder,
      "Positive_Control_CLASH_Hybrid_Counts.png"
    ),
    
    plot = p1,
    
    width = 11,
    height = 7,
    dpi = 400
  )
  
  
  ggsave(
    
    filename = file.path(
      figure_folder,
      "Positive_Control_CLASH_Hybrid_Counts.pdf"
    ),
    
    plot = p1,
    
    width = 11,
    height = 7
  )
}


# ------------------------------------------------------------
# 24. FIGURE — CURRENT DEG FILTER SURVIVAL
# ------------------------------------------------------------

survival_summary <- benchmark_table |>
  
  dplyr::count(
    benchmark_category,
    survives_target_DE_filter,
    name = "number_of_edges"
  ) |>
  
  dplyr::mutate(
    
    filter_status =
      ifelse(
        survives_target_DE_filter,
        "Retained",
        "Excluded"
      )
  )


if (nrow(survival_summary) > 0) {
  
  p2 <- ggplot(
    survival_summary,
    aes(
      x = benchmark_category,
      y = number_of_edges,
      fill = filter_status
    )
  ) +
    
    geom_col(
      position = "stack"
    ) +
    
    coord_flip() +
    
    labs(
      title =
        "Effect of the current target-DE filter on benchmark interactions",
      x =
        NULL,
      y =
        "Number of CLASH interactions",
      fill =
        "Current filter"
    ) +
    
    theme_bw(base_size = 12)
  
  
  ggsave(
    
    filename = file.path(
      figure_folder,
      "Benchmark_Current_DE_Filter_Survival.png"
    ),
    
    plot = p2,
    
    width = 10,
    height = 6,
    dpi = 400
  )
}


# ------------------------------------------------------------
# 25. WRITE HUMAN-READABLE AUDIT REPORT
# ------------------------------------------------------------

audit_file <- file.path(
  output_folder,
  "Positive_Control_Audit_Report.txt"
)


sink(audit_file)

cat(
  "BacRegRNA Project\n"
)

cat(
  "Positive-Control Benchmark Audit\n"
)

cat(
  "====================================\n\n"
)

cat(
  "Total collapsed CLASH edges: ",
  nrow(edges),
  "\n\n",
  sep = ""
)


cat(
  "RsaOI\n",
  "-----\n",
  "CLASH edges recovered: ",
  nrow(rsaoi_edges),
  "\n",
  sep = ""
)

cat(
  "Edges surviving current target-DE filter: ",
  sum(
    rsaoi_edges$target_is_significant_deg,
    na.rm = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "SprA2/RsaJ\n",
  "----------\n",
  "CLASH edges recovered: ",
  nrow(spra2_edges),
  "\n",
  sep = ""
)

cat(
  "Edges surviving current target-DE filter: ",
  sum(
    spra2_edges$target_is_significant_deg,
    na.rm = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Known target-name benchmark\n",
  "---------------------------\n",
  "hpr/atl matches across all edges: ",
  nrow(known_target_edges),
  "\n",
  sep = ""
)

cat(
  "RsaOI-hpr/atl direct matches: ",
  nrow(rsaoi_known_targets),
  "\n\n",
  sep = ""
)


cat(
  "nrdF candidates\n",
  "---------------\n",
  "CLASH edges recovered: ",
  nrow(nrdf_edges),
  "\n",
  sep = ""
)

cat(
  "Edges surviving current target-DE filter: ",
  sum(
    nrdf_edges$target_is_significant_deg,
    na.rm = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Interpretation rule\n",
  "-------------------\n"
)

cat(
  paste0(
    "This script is diagnostic only. ",
    "Failure of a known interaction to enter the ",
    "vancomycin-focused network does not demonstrate ",
    "absence of post-transcriptional regulation. ",
    "It identifies the pipeline step responsible for ",
    "exclusion, including missing target-DE mapping or ",
    "failure of the target to satisfy the current DEG criterion.\n"
  )
)

sink()


# ------------------------------------------------------------
# 26. SESSION INFORMATION
# ------------------------------------------------------------

capture.output(
  
  sessionInfo(),
  
  file = file.path(
    output_folder,
    "sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 27. FINAL CONSOLE MESSAGE
# ------------------------------------------------------------

cat(
  "\n========================================\n"
)

cat(
  "SCRIPT 13A COMPLETED SUCCESSFULLY\n"
)

cat(
  "========================================\n"
)

cat(
  "\nResults saved in:\n",
  output_folder,
  "\n"
)

cat(
  "\nKey files to inspect:\n"
)

cat(
  "1. tables/positive_control_diagnostic_summary.csv\n"
)

cat(
  "2. tables/positive_control_benchmark_all.csv\n"
)

cat(
  "3. tables/RsaOI_all_CLASH_edges.csv\n"
)

cat(
  "4. tables/SprA2_RsaJ_all_CLASH_edges.csv\n"
)

cat(
  "5. tables/RsaOI_hpr_atl_direct_search.csv\n"
)

cat(
  "6. tables/nrdF_candidate_CLASH_edges.csv\n"
)

cat(
  "7. Positive_Control_Audit_Report.txt\n"
)

cat(
  "\nIMPORTANT:\n",
  "Do not modify Script 14 yet.\n",
  "Review the benchmark outputs first.\n"
)