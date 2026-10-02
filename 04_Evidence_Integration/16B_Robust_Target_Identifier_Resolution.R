# ============================================================
# SCRIPT 16B
# Robust Target Identifier Resolution
# for the Evidence-Layered CLASH Network
#
# Purpose:
# - Resolve CLASH target identifiers against the JKD6008 master KB
# - Search across multiple identifier fields
# - Audit matched and unmatched targets
# - Avoid changing CLASH or transcriptomic evidence
# ============================================================


# ============================================================
# 01. PROJECT PATHS
# ============================================================

project_folder <- "D:/Bac-sRNA"

network_file <- file.path(
  project_folder,
  "14_revised_evidence_layered_network",
  "tables",
  "Evidence_Layered_CLASH_Network_ALL.csv"
)

master_candidates <- c(
  
  file.path(
    project_folder,
    "13_CLASH_network",
    "tables",
    "JKD6008_master_gene_knowledgebase_CLASH_updated.csv"
  ),
  
  file.path(
    project_folder,
    "12_master_knowledgebase",
    "tables",
    "JKD6008_master_gene_knowledgebase.csv"
  )
)

existing_master_files <- master_candidates[
  file.exists(master_candidates)
]

if (length(existing_master_files) == 0) {
  
  stop(
    "No master knowledgebase found."
  )
}

master_file <- existing_master_files[1]


# ============================================================
# 02. OUTPUT FOLDER
# ============================================================

output_folder <- file.path(
  project_folder,
  "16B_identifier_resolution"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 03. PACKAGES
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "tibble"
)

missing_packages <- setdiff(
  required_packages,
  rownames(installed.packages())
)

if (length(missing_packages) > 0) {
  
  install.packages(
    missing_packages,
    dependencies = TRUE
  )
}

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)


# ============================================================
# 04. HELPER FUNCTIONS
# ============================================================

safe_chr <- function(x) {
  
  x <- as.character(x)
  
  x[is.na(x)] <- ""
  
  stringr::str_squish(x)
}


normalize_identifier <- function(x) {
  
  x <- safe_chr(x)
  
  x <- stringr::str_to_lower(x)
  
  x <- stringr::str_replace_all(
    x,
    "\\s+",
    ""
  )
  
  x <- stringr::str_replace_all(
    x,
    "^gene:",
    ""
  )
  
  x
}


first_nonblank <- function(...) {
  
  xs <- list(...)
  
  n <- max(
    vapply(
      xs,
      length,
      integer(1)
    )
  )
  
  out <- rep(
    NA_character_,
    n
  )
  
  for (x in xs) {
    
    x <- rep_len(
      as.character(x),
      n
    )
    
    take <-
      (
        is.na(out) |
          out == ""
      ) &
      !is.na(x) &
      stringr::str_squish(x) != ""
    
    out[take] <- x[take]
  }
  
  out
}


# ============================================================
# 05. READ FILES
# ============================================================

if (!file.exists(network_file)) {
  
  stop(
    paste(
      "Network file not found:",
      network_file
    )
  )
}

network <- readr::read_csv(
  network_file,
  show_col_types = FALSE
)

master <- readr::read_csv(
  master_file,
  show_col_types = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 16B - IDENTIFIER RESOLUTION\n"
)

cat(
  "============================================\n\n"
)

cat(
  "Network interactions:",
  nrow(network),
  "\n"
)

cat(
  "Master KB rows:",
  nrow(master),
  "\n"
)

cat(
  "Master KB file:\n",
  master_file,
  "\n\n"
)


# ============================================================
# 06. REPORT AVAILABLE MASTER COLUMNS
# ============================================================

cat(
  "MASTER KB COLUMNS\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  names(master)
)


# ============================================================
# 07. ENSURE CORE MASTER COLUMNS EXIST
# ============================================================

required_master_columns <- c(
  "gene_id",
  "current_locus_tag",
  "gene_symbol",
  "gene_synonym",
  "master_label",
  "product",
  "product_short"
)


for (nm in required_master_columns) {
  
  if (!nm %in% names(master)) {
    
    master[[nm]] <- NA_character_
  }
}


# ============================================================
# 08. IDENTIFY ADDITIONAL ID COLUMNS
# ============================================================

candidate_id_patterns <- c(
  "refseq",
  "protein",
  "accession",
  "locus",
  "gene_id",
  "gene_symbol",
  "synonym",
  "master_label"
)


additional_id_columns <- names(master)[
  
  stringr::str_detect(
    
    stringr::str_to_lower(
      names(master)
    ),
    
    paste(
      candidate_id_patterns,
      collapse = "|"
    )
  )
]


additional_id_columns <- unique(
  c(
    required_master_columns[
      required_master_columns %in%
        names(master)
    ],
    additional_id_columns
  )
)


cat(
  "\nCandidate identifier columns:\n"
)

print(
  additional_id_columns
)


# ============================================================
# 09. CREATE MASTER ROW IDENTIFIER
# ============================================================

master <- master |>
  
  dplyr::mutate(
    
    master_row_id =
      dplyr::row_number(),
    
    gene_id =
      safe_chr(
        .data$gene_id
      ),
    
    current_locus_tag =
      safe_chr(
        .data$current_locus_tag
      ),
    
    gene_symbol =
      safe_chr(
        .data$gene_symbol
      ),
    
    gene_synonym =
      safe_chr(
        .data$gene_synonym
      ),
    
    master_label =
      safe_chr(
        .data$master_label
      ),
    
    product =
      safe_chr(
        .data$product
      ),
    
    product_short =
      safe_chr(
        .data$product_short
      )
  )


# ============================================================
# 10. BUILD LONG MASTER IDENTIFIER LOOKUP
# ============================================================

master_lookup <- master |>
  
  dplyr::select(
    
    "master_row_id",
    
    dplyr::all_of(
      additional_id_columns
    )
  ) |>
  
  tidyr::pivot_longer(
    
    cols =
      -master_row_id,
    
    names_to =
      "identifier_type",
    
    values_to =
      "identifier_value"
  ) |>
  
  dplyr::mutate(
    
    identifier_value =
      safe_chr(
        .data$identifier_value
      ),
    
    identifier_normalized =
      normalize_identifier(
        .data$identifier_value
      )
  ) |>
  
  dplyr::filter(
    
    .data$identifier_normalized != ""
  ) |>
  
  dplyr::distinct(
    
    .data$master_row_id,
    
    .data$identifier_type,
    
    .data$identifier_normalized,
    
    .keep_all = TRUE
  )


# ============================================================
# 11. IDENTIFIER TYPE PRIORITY
# ============================================================

master_lookup <- master_lookup |>
  
  dplyr::mutate(
    
    identifier_priority =
      dplyr::case_when(
        
        .data$identifier_type ==
          "gene_id" ~ 1L,
        
        stringr::str_detect(
          .data$identifier_type,
          stringr::regex(
            "refseq|protein.*accession",
            ignore_case = TRUE
          )
        ) ~ 2L,
        
        .data$identifier_type ==
          "current_locus_tag" ~ 3L,
        
        .data$identifier_type ==
          "gene_symbol" ~ 4L,
        
        .data$identifier_type ==
          "gene_synonym" ~ 5L,
        
        .data$identifier_type ==
          "master_label" ~ 6L,
        
        TRUE ~ 10L
      )
  )


# ============================================================
# 12. NETWORK TARGET IDENTIFIER COLUMNS
# ============================================================

candidate_network_columns <- c(
  "target_node",
  "target_label",
  "target_raw_names",
  "target_gene_id",
  "target_current_locus_tag",
  "target_gene_symbol",
  "target_refseq",
  "target_accession"
)


network_id_columns <- intersect(
  candidate_network_columns,
  names(network)
)


if (length(network_id_columns) == 0) {
  
  stop(
    "No usable target identifier columns found in network."
  )
}


cat(
  "\nNetwork target identifier columns:\n"
)

print(
  network_id_columns
)


# ============================================================
# 13. CREATE NETWORK EDGE ID
# ============================================================

network <- network |>
  
  dplyr::mutate(
    
    resolution_edge_id =
      dplyr::row_number()
  )


# ============================================================
# 14. BUILD LONG NETWORK IDENTIFIER TABLE
# ============================================================

network_lookup <- network |>
  
  dplyr::select(
    
    "resolution_edge_id",
    
    dplyr::all_of(
      network_id_columns
    )
  ) |>
  
  tidyr::pivot_longer(
    
    cols =
      -resolution_edge_id,
    
    names_to =
      "network_identifier_type",
    
    values_to =
      "network_identifier_value"
  ) |>
  
  dplyr::mutate(
    
    network_identifier_value =
      safe_chr(
        .data$network_identifier_value
      )
  )


# ============================================================
# 15. SPLIT MULTI-VALUE IDENTIFIERS
#
# target_raw_names may contain multiple aliases.
# ============================================================

network_lookup <- network_lookup |>
  
  tidyr::separate_rows(
    
    network_identifier_value,
    
    sep = "[;,|]"
  ) |>
  
  dplyr::mutate(
    
    network_identifier_value =
      safe_chr(
        .data$network_identifier_value
      ),
    
    identifier_normalized =
      normalize_identifier(
        .data$network_identifier_value
      )
  ) |>
  
  dplyr::filter(
    
    .data$identifier_normalized != ""
  ) |>
  
  dplyr::distinct()


# ============================================================
# 16. NETWORK IDENTIFIER PRIORITY
# ============================================================

network_lookup <- network_lookup |>
  
  dplyr::mutate(
    
    network_identifier_priority =
      dplyr::case_when(
        
        .data$network_identifier_type ==
          "target_node" ~ 1L,
        
        .data$network_identifier_type ==
          "target_gene_id" ~ 2L,
        
        .data$network_identifier_type ==
          "target_refseq" ~ 3L,
        
        .data$network_identifier_type ==
          "target_current_locus_tag" ~ 4L,
        
        .data$network_identifier_type ==
          "target_gene_symbol" ~ 5L,
        
        .data$network_identifier_type ==
          "target_label" ~ 6L,
        
        .data$network_identifier_type ==
          "target_raw_names" ~ 7L,
        
        TRUE ~ 10L
      )
  )


# ============================================================
# 17. MATCH NETWORK IDENTIFIERS TO MASTER LOOKUP
# ============================================================

candidate_matches <- network_lookup |>
  
  dplyr::inner_join(
    
    master_lookup,
    
    by =
      "identifier_normalized"
  )


# ============================================================
# 18. CALCULATE MATCH PRIORITY
# ============================================================

candidate_matches <- candidate_matches |>
  
  dplyr::mutate(
    
    match_priority =
      .data$network_identifier_priority *
      100L +
      .data$identifier_priority
  )


# ============================================================
# 19. COUNT UNIQUE MASTER ROWS PER EDGE
# ============================================================

match_ambiguity <- candidate_matches |>
  
  dplyr::group_by(
    .data$resolution_edge_id
  ) |>
  
  dplyr::summarise(
    
    number_of_candidate_master_rows =
      dplyr::n_distinct(
        .data$master_row_id
      ),
    
    .groups =
      "drop"
  )


# ============================================================
# 20. SELECT BEST MATCH
# ============================================================

best_match <- candidate_matches |>
  
  dplyr::arrange(
    
    .data$resolution_edge_id,
    
    .data$match_priority,
    
    .data$master_row_id
  ) |>
  
  dplyr::group_by(
    .data$resolution_edge_id
  ) |>
  
  dplyr::slice_head(
    n = 1
  ) |>
  
  dplyr::ungroup() |>
  
  dplyr::left_join(
    
    match_ambiguity,
    
    by =
      "resolution_edge_id"
  )


# ============================================================
# 21. ATTACH MASTER ANNOTATION
# ============================================================

master_annotation_columns <- master |>
  
  dplyr::select(
    
    "master_row_id",
    
    dplyr::everything()
  )


resolved_network <- network |>
  
  dplyr::left_join(
    
    best_match |>
      
      dplyr::select(
        
        "resolution_edge_id",
        
        "master_row_id",
        
        "network_identifier_type",
        
        "network_identifier_value",
        
        "identifier_type",
        
        "identifier_value",
        
        "identifier_normalized",
        
        "number_of_candidate_master_rows"
      ),
    
    by =
      "resolution_edge_id"
  ) |>
  
  dplyr::left_join(
    
    master_annotation_columns,
    
    by =
      "master_row_id",
    
    suffix =
      c(
        "",
        "_master"
      )
  )


# ============================================================
# 22. MATCH STATUS
# ============================================================

resolved_network <- resolved_network |>
  
  dplyr::mutate(
    
    identifier_match_status =
      dplyr::case_when(
        
        is.na(
          .data$master_row_id
        ) ~
          
          "Unmatched",
        
        .data$number_of_candidate_master_rows > 1 ~
          
          "Matched - ambiguous candidates",
        
        TRUE ~
          
          "Matched - unique"
      ),
    
    identifier_match_method =
      dplyr::case_when(
        
        is.na(
          .data$master_row_id
        ) ~
          
          "No master KB match",
        
        TRUE ~
          
          paste0(
            
            .data$network_identifier_type,
            
            " -> ",
            
            .data$identifier_type
          )
      )
  )


# ============================================================
# 23. RESOLVED GENE LABEL
# ============================================================

resolved_network <- resolved_network |>
  
  dplyr::mutate(
    
    resolved_gene_symbol =
      first_nonblank(
        
        .data$gene_symbol,
        
        if (
          "gene_symbol_master" %in%
          names(resolved_network)
        ) {
          .data$gene_symbol_master
        } else {
          rep(
            NA_character_,
            dplyr::n()
          )
        },
        
        .data$target_label,
        
        .data$target_node
      ),
    
    resolved_product =
      first_nonblank(
        
        .data$product_short,
        
        .data$product,
        
        .data$target_label,
        
        rep(
          "Uncharacterized product",
          dplyr::n()
        )
      )
  )


# ============================================================
# 24. MATCH SUMMARY
# ============================================================

match_summary <- resolved_network |>
  
  dplyr::count(
    
    .data$identifier_match_status,
    
    name =
      "number_of_interactions"
  ) |>
  
  dplyr::mutate(
    
    percent_of_network =
      round(
        
        100 *
          .data$number_of_interactions /
          nrow(resolved_network),
        
        2
      )
  )


# ============================================================
# 25. UNIQUE TARGET MATCH SUMMARY
# ============================================================

target_match_summary <- resolved_network |>
  
  dplyr::group_by(
    .data$target_node
  ) |>
  
  dplyr::summarise(
    
    target_label =
      dplyr::first(
        .data$target_label
      ),
    
    matched =
      any(
        !is.na(
          .data$master_row_id
        )
      ),
    
    master_row_id =
      dplyr::first(
        .data$master_row_id[
          !is.na(
            .data$master_row_id
          )
        ]
      ),
    
    resolved_gene_symbol =
      dplyr::first(
        .data$resolved_gene_symbol
      ),
    
    resolved_product =
      dplyr::first(
        .data$resolved_product
      ),
    
    .groups =
      "drop"
  )


# ============================================================
# 26. OVERALL RESOLUTION STATISTICS
# ============================================================

total_interactions <-
  nrow(
    resolved_network
  )


matched_interactions <-
  sum(
    !is.na(
      resolved_network$
        master_row_id
    )
  )


unique_targets <-
  dplyr::n_distinct(
    resolved_network$
      target_node
  )


matched_targets <-
  sum(
    target_match_summary$
      matched
  )


interaction_match_rate <-
  round(
    100 *
      matched_interactions /
      total_interactions,
    2
  )


target_match_rate <-
  round(
    100 *
      matched_targets /
      unique_targets,
    2
  )


resolution_summary <- tibble::tibble(
  
  metric =
    c(
      
      "Total network interactions",
      
      "Matched interactions",
      
      "Unmatched interactions",
      
      "Interaction match rate (%)",
      
      "Unique network targets",
      
      "Matched unique targets",
      
      "Unmatched unique targets",
      
      "Unique target match rate (%)"
    ),
  
  value =
    c(
      
      total_interactions,
      
      matched_interactions,
      
      total_interactions -
        matched_interactions,
      
      interaction_match_rate,
      
      unique_targets,
      
      matched_targets,
      
      unique_targets -
        matched_targets,
      
      target_match_rate
    )
)


# ============================================================
# 27. UNMATCHED TARGET AUDIT
# ============================================================

unmatched_targets <- target_match_summary |>
  
  dplyr::filter(
    !.data$matched
  ) |>
  
  dplyr::arrange(
    .data$target_node
  )


# ============================================================
# 28. AMBIGUOUS MATCH AUDIT
# ============================================================

ambiguous_matches <- resolved_network |>
  
  dplyr::filter(
    
    .data$identifier_match_status ==
      "Matched - ambiguous candidates"
  ) |>
  
  dplyr::select(
    
    "resolution_edge_id",
    
    "source_label",
    
    "target_node",
    
    "target_label",
    
    "network_identifier_type",
    
    "network_identifier_value",
    
    "identifier_type",
    
    "identifier_value",
    
    "number_of_candidate_master_rows",
    
    "resolved_gene_symbol",
    
    "resolved_product"
  )


# ============================================================
# 29. MATCH METHOD SUMMARY
# ============================================================

match_method_summary <- resolved_network |>
  
  dplyr::filter(
    
    !is.na(
      .data$master_row_id
    )
  ) |>
  
  dplyr::count(
    
    .data$identifier_match_method,
    
    name =
      "number_of_interactions"
  ) |>
  
  dplyr::arrange(
    
    dplyr::desc(
      .data$number_of_interactions
    )
  )


# ============================================================
# 30. CRITICAL INTERACTION AUDIT
# ============================================================

critical_resolution <- resolved_network |>
  
  dplyr::filter(
    
    stringr::str_detect(
      
      stringr::str_to_lower(
        
        paste(
          
          .data$source_label,
          
          .data$source_raw_names
        )
      ),
      
      "rsaoi|sau[-_ ]?6477|spra2|rsaj"
    ) |
      
      stringr::str_detect(
        
        stringr::str_to_lower(
          
          paste(
            
            .data$target_label,
            
            .data$target_node,
            
            .data$resolved_gene_symbol
          )
        ),
        
        "\\bnrdf\\b"
      )
  ) |>
  
  dplyr::select(
    
    "source_label",
    
    "target_node",
    
    "target_label",
    
    "identifier_match_status",
    
    "identifier_match_method",
    
    "resolved_gene_symbol",
    
    "current_locus_tag",
    
    "gene_id",
    
    "resolved_product",
    
    "number_of_candidate_master_rows"
  )


# ============================================================
# 31. SAVE RESOLVED NETWORK
# ============================================================

readr::write_csv(
  
  resolved_network,
  
  file.path(
    
    output_folder,
    
    "Evidence_Layered_Network_ID_Resolved.csv"
  )
)


# ============================================================
# 32. SAVE TARGET RESOLUTION TABLE
# ============================================================

readr::write_csv(
  
  target_match_summary,
  
  file.path(
    
    output_folder,
    
    "Unique_Target_Identifier_Resolution.csv"
  )
)


# ============================================================
# 33. SAVE UNMATCHED TARGETS
# ============================================================

readr::write_csv(
  
  unmatched_targets,
  
  file.path(
    
    output_folder,
    
    "Unmatched_Targets.csv"
  )
)


# ============================================================
# 34. SAVE AMBIGUOUS MATCHES
# ============================================================

readr::write_csv(
  
  ambiguous_matches,
  
  file.path(
    
    output_folder,
    
    "Ambiguous_Target_Matches.csv"
  )
)


# ============================================================
# 35. SAVE MATCH SUMMARIES
# ============================================================

readr::write_csv(
  
  resolution_summary,
  
  file.path(
    
    output_folder,
    
    "Identifier_Resolution_Summary.csv"
  )
)


readr::write_csv(
  
  match_summary,
  
  file.path(
    
    output_folder,
    
    "Identifier_Match_Status_Summary.csv"
  )
)


readr::write_csv(
  
  match_method_summary,
  
  file.path(
    
    output_folder,
    
    "Identifier_Match_Method_Summary.csv"
  )
)


readr::write_csv(
  
  critical_resolution,
  
  file.path(
    
    output_folder,
    
    "Critical_Interaction_Identifier_Audit.csv"
  )
)


# ============================================================
# 36. SAVE MASTER LOOKUP FOR REPRODUCIBILITY
# ============================================================

readr::write_csv(
  
  master_lookup,
  
  file.path(
    
    output_folder,
    
    "Master_KB_Identifier_Lookup.csv"
  )
)


# ============================================================
# 37. OUTPUT VERIFICATION
# ============================================================

expected_outputs <- c(
  
  file.path(
    output_folder,
    "Evidence_Layered_Network_ID_Resolved.csv"
  ),
  
  file.path(
    output_folder,
    "Unique_Target_Identifier_Resolution.csv"
  ),
  
  file.path(
    output_folder,
    "Unmatched_Targets.csv"
  ),
  
  file.path(
    output_folder,
    "Ambiguous_Target_Matches.csv"
  ),
  
  file.path(
    output_folder,
    "Identifier_Resolution_Summary.csv"
  ),
  
  file.path(
    output_folder,
    "Identifier_Match_Method_Summary.csv"
  ),
  
  file.path(
    output_folder,
    "Critical_Interaction_Identifier_Audit.csv"
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
    
    output_folder,
    
    "Script16B_Output_Verification.csv"
  )
)


# ============================================================
# 38. SESSION INFORMATION
# ============================================================

writeLines(
  
  capture.output(
    sessionInfo()
  ),
  
  con =
    file.path(
      
      output_folder,
      
      "sessionInfo.txt"
    )
)


# ============================================================
# 39. FINAL REPORT
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 16B COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================\n\n"
)


cat(
  "IDENTIFIER RESOLUTION SUMMARY\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  resolution_summary,
  n = Inf
)


cat(
  "\nMATCH STATUS SUMMARY\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  match_summary,
  n = Inf
)


cat(
  "\nMATCH METHOD SUMMARY\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  match_method_summary,
  n = Inf
)


cat(
  "\nCRITICAL INTERACTION IDENTIFIER AUDIT\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  critical_resolution,
  n = Inf,
  width = Inf
)


cat(
  "\nUNMATCHED TARGETS\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  unmatched_targets,
  n = Inf,
  width = Inf
)


cat(
  "\nOUTPUT VERIFICATION\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  output_verification,
  n = Inf,
  width = Inf
)


cat(
  "\n============================================\n"
)

cat(
  "IMPORTANT\n"
)

cat(
  "============================================\n"
)

cat(
  paste(
    "\nThis script performs identifier resolution only.",
    "\nIt does not alter CLASH evidence scores.",
    "\nIt does not alter transcriptomic evidence.",
    "\nIt does not create biological-priority scores.",
    "\nReview the identifier match rate before rebuilding",
    "functional annotation.\n"
  )
)


cat(
  "\nMain resolved network:\n",
  file.path(
    output_folder,
    "Evidence_Layered_Network_ID_Resolved.csv"
  ),
  "\n"
)


cat(
  "\n============================================\n"
)

cat(
  "END OF SCRIPT 16B\n"
)

cat(
  "============================================\n"
)