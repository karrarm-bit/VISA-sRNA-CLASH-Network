# ============================================================
# SCRIPT 17D v2
# Final Evidence Integration for Putative 3'UTR-associated RNAs
#
# Staphylococcus aureus JKD6008
#
# Candidates:
#   3UTR-RS04120
#   3UTR-RS05585
#
# Evidence layers:
#   1. Genomic context
#   2. RNase III-CLASH interaction evidence
#   3. Target transcriptomic response
#   4. Term-seq boundary evidence
#
# IMPORTANT:
#   - No integrated biological-priority score
#   - No claim of independent RNA processing
#   - No claim of direct regulation
#   - No causal claim regarding vancomycin response
#   - Preferred terminology:
#       "putative 3'UTR-associated RNA candidate"
# ============================================================


# ============================================================
# 01. PACKAGES
# ============================================================

library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)


# ============================================================
# 02. PROJECT PATHS
# ============================================================

project <- "D:/Bac-sRNA"

outdir <- file.path(
  project,
  "17D_final_3UTR_evidence"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 03. INPUT FILES
# ============================================================

network_file <- file.path(
  project,
  "14_revised_evidence_layered_network",
  "tables",
  "nrdF_Candidate_Interactions.csv"
)

term_file <- file.path(
  project,
  "17A_TermSeq_3UTR_boundary_audit",
  "tables",
  "TermSeq_Candidate_Boundary_Summary.csv"
)

neighborhood_file <- file.path(
  project,
  "17C_true_genomic_neighborhood",
  "Candidate_Closest_Features.csv"
)


# ============================================================
# 04. CHECK INPUT FILES
# ============================================================

input_check <- tibble(
  input = c(
    "Revised network",
    "Term-seq audit",
    "Genomic neighborhood"
  ),
  path = c(
    network_file,
    term_file,
    neighborhood_file
  ),
  exists = c(
    file.exists(network_file),
    file.exists(term_file),
    file.exists(neighborhood_file)
  )
)

cat("\n============================================\n")
cat("INPUT FILE CHECK\n")
cat("============================================\n")

print(
  input_check,
  n = Inf,
  width = Inf
)

if (any(!input_check$exists)) {
  
  stop(
    "One or more required input files are missing."
  )
}


# ============================================================
# 05. READ INPUT FILES
# ============================================================

net <- read_csv(
  network_file,
  show_col_types = FALSE
)

term <- read_csv(
  term_file,
  show_col_types = FALSE
)

neigh <- read_csv(
  neighborhood_file,
  show_col_types = FALSE
)


cat("\nInput dimensions:\n")
cat("Network:", nrow(net), "rows\n")
cat("Term-seq:", nrow(term), "rows\n")
cat("Neighborhood:", nrow(neigh), "rows\n")


# ============================================================
# 06. VERIFY REQUIRED NETWORK COLUMNS
# ============================================================

required_network_columns <- c(
  
  "source_node",
  "target_node",
  "source_label",
  "target_label",
  
  "number_of_supporting_rows",
  "total_hybrid_count",
  
  "best_raw_p_value",
  "best_adjusted_p_value",
  "best_connection_score",
  
  "maximum_number_of_experiments",
  
  "target_vancomycin_response",
  "target_log2_fold_change",
  "target_deseq_adjusted_p",
  "target_is_significant_deg",
  "target_strong_response",
  
  "target_transcriptomic_layer",
  
  "clash_evidence_score",
  "clash_support_category",
  
  "evidence_profile"
)


missing_network_columns <- setdiff(
  required_network_columns,
  names(net)
)


if (length(missing_network_columns) > 0) {
  
  cat("\nMissing network columns:\n")
  
  print(
    missing_network_columns
  )
  
  stop(
    "Required columns are missing from nrdF_Candidate_Interactions.csv."
  )
}


cat("\nNetwork column verification: PASS\n")


# ============================================================
# 07. DEFINE CANDIDATE GENOMIC CONTEXT
#
# IMPORTANT:
# Parent assignments derive from Script 17C genomic
# neighborhood inspection.
#
# RS04120 / RS05585 embedded in the CLASH labels are NOT
# interpreted as literal parent locus tags.
# ============================================================

context <- tribble(
  
  ~candidate_short,
  ~candidate_pattern,
  ~candidate_start,
  ~candidate_end,
  ~candidate_strand,
  ~putative_parent_gene,
  ~putative_parent_locus,
  ~parent_strand,
  ~parent_3prime_coordinate,
  ~distance_from_parent_3prime_nt,
  
  "3UTR-RS04120",
  "RS04120",
  811819,
  811847,
  "+",
  "nrdF",
  "SAA6008_RS03945",
  "+",
  811766,
  53,
  
  "3UTR-RS05585",
  "RS05585",
  1084557,
  1084613,
  "-",
  "folD",
  "SAA6008_RS05375",
  "-",
  1084790,
  177
)


# ============================================================
# 08. VERIFY NEIGHBORHOOD COLUMNS
# ============================================================

required_neighborhood_columns <- c(
  "candidate",
  "current_locus_tag",
  "strand"
)


missing_neighborhood_columns <- setdiff(
  required_neighborhood_columns,
  names(neigh)
)


if (length(missing_neighborhood_columns) > 0) {
  
  cat("\nMissing neighborhood columns:\n")
  
  print(
    missing_neighborhood_columns
  )
  
  stop(
    "Required columns are missing from Script 17C output."
  )
}


# ============================================================
# 09. VERIFY PARENT ASSIGNMENTS AGAINST 17C
# ============================================================

context <- context %>%
  
  rowwise() %>%
  
  mutate(
    
    parent_found_in_17C =
      any(
        neigh$candidate == candidate_short &
          neigh$current_locus_tag == putative_parent_locus,
        na.rm = TRUE
      ),
    
    parent_strand_confirmed =
      any(
        neigh$candidate == candidate_short &
          neigh$current_locus_tag == putative_parent_locus &
          neigh$strand == parent_strand,
        na.rm = TRUE
      )
  ) %>%
  
  ungroup()


cat("\n============================================\n")
cat("GENOMIC PARENT VERIFICATION\n")
cat("============================================\n")

print(
  context %>%
    select(
      candidate_short,
      putative_parent_gene,
      putative_parent_locus,
      candidate_strand,
      parent_strand,
      distance_from_parent_3prime_nt,
      parent_found_in_17C,
      parent_strand_confirmed
    ),
  n = Inf,
  width = Inf
)


if (any(!context$parent_found_in_17C)) {
  
  warning(
    "One or more proposed parent genes were not recovered in Script 17C."
  )
}


# ============================================================
# 10. IDENTIFY CANDIDATES IN REVISED NETWORK
# ============================================================

net_candidates <- net %>%
  
  mutate(
    
    candidate_short =
      case_when(
        
        str_detect(
          source_label,
          regex(
            "RS04120",
            ignore_case = TRUE
          )
        ) ~ "3UTR-RS04120",
        
        str_detect(
          source_label,
          regex(
            "RS05585",
            ignore_case = TRUE
          )
        ) ~ "3UTR-RS05585",
        
        TRUE ~ NA_character_
      )
  ) %>%
  
  filter(
    !is.na(candidate_short)
  )


cat("\n============================================\n")
cat("NETWORK CANDIDATE RECOVERY\n")
cat("============================================\n")

network_recovery <- net_candidates %>%
  
  count(
    candidate_short,
    name = "number_of_network_edges"
  )


print(
  network_recovery,
  n = Inf,
  width = Inf
)


expected_candidates <- c(
  "3UTR-RS04120",
  "3UTR-RS05585"
)


missing_candidates <- setdiff(
  expected_candidates,
  network_recovery$candidate_short
)


if (length(missing_candidates) > 0) {
  
  stop(
    paste(
      "Candidate(s) missing from revised network:",
      paste(
        missing_candidates,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 11. RETAIN ACTUAL SCRIPT-14 NETWORK COLUMNS
# ============================================================

net_selected <- net_candidates %>%
  
  select(
    
    candidate_short,
    
    source_node,
    source_label,
    
    target_node,
    target_label,
    
    number_of_supporting_rows,
    total_hybrid_count,
    
    best_raw_p_value,
    best_adjusted_p_value,
    best_connection_score,
    
    maximum_number_of_experiments,
    
    target_vancomycin_response,
    target_log2_fold_change,
    target_deseq_adjusted_p,
    target_is_significant_deg,
    target_strong_response,
    
    target_transcriptomic_layer,
    
    clash_evidence_score,
    clash_support_category,
    
    evidence_profile
  )


# ============================================================
# 12. VERIFY THAT BOTH CANDIDATES TARGET nrdF
# ============================================================

cat("\n============================================\n")
cat("CANDIDATE TARGET CHECK\n")
cat("============================================\n")

print(
  net_selected %>%
    select(
      candidate_short,
      source_label,
      target_node,
      target_label,
      total_hybrid_count,
      maximum_number_of_experiments,
      best_connection_score
    ),
  n = Inf,
  width = Inf
)


# ============================================================
# 13. TERM-SEQ COLUMN CHECK
# ============================================================

cat("\n============================================\n")
cat("AVAILABLE TERM-SEQ COLUMNS\n")
cat("============================================\n")

print(
  names(term)
)


if (!"candidate" %in% names(term)) {
  
  stop(
    "Column 'candidate' is missing from Term-seq summary."
  )
}


# ============================================================
# 14. STANDARDIZE TERM-SEQ CANDIDATE NAMES
# ============================================================

term_candidates <- term %>%
  
  mutate(
    
    candidate_short =
      case_when(
        
        str_detect(
          candidate,
          regex(
            "RS04120",
            ignore_case = TRUE
          )
        ) ~ "3UTR-RS04120",
        
        str_detect(
          candidate,
          regex(
            "RS05585",
            ignore_case = TRUE
          )
        ) ~ "3UTR-RS05585",
        
        TRUE ~ NA_character_
      )
  ) %>%
  
  filter(
    !is.na(candidate_short)
  )


# ============================================================
# 15. TERM-SEQ HELPER
# ============================================================

get_term_column <- function(
    df,
    candidates
) {
  
  hit <- candidates[
    candidates %in% names(df)
  ]
  
  if (length(hit) == 0) {
    return(NULL)
  }
  
  hit[1]
}


signal_col <- get_term_column(
  term_candidates,
  c(
    "replicates_with_signal",
    "n_replicates_with_signal",
    "signal_replicates"
  )
)


enrichment_col <- get_term_column(
  term_candidates,
  c(
    "replicates_with_enrichment",
    "n_replicates_with_enrichment",
    "enriched_replicates"
  )
)


strong_peak_col <- get_term_column(
  term_candidates,
  c(
    "replicates_with_strong_local_peak",
    "n_replicates_with_strong_local_peak",
    "strong_peak_replicates"
  )
)


median_signal_col <- get_term_column(
  term_candidates,
  c(
    "median_boundary_signal",
    "median_signal"
  )
)


median_ratio_col <- get_term_column(
  term_candidates,
  c(
    "median_enrichment_ratio",
    "median_ratio"
  )
)


boundary_class_col <- get_term_column(
  term_candidates,
  c(
    "TermSeq_boundary_class",
    "termseq_boundary_class",
    "boundary_class"
  )
)


# ============================================================
# 16. CREATE STANDARD TERM-SEQ VARIABLES
# ============================================================

term_candidates$termseq_signal_replicates <-
  if (!is.null(signal_col)) {
    suppressWarnings(
      as.numeric(
        term_candidates[[signal_col]]
      )
    )
  } else {
    NA_real_
  }


term_candidates$termseq_enriched_replicates <-
  if (!is.null(enrichment_col)) {
    suppressWarnings(
      as.numeric(
        term_candidates[[enrichment_col]]
      )
    )
  } else {
    NA_real_
  }


term_candidates$termseq_strong_peak_replicates <-
  if (!is.null(strong_peak_col)) {
    suppressWarnings(
      as.numeric(
        term_candidates[[strong_peak_col]]
      )
    )
  } else {
    NA_real_
  }


term_candidates$termseq_median_boundary_signal <-
  if (!is.null(median_signal_col)) {
    suppressWarnings(
      as.numeric(
        term_candidates[[median_signal_col]]
      )
    )
  } else {
    NA_real_
  }


term_candidates$termseq_median_enrichment_ratio <-
  if (!is.null(median_ratio_col)) {
    suppressWarnings(
      as.numeric(
        term_candidates[[median_ratio_col]]
      )
    )
  } else {
    NA_real_
  }


term_candidates$termseq_boundary_class <-
  if (!is.null(boundary_class_col)) {
    as.character(
      term_candidates[[boundary_class_col]]
    )
  } else {
    NA_character_
  }


# ============================================================
# 17. REDUCE TERM-SEQ TABLE
# ============================================================

term_selected <- term_candidates %>%
  
  select(
    candidate_short,
    termseq_signal_replicates,
    termseq_enriched_replicates,
    termseq_strong_peak_replicates,
    termseq_median_boundary_signal,
    termseq_median_enrichment_ratio,
    termseq_boundary_class
  )


cat("\n============================================\n")
cat("TERM-SEQ CANDIDATE EVIDENCE\n")
cat("============================================\n")

print(
  term_selected,
  n = Inf,
  width = Inf
)


# ============================================================
# 18. INTEGRATE ALL EVIDENCE
# ============================================================

final <- context %>%
  
  left_join(
    net_selected,
    by = "candidate_short"
  ) %>%
  
  left_join(
    term_selected,
    by = "candidate_short"
  )


# ============================================================
# 19. CHECK FOR DUPLICATION AFTER JOINS
# ============================================================

duplication_check <- final %>%
  
  count(
    candidate_short,
    name = "rows_after_join"
  )


cat("\n============================================\n")
cat("JOIN CHECK\n")
cat("============================================\n")

print(
  duplication_check,
  n = Inf,
  width = Inf
)


if (any(duplication_check$rows_after_join > 1)) {
  
  warning(
    paste(
      "One or more candidates have multiple rows after integration.",
      "Inspect network or Term-seq inputs before manuscript use."
    )
  )
}


# ============================================================
# 20. GENOMIC-CONTEXT INTERPRETATION
# ============================================================

final <- final %>%
  
  mutate(
    
    genomic_context_class =
      case_when(
        
        parent_found_in_17C &
          parent_strand_confirmed &
          distance_from_parent_3prime_nt <= 100 ~
          
          "Strongly compatible with 3′UTR-associated location",
        
        parent_found_in_17C &
          parent_strand_confirmed &
          distance_from_parent_3prime_nt <= 250 ~
          
          "Compatible with extended 3′UTR-associated location",
        
        TRUE ~
          
          "Uncertain genomic-context support"
      ),
    
    genomic_context_interpretation =
      case_when(
        
        candidate_short == "3UTR-RS04120" ~
          
          paste0(
            "The candidate lies ",
            distance_from_parent_3prime_nt,
            " nt from the annotated nrdF 3′ end ",
            "on the same strand."
          ),
        
        candidate_short == "3UTR-RS05585" ~
          
          paste0(
            "The candidate lies ",
            distance_from_parent_3prime_nt,
            " nt from the annotated folD 3′ end ",
            "on the same strand."
          ),
        
        TRUE ~
          
          NA_character_
      )
  )


# ============================================================
# 21. CLASH INTERPRETATION
# ============================================================

final <- final %>%
  
  mutate(
    
    clash_interpretation =
      paste0(
        "RNase III-CLASH recovered the candidate-nrdF interaction ",
        "with ",
        total_hybrid_count,
        " hybrid(s), ",
        maximum_number_of_experiments,
        " experiment(s), and a best connection score of ",
        round(
          best_connection_score,
          3
        ),
        "."
      )
  )


# ============================================================
# 22. TRANSCRIPTOMIC INTERPRETATION
# ============================================================

final <- final %>%
  
  mutate(
    
    target_transcriptomic_interpretation =
      case_when(
        
        !is.na(target_deseq_adjusted_p) &
          target_deseq_adjusted_p < 0.05 &
          !is.na(target_log2_fold_change) &
          target_log2_fold_change < 0 ~
          
          paste0(
            "The target transcript is significantly decreased ",
            "under vancomycin exposure ",
            "(log2FC = ",
            round(
              target_log2_fold_change,
              3
            ),
            "; adjusted P = ",
            format(
              target_deseq_adjusted_p,
              scientific = TRUE,
              digits = 3
            ),
            ")."
          ),
        
        !is.na(target_deseq_adjusted_p) &
          target_deseq_adjusted_p < 0.05 &
          !is.na(target_log2_fold_change) &
          target_log2_fold_change > 0 ~
          
          paste0(
            "The target transcript is significantly increased ",
            "under vancomycin exposure ",
            "(log2FC = ",
            round(
              target_log2_fold_change,
              3
            ),
            "; adjusted P = ",
            format(
              target_deseq_adjusted_p,
              scientific = TRUE,
              digits = 3
            ),
            ")."
          ),
        
        !is.na(target_deseq_adjusted_p) &
          target_deseq_adjusted_p < 0.05 ~
          
          paste(
            "The target transcript shows a statistically",
            "significant response under vancomycin exposure."
          ),
        
        !is.na(target_deseq_adjusted_p) ~
          
          paste(
            "The target transcript was measured but did not",
            "meet the predefined significance threshold."
          ),
        
        TRUE ~
          
          paste(
            "Mapped transcriptomic evidence was not available."
          )
      )
  )


# ============================================================
# 23. TERM-SEQ INTERPRETATION
# ============================================================

final <- final %>%
  
  mutate(
    
    termseq_interpretation =
      case_when(
        
        !is.na(termseq_strong_peak_replicates) &
          termseq_strong_peak_replicates >= 2 ~
          
          paste(
            "A reproducible sharp Term-seq peak is present",
            "near the candidate boundary.",
            "This supports a discrete 3′ end but does not",
            "alone establish independent RNA processing."
          ),
        
        !is.na(termseq_enriched_replicates) &
          termseq_enriched_replicates >= 2 ~
          
          paste(
            "Reproducible local Term-seq enrichment is present",
            "near the candidate boundary, but a reproducible",
            "sharp 3′-end peak is not established."
          ),
        
        !is.na(termseq_signal_replicates) &
          termseq_signal_replicates >= 2 ~
          
          paste(
            "Term-seq signal is detectable in both replicates,",
            "but reproducible enrichment or a sharp 3′ boundary",
            "is not established."
          ),
        
        !is.na(termseq_signal_replicates) &
          termseq_signal_replicates == 1 ~
          
          paste(
            "Term-seq signal is detected in only one replicate;",
            "this is insufficient evidence for a reproducible",
            "discrete 3′ boundary."
          ),
        
        TRUE ~
          
          paste(
            "No reproducible Term-seq evidence for a discrete",
            "candidate 3′ boundary was detected."
          )
      )
  )


# ============================================================
# 24. CONSERVATIVE BIOLOGICAL STATUS
# ============================================================

final <- final %>%
  
  mutate(
    
    recommended_nomenclature =
      "putative 3′UTR-associated RNA candidate",
    
    independent_processing_status =
      "Not established",
    
    direct_regulation_status =
      "Not established",
    
    causal_vancomycin_role =
      "Not established"
  )


# ============================================================
# 25. MANUSCRIPT-SAFE CANDIDATE CONCLUSIONS
# ============================================================

final <- final %>%
  
  mutate(
    
    manuscript_conclusion =
      case_when(
        
        candidate_short == "3UTR-RS04120" ~
          
          paste(
            "The RS04120-labelled CLASH RNA fragment lies",
            "53 nt from the annotated nrdF 3′ end on the same",
            "strand, supporting its description as a putative",
            "nrdF 3′UTR-associated RNA candidate.",
            "RNase III-CLASH recovered an interaction with nrdF,",
            "while nrdF showed a significant transcriptional",
            "decrease under vancomycin exposure.",
            "Term-seq did not demonstrate a reproducible sharp",
            "candidate 3′ boundary; therefore independent",
            "processing and direct regulation remain unproven."
          ),
        
        candidate_short == "3UTR-RS05585" ~
          
          paste(
            "The RS05585-labelled CLASH RNA fragment lies",
            "177 nt from the annotated folD 3′ end on the same",
            "strand, consistent with a putative extended folD",
            "3′UTR-associated RNA candidate.",
            "RNase III-CLASH recovered an interaction with nrdF,",
            "while nrdF showed a significant transcriptional",
            "decrease under vancomycin exposure.",
            "Term-seq did not provide reproducible evidence",
            "for a discrete candidate 3′ boundary; therefore",
            "independent processing and direct regulation",
            "remain unproven."
          ),
        
        TRUE ~ NA_character_
      )
  )


# ============================================================
# 26. CLAIM AUDIT
# ============================================================

claim_audit <- final %>%
  
  transmute(
    
    candidate =
      candidate_short,
    
    proposed_parent =
      putative_parent_gene,
    
    target =
      target_label,
    
    genomic_3UTR_context_supported =
      parent_found_in_17C &
      parent_strand_confirmed &
      distance_from_parent_3prime_nt <= 250,
    
    CLASH_interaction_supported =
      !is.na(total_hybrid_count) &
      total_hybrid_count > 0,
    
    target_significant_transcriptional_response =
      !is.na(target_deseq_adjusted_p) &
      target_deseq_adjusted_p < 0.05,
    
    reproducible_sharp_TermSeq_boundary =
      !is.na(termseq_strong_peak_replicates) &
      termseq_strong_peak_replicates >= 2,
    
    independent_processing_demonstrated =
      FALSE,
    
    direct_regulation_demonstrated =
      FALSE,
    
    causal_vancomycin_role_demonstrated =
      FALSE
  )


# ============================================================
# 27. EVIDENCE-LAYER SUMMARY
# ============================================================

evidence_summary <- final %>%
  
  transmute(
    
    candidate =
      candidate_short,
    
    putative_parent =
      putative_parent_gene,
    
    distance_from_parent_3prime_nt =
      distance_from_parent_3prime_nt,
    
    genomic_context =
      genomic_context_class,
    
    CLASH_support =
      clash_support_category,
    
    CLASH_evidence_score =
      clash_evidence_score,
    
    hybrid_count =
      total_hybrid_count,
    
    experiments =
      maximum_number_of_experiments,
    
    connection_score =
      best_connection_score,
    
    target =
      target_label,
    
    target_log2FC =
      target_log2_fold_change,
    
    target_adjusted_P =
      target_deseq_adjusted_p,
    
    TermSeq_signal_replicates =
      termseq_signal_replicates,
    
    TermSeq_enriched_replicates =
      termseq_enriched_replicates,
    
    TermSeq_strong_peak_replicates =
      termseq_strong_peak_replicates,
    
    TermSeq_boundary_class =
      termseq_boundary_class,
    
    independent_processing =
      "Not established",
    
    direct_regulation =
      "Not established"
  )


# ============================================================
# 28. MANUSCRIPT-READY TABLE
# ============================================================

manuscript_table <- final %>%
  
  transmute(
    
    Candidate =
      candidate_short,
    
    Putative_parent_gene =
      putative_parent_gene,
    
    Parent_locus =
      putative_parent_locus,
    
    Candidate_coordinates =
      paste0(
        candidate_start,
        "-",
        candidate_end
      ),
    
    Strand =
      candidate_strand,
    
    Distance_from_parent_3prime_nt =
      distance_from_parent_3prime_nt,
    
    CLASH_target =
      target_label,
    
    Supporting_CLASH_rows =
      number_of_supporting_rows,
    
    Hybrid_count =
      total_hybrid_count,
    
    Experiments =
      maximum_number_of_experiments,
    
    Best_connection_score =
      best_connection_score,
    
    Best_CLASH_adjusted_P =
      best_adjusted_p_value,
    
    CLASH_evidence_score =
      clash_evidence_score,
    
    CLASH_support_category =
      clash_support_category,
    
    Target_log2FC =
      target_log2_fold_change,
    
    Target_DESeq2_adjusted_P =
      target_deseq_adjusted_p,
    
    Target_significant_DEG =
      target_is_significant_deg,
    
    TermSeq_signal_replicates =
      termseq_signal_replicates,
    
    TermSeq_enriched_replicates =
      termseq_enriched_replicates,
    
    TermSeq_strong_peak_replicates =
      termseq_strong_peak_replicates,
    
    TermSeq_boundary_class =
      termseq_boundary_class,
    
    Recommended_nomenclature =
      recommended_nomenclature,
    
    Independent_processing =
      independent_processing_status,
    
    Direct_regulation =
      direct_regulation_status
  )


# ============================================================
# 29. SAVE OUTPUT TABLES
# ============================================================

write_csv(
  final,
  file.path(
    outdir,
    "FINAL_3UTR_Candidate_Evidence_Integration.csv"
  )
)


write_csv(
  manuscript_table,
  file.path(
    outdir,
    "Manuscript_Ready_3UTR_Evidence_Table.csv"
  )
)


write_csv(
  claim_audit,
  file.path(
    outdir,
    "3UTR_Claim_Audit.csv"
  )
)


write_csv(
  evidence_summary,
  file.path(
    outdir,
    "3UTR_Evidence_Layer_Summary.csv"
  )
)


# ============================================================
# 30. SAVE METHODOLOGICAL NOTES
# ============================================================

method_notes <- c(
  
  "SCRIPT 17D v2",
  "FINAL 3'UTR EVIDENCE INTEGRATION",
  
  "",
  
  "Evidence layers were kept analytically separate:",
  
  "1. genomic context",
  "2. RNase III-CLASH interaction evidence",
  "3. target transcriptomic response",
  "4. Term-seq boundary evidence",
  
  "",
  
  paste(
    "The RS04120 and RS05585 components of the candidate",
    "labels were not treated as literal parent locus tags."
  ),
  
  "",
  
  paste(
    "Parent-gene assignments were based on the genomic",
    "neighborhood analysis from Script 17C."
  ),
  
  "",
  
  paste(
    "3UTR-RS04120 is located 53 nt from the annotated",
    "nrdF 3-prime end on the same strand."
  ),
  
  "",
  
  paste(
    "3UTR-RS05585 is located 177 nt from the annotated",
    "folD 3-prime end on the same strand."
  ),
  
  "",
  
  paste(
    "CLASH interaction evidence does not by itself establish",
    "direct regulatory function."
  ),
  
  "",
  
  paste(
    "Target differential expression is treated as an",
    "independent evidence layer and is not incorporated",
    "into the CLASH evidence score."
  ),
  
  "",
  
  paste(
    "Term-seq signal is interpreted as boundary evidence",
    "and not as proof of independent RNA processing."
  ),
  
  "",
  
  paste(
    "Preferred terminology:",
    "putative 3-prime-UTR-associated RNA candidate."
  ),
  
  "",
  
  "No integrated biological-priority score was calculated."
)


writeLines(
  method_notes,
  file.path(
    outdir,
    "Script17D_Methodological_Notes.txt"
  )
)


# ============================================================
# 31. SESSION INFORMATION
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    outdir,
    "sessionInfo.txt"
  )
)


# ============================================================
# 32. FINAL CONSOLE REPORT
# ============================================================

cat("\n\n")
cat("============================================\n")
cat("SCRIPT 17D v2 COMPLETED\n")
cat("============================================\n")


cat("\nFINAL 3′UTR EVIDENCE INTEGRATION\n")
cat("--------------------------------------------\n")


print(
  final %>%
    
    select(
      
      candidate_short,
      
      putative_parent_gene,
      putative_parent_locus,
      
      candidate_strand,
      parent_strand,
      
      distance_from_parent_3prime_nt,
      
      genomic_context_class,
      
      target_label,
      
      number_of_supporting_rows,
      total_hybrid_count,
      maximum_number_of_experiments,
      
      best_connection_score,
      best_adjusted_p_value,
      
      clash_evidence_score,
      clash_support_category,
      
      target_log2_fold_change,
      target_deseq_adjusted_p,
      
      target_is_significant_deg,
      target_strong_response,
      
      termseq_signal_replicates,
      termseq_enriched_replicates,
      termseq_strong_peak_replicates,
      
      termseq_boundary_class,
      
      independent_processing_status,
      direct_regulation_status
    ),
  
  n = Inf,
  width = Inf
)


cat("\n\n")
cat("CLAIM AUDIT\n")
cat("--------------------------------------------\n")


print(
  claim_audit,
  n = Inf,
  width = Inf
)


cat("\n\n")
cat("EVIDENCE-LAYER SUMMARY\n")
cat("--------------------------------------------\n")


print(
  evidence_summary,
  n = Inf,
  width = Inf
)


cat("\n\n")
cat("MANUSCRIPT-SAFE CONCLUSIONS\n")
cat("--------------------------------------------\n")


for (i in seq_len(nrow(final))) {
  
  cat(
    "\n",
    final$candidate_short[i],
    "\n",
    final$manuscript_conclusion[i],
    "\n",
    sep = ""
  )
}


cat("\n\n")
cat("============================================\n")
cat("FINAL INTERPRETATION RULES\n")
cat("============================================\n")


cat(
  paste0(
    
    "\nPreferred terminology:\n",
    "  putative 3′UTR-associated RNA candidate\n",
    
    "\nEvidence supported separately:\n",
    "  - genomic proximity/context\n",
    "  - RNase III-CLASH interaction\n",
    "  - target transcriptomic response\n",
    "  - descriptive Term-seq boundary evidence\n",
    
    "\nNot established:\n",
    "  - independent RNA processing\n",
    "  - direct regulation of nrdF\n",
    "  - causal role in vancomycin response\n",
    
    "\nNo integrated biological-priority score was used.\n"
  )
)


cat("\nOutputs saved to:\n")
cat(outdir, "\n")


cat("\nExpected files:\n")

cat(
  "1. FINAL_3UTR_Candidate_Evidence_Integration.csv\n"
)

cat(
  "2. Manuscript_Ready_3UTR_Evidence_Table.csv\n"
)

cat(
  "3. 3UTR_Claim_Audit.csv\n"
)

cat(
  "4. 3UTR_Evidence_Layer_Summary.csv\n"
)

cat(
  "5. Script17D_Methodological_Notes.txt\n"
)

cat(
  "6. sessionInfo.txt\n"
)


cat("\n============================================\n")
cat("END SCRIPT 17D v2\n")
cat("============================================\n")
list.files(
  "D:/Bac-sRNA/17A_TermSeq_3UTR_boundary_audit",
  recursive = TRUE,
  full.names = TRUE
)