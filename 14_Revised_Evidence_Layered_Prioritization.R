# ============================================================
# BacRegRNA Project
# Script 14 REVISED
# Evidence-Layered sRNA-mRNA Prioritization
#
# PRINCIPLE
# ------------------------------------------------------------
# CLASH evidence defines the interaction network.
#
# Target differential expression is an INDEPENDENT evidence
# layer and is NOT an obligatory gate for network inclusion.
#
# This prevents exclusion of plausible post-transcriptional
# interactions whose target mRNA abundance does not change
# significantly.
#
# IMPORTANT:
# - Original Script 14 is NOT overwritten.
# - Original scoring is preserved as a legacy comparison.
# - New evidence-layered outputs are written separately.
# - Known controls are annotated but receive NO bonus points.
# ============================================================


# ------------------------------------------------------------
# 1. PROJECT PATHS
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

input_folder <- file.path(
  project_folder,
  "13_CLASH_network",
  "tables"
)

collapsed_file <- file.path(
  input_folder,
  "GSE254532_sRNA_mRNA_edges_collapsed.csv"
)

regulatory_summary_file <- file.path(
  input_folder,
  "GSE254532_regulatory_RNA_summary.csv"
)

kb_file <- file.path(
  project_folder,
  "12_master_knowledgebase",
  "tables",
  "JKD6008_master_gene_knowledgebase_CLASH_updated.csv"
)

# fallback if updated KB is absent
kb_fallback <- file.path(
  project_folder,
  "12_master_knowledgebase",
  "tables",
  "JKD6008_master_gene_knowledgebase.csv"
)

output_folder <- file.path(
  project_folder,
  "14_revised_evidence_layered_network"
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
  install.packages(missing_packages)
}

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)


# ------------------------------------------------------------
# 3. VERIFY INPUTS
# ------------------------------------------------------------

if (!file.exists(collapsed_file)) {
  stop(
    paste0(
      "Required collapsed CLASH file not found:\n",
      collapsed_file
    )
  )
}

if (!file.exists(kb_file)) {
  
  if (file.exists(kb_fallback)) {
    
    kb_file <- kb_fallback
    
    warning(
      "CLASH-updated knowledgebase not found. ",
      "Using base master knowledgebase."
    )
    
  } else {
    
    stop(
      "Neither updated nor base master knowledgebase was found."
    )
  }
}


# ------------------------------------------------------------
# 4. READ DATA
# ------------------------------------------------------------

edges <- readr::read_csv(
  collapsed_file,
  show_col_types = FALSE
)

kb <- readr::read_csv(
  kb_file,
  show_col_types = FALSE
)

cat(
  "\n========================================\n",
  "SCRIPT 14 REVISED\n",
  "EVIDENCE-LAYERED PRIORITIZATION\n",
  "========================================\n"
)

cat(
  "\nCollapsed CLASH edges:",
  nrow(edges),
  "\n"
)

cat(
  "Knowledgebase entries:",
  nrow(kb),
  "\n"
)


# ------------------------------------------------------------
# 5. HELPER FUNCTIONS
# ------------------------------------------------------------

safe_chr <- function(x) {
  
  x <- as.character(x)
  
  x[is.na(x)] <- ""
  
  x
}


safe_num <- function(x) {
  
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
  
  y <- str_to_lower(
    str_squish(
      as.character(x)
    )
  )
  
  y %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}


ensure_column <- function(
    data,
    column_name,
    default_value
) {
  
  if (!column_name %in% names(data)) {
    
    data[[column_name]] <-
      rep(
        default_value,
        nrow(data)
      )
  }
  
  data
}


# ------------------------------------------------------------
# 6. ENSURE REQUIRED COLUMNS
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


for (nm in names(required_columns)) {
  
  edges <- ensure_column(
    edges,
    nm,
    required_columns[[nm]]
  )
}


# ------------------------------------------------------------
# 7. STANDARDIZE TYPES
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

for (nm in character_columns) {
  
  edges[[nm]] <-
    safe_chr(
      edges[[nm]]
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

for (nm in numeric_columns) {
  
  edges[[nm]] <-
    safe_num(
      edges[[nm]]
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
# 8. SEARCHABLE TEXT
# ------------------------------------------------------------

edges <- edges |>
  mutate(
    
    source_search_text =
      str_to_lower(
        paste(
          source_node,
          source_label,
          source_raw_names,
          source_type,
          sep = " | "
        )
      ),
    
    target_search_text =
      str_to_lower(
        paste(
          target_node,
          target_label,
          target_raw_names,
          target_type,
          sep = " | "
        )
      )
  )


# ------------------------------------------------------------
# 9. CONSERVATIVE SOURCE-RNA CLASSIFICATION
# ------------------------------------------------------------

edges <- edges |>
  mutate(
    
    source_RNA_class =
      case_when(
        
        str_detect(
          source_search_text,
          regex(
            "rnaIII|rsa[a-z0-9/]*|spr[a-z0-9/]*|teg[0-9a-z_-]*",
            ignore_case = TRUE
          )
        ) ~
          "Named sRNA",
        
        str_detect(
          source_search_text,
          regex(
            "3utr|3'utr|3′utr",
            ignore_case = TRUE
          )
        ) ~
          "Putative 3UTR-associated RNA",
        
        str_detect(
          source_search_text,
          regex(
            "5utr|5'utr|5′utr",
            ignore_case = TRUE
          )
        ) ~
          "Putative 5UTR-associated RNA",
        
        str_detect(
          source_search_text,
          regex(
            "intergenic|igr",
            ignore_case = TRUE
          )
        ) ~
          "Intergenic RNA candidate",
        
        str_detect(
          source_search_text,
          regex(
            "srna[-_ ]?[0-9]+",
            ignore_case = TRUE
          )
        ) ~
          "Numbered sRNA candidate",
        
        TRUE ~
          "Other regulatory RNA candidate"
      )
  )


# ------------------------------------------------------------
# 10. BENCHMARK ANNOTATION
#
# NO score is added for being a known RNA.
# ------------------------------------------------------------

edges <- edges |>
  mutate(
    
    benchmark_annotation =
      case_when(
        
        str_detect(
          source_search_text,
          regex(
            "rsaoi|sau[-_ ]?6477",
            ignore_case = TRUE
          )
        ) ~
          "Known vancomycin-responsive sRNA: RsaOI",
        
        str_detect(
          source_search_text,
          regex(
            "spra2|rsaj",
            ignore_case = TRUE
          )
        ) ~
          "Known sRNA: SprA2/RsaJ",
        
        TRUE ~
          "No benchmark label"
      )
  )


# ------------------------------------------------------------
# 11. TARGET TRANSCRIPTOMIC EVIDENCE LAYER
#
# CRITICAL REVISION:
# target DEG status is annotated but is NOT a network gate.
# ------------------------------------------------------------

edges <- edges |>
  mutate(
    
    target_transcriptomic_layer =
      case_when(
        
        target_strong_response ~
          "Strong vancomycin-responsive target",
        
        target_is_significant_deg ~
          "Significant vancomycin-responsive target",
        
        !is.na(target_log2_fold_change) ~
          "Measured but below DEG threshold",
        
        TRUE ~
          "No mapped transcriptomic evidence"
      )
  )


# ------------------------------------------------------------
# 12. CLASH REPRODUCIBILITY COMPONENTS
#
# Preserve original evidence concepts.
# ------------------------------------------------------------

edges <- edges |>
  mutate(
    
    support_score =
      case_when(
        
        number_of_supporting_rows >= 3 ~ 2,
        
        number_of_supporting_rows >= 2 ~ 1,
        
        TRUE ~ 0
      ),
    
    
    experiment_score =
      case_when(
        
        maximum_number_of_experiments >= 3 ~ 2,
        
        maximum_number_of_experiments >= 2 ~ 1,
        
        TRUE ~ 0
      ),
    
    
    hybrid_score =
      case_when(
        
        total_hybrid_count >= 50 ~ 3,
        
        total_hybrid_count >= 20 ~ 2,
        
        total_hybrid_count >= 5 ~ 1,
        
        TRUE ~ 0
      ),
    
    
    connection_score_component =
      case_when(
        
        best_connection_score >= 0.50 ~ 3,
        
        best_connection_score >= 0.20 ~ 2,
        
        best_connection_score >= 0.05 ~ 1,
        
        TRUE ~ 0
      ),
    
    
    clash_significance_score =
      case_when(
        
        !is.na(best_adjusted_p_value) &
          best_adjusted_p_value <= 1e-10 ~ 3,
        
        !is.na(best_adjusted_p_value) &
          best_adjusted_p_value <= 0.001 ~ 2,
        
        !is.na(best_adjusted_p_value) &
          best_adjusted_p_value < 0.05 ~ 1,
        
        TRUE ~ 0
      )
  )


# ------------------------------------------------------------
# 13. CLASH-ONLY SCORE
#
# This score is deliberately independent of target DEG.
# Maximum = 13
# ------------------------------------------------------------

edges <- edges |>
  mutate(
    
    clash_evidence_score =
      support_score +
      experiment_score +
      hybrid_score +
      connection_score_component +
      clash_significance_score
  )


# ------------------------------------------------------------
# 14. CLASH EVIDENCE CATEGORY
#
# Descriptive evidence bins.
# These do NOT claim functional validation.
# ------------------------------------------------------------

edges <- edges |>
  mutate(
    
    clash_support_category =
      case_when(
        
        clash_evidence_score >= 9 ~
          "Higher CLASH support",
        
        clash_evidence_score >= 6 ~
          "Intermediate CLASH support",
        
        clash_evidence_score >= 3 ~
          "Limited CLASH support",
        
        TRUE ~
          "Minimal CLASH support"
      )
  )


# ------------------------------------------------------------
# 15. TRANSCRIPTOMIC SUPPORT SCORE
#
# Kept separate from CLASH score.
# ------------------------------------------------------------

edges <- edges |>
  mutate(
    
    transcriptomic_support_score =
      case_when(
        
        target_strong_response ~ 4,
        
        target_is_significant_deg ~ 3,
        
        !is.na(target_log2_fold_change) ~ 1,
        
        TRUE ~ 0
      )
  )


# ------------------------------------------------------------
# 16. ANNOTATION SUPPORT
#
# Used descriptively only.
# UTR candidates receive no artificial bonus.
# ------------------------------------------------------------

edges <- edges |>
  mutate(
    
    annotation_support_score =
      case_when(
        
        source_RNA_class == "Named sRNA" ~ 2,
        
        source_RNA_class ==
          "Numbered sRNA candidate" ~ 1,
        
        source_RNA_class ==
          "Intergenic RNA candidate" ~ 1,
        
        TRUE ~ 0
      )
  )


# ------------------------------------------------------------
# 17. LEGACY-COMPARABLE TOTAL SCORE
#
# IMPORTANT:
# Retained only to compare with the previous Script 14.
#
# It must NOT be used as the sole definition of the network.
# ------------------------------------------------------------

edges <- edges |>
  mutate(
    
    legacy_comparable_total_score =
      clash_evidence_score +
      transcriptomic_support_score +
      annotation_support_score,
    
    
    legacy_confidence_category =
      case_when(
        
        legacy_comparable_total_score >= 14 ~
          "High",
        
        legacy_comparable_total_score >= 9 ~
          "Moderate",
        
        legacy_comparable_total_score >= 5 ~
          "Supported",
        
        TRUE ~
          "Exploratory"
      )
  )


# ------------------------------------------------------------
# 18. EVIDENCE-LAYER PROFILE
#
# This becomes the preferred interpretation.
# ------------------------------------------------------------

edges <- edges |>
  mutate(
    
    evidence_profile =
      case_when(
        
        clash_evidence_score >= 6 &
          target_is_significant_deg ~
          
          paste0(
            "CLASH-supported + ",
            "transcriptomically responsive target"
          ),
        
        
        clash_evidence_score >= 6 &
          !target_is_significant_deg &
          !is.na(target_log2_fold_change) ~
          
          paste0(
            "CLASH-supported + ",
            "target measured without qualifying DEG"
          ),
        
        
        clash_evidence_score >= 6 &
          is.na(target_log2_fold_change) ~
          
          paste0(
            "CLASH-supported + ",
            "no mapped target transcriptomic evidence"
          ),
        
        
        target_is_significant_deg ~
          
          paste0(
            "Limited CLASH support + ",
            "transcriptomically responsive target"
          ),
        
        
        !is.na(target_log2_fold_change) ~
          
          paste0(
            "Limited CLASH support + ",
            "target measured without qualifying DEG"
          ),
        
        
        TRUE ~
          
          paste0(
            "CLASH interaction with limited ",
            "integrated supporting evidence"
          )
      )
  )


# ------------------------------------------------------------
# 19. PRIMARY NETWORK
#
# ALL collapsed sRNA-mRNA CLASH edges.
#
# No DEG gating.
# ------------------------------------------------------------

clash_supported_network <- edges


# ------------------------------------------------------------
# 20. VANCOMYCIN-RESPONSIVE TARGET SUBSET
#
# This reproduces the biologically narrower subset,
# but it is explicitly a SUBSET rather than the whole network.
# ------------------------------------------------------------

responsive_target_subset <- edges |>
  filter(
    target_is_significant_deg
  )


# ------------------------------------------------------------
# 21. STRONG-RESPONSE TARGET SUBSET
# ------------------------------------------------------------

strong_response_subset <- edges |>
  filter(
    target_strong_response
  )


# ------------------------------------------------------------
# 22. KNOWN-RNA BENCHMARK SUBSET
# ------------------------------------------------------------

benchmark_subset <- edges |>
  filter(
    benchmark_annotation !=
      "No benchmark label"
  )


# ------------------------------------------------------------
# 23. 3'UTR-ASSOCIATED CANDIDATES
# ------------------------------------------------------------

utr3_candidates <- edges |>
  filter(
    source_RNA_class ==
      "Putative 3UTR-associated RNA"
  ) |>
  arrange(
    desc(target_is_significant_deg),
    desc(clash_evidence_score),
    desc(total_hybrid_count)
  )


# ------------------------------------------------------------
# 24. nrdF CANDIDATE SUBSET
# ------------------------------------------------------------

nrdf_candidates <- edges |>
  filter(
    str_detect(
      target_search_text,
      regex(
        "\\bnrdf\\b",
        ignore_case = TRUE
      )
    )
  ) |>
  arrange(
    desc(clash_evidence_score),
    desc(total_hybrid_count)
  )


# ------------------------------------------------------------
# 25. RsaOI SUBSET
# ------------------------------------------------------------

rsaoi_subset <- edges |>
  filter(
    str_detect(
      source_search_text,
      regex(
        "rsaoi|sau[-_ ]?6477",
        ignore_case = TRUE
      )
    )
  )


# ------------------------------------------------------------
# 26. SprA2/RsaJ SUBSET
# ------------------------------------------------------------

spra2_rsaj_subset <- edges |>
  filter(
    str_detect(
      source_search_text,
      regex(
        "spra2|rsaj",
        ignore_case = TRUE
      )
    )
  )


# ------------------------------------------------------------
# 27. CORE OUTPUT TABLE
# ------------------------------------------------------------

core_columns <- c(
  
  "source_node",
  "source_label",
  "source_raw_names",
  "source_RNA_class",
  
  "target_node",
  "target_label",
  "target_raw_names",
  
  "number_of_supporting_rows",
  "total_hybrid_count",
  "maximum_number_of_experiments",
  
  "best_adjusted_p_value",
  "best_connection_score",
  
  "clash_evidence_score",
  "clash_support_category",
  
  "target_log2_fold_change",
  "target_deseq_adjusted_p",
  "target_is_significant_deg",
  "target_strong_response",
  
  "target_transcriptomic_layer",
  "transcriptomic_support_score",
  
  "benchmark_annotation",
  
  "annotation_support_score",
  
  "legacy_comparable_total_score",
  "legacy_confidence_category",
  
  "evidence_profile"
)


core_columns <- intersect(
  core_columns,
  names(edges)
)


core_network_table <- edges |>
  select(
    all_of(core_columns)
  ) |>
  arrange(
    desc(clash_evidence_score),
    desc(target_is_significant_deg),
    desc(total_hybrid_count)
  )


# ------------------------------------------------------------
# 28. WRITE PRIMARY OUTPUTS
# ------------------------------------------------------------

write_csv(
  
  core_network_table,
  
  file.path(
    table_folder,
    "Evidence_Layered_CLASH_Network_ALL.csv"
  )
)


write_csv(
  
  responsive_target_subset,
  
  file.path(
    table_folder,
    "Vancomycin_Responsive_Target_Subset.csv"
  )
)


write_csv(
  
  strong_response_subset,
  
  file.path(
    table_folder,
    "Strong_Response_Target_Subset.csv"
  )
)


write_csv(
  
  benchmark_subset,
  
  file.path(
    table_folder,
    "Known_sRNA_Benchmark_Subset.csv"
  )
)


write_csv(
  
  utr3_candidates,
  
  file.path(
    table_folder,
    "Putative_3UTR_Associated_RNA_Interactions.csv"
  )
)


write_csv(
  
  nrdf_candidates,
  
  file.path(
    table_folder,
    "nrdF_Candidate_Interactions.csv"
  )
)


write_csv(
  
  rsaoi_subset,
  
  file.path(
    table_folder,
    "RsaOI_CLASH_Interactions.csv"
  )
)


write_csv(
  
  spra2_rsaj_subset,
  
  file.path(
    table_folder,
    "SprA2_RsaJ_CLASH_Interactions.csv"
  )
)


# ------------------------------------------------------------
# 29. EVIDENCE-LAYER SUMMARY
# ------------------------------------------------------------

evidence_summary <- tibble(
  
  metric = c(
    
    "Total collapsed CLASH sRNA-mRNA edges",
    
    "Edges with significant target DEG",
    
    "Edges with strong target response",
    
    "Edges with measured target but below DEG threshold",
    
    "Edges without mapped target transcriptomic evidence",
    
    "Named sRNA edges",
    
    "Putative 3UTR-associated RNA edges",
    
    "RsaOI edges",
    
    "SprA2/RsaJ edges",
    
    "nrdF-targeting edges"
  ),
  
  value = c(
    
    nrow(edges),
    
    sum(
      edges$target_is_significant_deg,
      na.rm = TRUE
    ),
    
    sum(
      edges$target_strong_response,
      na.rm = TRUE
    ),
    
    sum(
      !edges$target_is_significant_deg &
        !is.na(
          edges$target_log2_fold_change
        ),
      na.rm = TRUE
    ),
    
    sum(
      is.na(
        edges$target_log2_fold_change
      )
    ),
    
    sum(
      edges$source_RNA_class ==
        "Named sRNA",
      na.rm = TRUE
    ),
    
    sum(
      edges$source_RNA_class ==
        "Putative 3UTR-associated RNA",
      na.rm = TRUE
    ),
    
    nrow(
      rsaoi_subset
    ),
    
    nrow(
      spra2_rsaj_subset
    ),
    
    nrow(
      nrdf_candidates
    )
  )
)


write_csv(
  
  evidence_summary,
  
  file.path(
    table_folder,
    "Evidence_Layer_Summary.csv"
  )
)


# ------------------------------------------------------------
# 30. CLASH SUPPORT DISTRIBUTION
# ------------------------------------------------------------

support_distribution <- edges |>
  count(
    clash_support_category,
    name = "number_of_edges"
  )


write_csv(
  
  support_distribution,
  
  file.path(
    table_folder,
    "CLASH_Support_Distribution.csv"
  )
)


# ------------------------------------------------------------
# 31. TRANSCRIPTOMIC-LAYER DISTRIBUTION
# ------------------------------------------------------------

transcriptomic_distribution <- edges |>
  count(
    target_transcriptomic_layer,
    name = "number_of_edges"
  )


write_csv(
  
  transcriptomic_distribution,
  
  file.path(
    table_folder,
    "Target_Transcriptomic_Layer_Distribution.csv"
  )
)


# ------------------------------------------------------------
# 32. BENCHMARK VS nrdF COMPARISON
#
# Descriptive only.
# Does NOT claim stronger biological evidence.
# ------------------------------------------------------------

benchmark_nrdf_comparison <- bind_rows(
  
  rsaoi_subset |>
    mutate(
      comparison_group = "RsaOI"
    ),
  
  spra2_rsaj_subset |>
    mutate(
      comparison_group = "SprA2/RsaJ"
    ),
  
  nrdf_candidates |>
    mutate(
      comparison_group = "nrdF candidates"
    )
  
) |>
  
  select(
    comparison_group,
    all_of(core_columns)
  )


write_csv(
  
  benchmark_nrdf_comparison,
  
  file.path(
    table_folder,
    "Benchmark_vs_nrdF_Descriptive_Comparison.csv"
  )
)


# ------------------------------------------------------------
# 33. FIGURE 1
# CLASH support vs target transcriptomic evidence
# ------------------------------------------------------------

plot1_data <- edges |>
  mutate(
    
    target_DE_status =
      ifelse(
        target_is_significant_deg,
        "Significant target DEG",
        "Not significant / unmapped"
      )
  )


p1 <- ggplot(
  plot1_data,
  aes(
    x = clash_evidence_score,
    y = total_hybrid_count,
    shape = target_DE_status
  )
) +
  
  geom_point(
    alpha = 0.65,
    size = 2.5
  ) +
  
  labs(
    
    title =
      "CLASH support and target transcriptomic response",
    
    subtitle =
      paste0(
        "Target differential expression is represented ",
        "as an independent evidence layer"
      ),
    
    x =
      "CLASH evidence score",
    
    y =
      "Total hybrid count",
    
    shape =
      "Target transcriptomic status"
  ) +
  
  theme_bw(
    base_size = 12
  )


ggsave(
  
  file.path(
    figure_folder,
    "CLASH_vs_Target_Transcriptomic_Evidence.png"
  ),
  
  p1,
  
  width = 9,
  height = 6,
  dpi = 400
)


# ------------------------------------------------------------
# 34. FIGURE 2
# Evidence-layer distribution
# ------------------------------------------------------------

p2 <- ggplot(
  transcriptomic_distribution,
  aes(
    x = reorder(
      target_transcriptomic_layer,
      number_of_edges
    ),
    y = number_of_edges
  )
) +
  
  geom_col() +
  
  coord_flip() +
  
  labs(
    
    title =
      "Transcriptomic evidence across CLASH-supported interactions",
    
    x =
      NULL,
    
    y =
      "Number of interactions"
  ) +
  
  theme_bw(
    base_size = 12
  )


ggsave(
  
  file.path(
    figure_folder,
    "Target_Transcriptomic_Evidence_Distribution.png"
  ),
  
  p2,
  
  width = 9,
  height = 6,
  dpi = 400
)


# ------------------------------------------------------------
# 35. FIGURE 3
# Benchmark and nrdF descriptive comparison
# ------------------------------------------------------------

plot3_data <- benchmark_nrdf_comparison |>
  mutate(
    
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


if (nrow(plot3_data) > 0) {
  
  p3 <- ggplot(
    plot3_data,
    aes(
      x = reorder(
        interaction_label,
        total_hybrid_count
      ),
      y = total_hybrid_count
    )
  ) +
    
    geom_col() +
    
    coord_flip() +
    
    facet_wrap(
      ~ comparison_group,
      scales = "free_y"
    ) +
    
    labs(
      
      title =
        "Descriptive CLASH support for benchmark and nrdF interactions",
      
      subtitle =
        "Hybrid counts are not interpreted as functional validation",
      
      x =
        "Interaction",
      
      y =
        "Total CLASH hybrid count"
    ) +
    
    theme_bw(
      base_size = 11
    )
  
  
  ggsave(
    
    file.path(
      figure_folder,
      "Benchmark_vs_nrdF_CLASH_Support.png"
    ),
    
    p3,
    
    width = 11,
    height = 7,
    dpi = 400
  )
}


# ------------------------------------------------------------
# 36. CREATE MANUSCRIPT-READY AUDIT TABLE
# ------------------------------------------------------------

manuscript_audit <- edges |>
  transmute(
    
    regulatory_RNA =
      ifelse(
        source_label == "",
        source_raw_names,
        source_label
      ),
    
    RNA_class =
      source_RNA_class,
    
    target =
      ifelse(
        target_label == "",
        target_node,
        target_label
      ),
    
    supporting_rows =
      number_of_supporting_rows,
    
    CLASH_hybrids =
      total_hybrid_count,
    
    experiments =
      maximum_number_of_experiments,
    
    CLASH_adjusted_P =
      best_adjusted_p_value,
    
    connection_score =
      best_connection_score,
    
    CLASH_evidence_score =
      clash_evidence_score,
    
    CLASH_support =
      clash_support_category,
    
    target_log2FC =
      target_log2_fold_change,
    
    target_DE_adjusted_P =
      target_deseq_adjusted_p,
    
    target_DE =
      target_is_significant_deg,
    
    transcriptomic_layer =
      target_transcriptomic_layer,
    
    benchmark =
      benchmark_annotation,
    
    integrated_evidence_profile =
      evidence_profile
  ) |>
  
  arrange(
    desc(CLASH_evidence_score),
    desc(target_DE),
    desc(CLASH_hybrids)
  )


write_csv(
  
  manuscript_audit,
  
  file.path(
    table_folder,
    "Manuscript_Ready_Evidence_Audit_Table.csv"
  )
)


# ------------------------------------------------------------
# 37. WRITE INTERPRETATION REPORT
# ------------------------------------------------------------

report_file <- file.path(
  output_folder,
  "Evidence_Layered_Analysis_Report.txt"
)


sink(
  report_file
)


cat(
  "BacRegRNA Project\n"
)

cat(
  "Evidence-Layered Regulatory RNA Analysis\n"
)

cat(
  "=========================================\n\n"
)


cat(
  "CORE ANALYTICAL PRINCIPLE\n"
)

cat(
  "-------------------------\n"
)

cat(
  paste0(
    "CLASH-supported RNA-mRNA interactions define the ",
    "interaction network. Target differential expression ",
    "is treated as an independent supporting evidence ",
    "layer rather than an obligatory inclusion criterion. ",
    "This distinction is important because regulatory RNAs ",
    "may influence translation without producing a qualifying ",
    "change in steady-state target mRNA abundance.\n\n"
  )
)


cat(
  "NETWORK SUMMARY\n"
)

cat(
  "---------------\n"
)

print(
  evidence_summary
)


cat(
  "\n\nCLASH SUPPORT DISTRIBUTION\n"
)

cat(
  "--------------------------\n"
)

print(
  support_distribution
)


cat(
  "\n\nTARGET TRANSCRIPTOMIC LAYERS\n"
)

cat(
  "----------------------------\n"
)

print(
  transcriptomic_distribution
)


cat(
  "\n\nBENCHMARK INTERPRETATION\n"
)

cat(
  "------------------------\n"
)

cat(
  paste0(
    "Known RNA identities are used for benchmarking and ",
    "biological context only. No score bonus is assigned ",
    "because an RNA is previously characterized. Failure ",
    "to recover a previously reported RNA-target pair from ",
    "the analyzed raw CLASH table is treated as a dataset ",
    "or recovery limitation and not as evidence against the ",
    "published interaction.\n"
  )
)


cat(
  "\n\nnrdF INTERPRETATION\n"
)

cat(
  "-------------------\n"
)

cat(
  paste0(
    "The two putative 3UTR-associated RNA interactions ",
    "targeting nrdF are retained as candidate interactions ",
    "because they are directly represented in the analyzed ",
    "CLASH data and are accompanied by target transcriptomic ",
    "response. Their inclusion does not establish direct ",
    "regulatory activity, regulatory direction, or causality.\n"
  )
)


sink()


# ------------------------------------------------------------
# 38. SESSION INFO
# ------------------------------------------------------------

capture.output(
  
  sessionInfo(),
  
  file = file.path(
    output_folder,
    "sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 39. FINAL CONSOLE OUTPUT
# ------------------------------------------------------------

cat(
  "\n========================================\n"
)

cat(
  "SCRIPT 14 REVISED COMPLETED\n"
)

cat(
  "========================================\n"
)

cat(
  "\nOutput folder:\n",
  output_folder,
  "\n"
)


cat(
  "\nMOST IMPORTANT OUTPUTS:\n"
)

cat(
  "1. tables/Evidence_Layer_Summary.csv\n"
)

cat(
  "2. tables/Evidence_Layered_CLASH_Network_ALL.csv\n"
)

cat(
  "3. tables/Benchmark_vs_nrdF_Descriptive_Comparison.csv\n"
)

cat(
  "4. tables/nrdF_Candidate_Interactions.csv\n"
)

cat(
  "5. tables/Known_sRNA_Benchmark_Subset.csv\n"
)

cat(
  "6. tables/Manuscript_Ready_Evidence_Audit_Table.csv\n"
)

cat(
  "7. Evidence_Layered_Analysis_Report.txt\n"
)


cat(
  "\nIMPORTANT:\n"
)

cat(
  paste0(
    "Do not delete the original Script 14 outputs. ",
    "The next step is to compare OLD vs REVISED ",
    "prioritization before modifying downstream scripts.\n"
  )
)