# ============================================================
# BacRegRNA Project
# Script 12 — FIXED
#
# Build the JKD6008 Master Gene Knowledgebase
#
# Central objective:
# Reconstructing the sRNA regulatory network underlying the
# vancomycin response in VISA Staphylococcus aureus
#
# This script:
# 1. Integrates strain-specific annotation with DESeq2 results.
# 2. Parses genomic coordinates.
# 3. Calculates expression summaries safely.
# 4. Classifies vancomycin-responsive genes.
# 5. Creates node and edge templates for sRNA network analysis.
# 6. Preserves the same project folder structure.
# ============================================================


# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

annotation_file <- file.path(
  project_folder,
  "10_annotation",
  "tables",
  "JKD6008_genome_annotation.csv"
)

annotated_results_file <- file.path(
  project_folder,
  "10_annotation",
  "tables",
  "GSE254530_DESeq2_annotated_all_genes.csv"
)

significant_results_file <- file.path(
  project_folder,
  "10_annotation",
  "tables",
  "GSE254530_DESeq2_annotated_significant_genes.csv"
)

normalized_counts_file <- file.path(
  project_folder,
  "09_RNAseq_analysis",
  "tables",
  "GSE254530_normalized_counts.csv"
)

output_folder <- file.path(
  project_folder,
  "12_master_knowledgebase"
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
# 2. Install required packages
# ------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "tidyr",
  "tibble",
  "purrr"
)

missing_packages <- required_packages[
  !required_packages %in%
    rownames(installed.packages())
]

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
library(tidyr)
library(tibble)
library(purrr)


# ------------------------------------------------------------
# 4. Verify required input files
# ------------------------------------------------------------

required_files <- c(
  annotation_file,
  annotated_results_file,
  significant_results_file,
  normalized_counts_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    "Required input files were not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

cat(
  "\nAll required input files were found.\n"
)


# ------------------------------------------------------------
# 5. Read input data
# ------------------------------------------------------------

genome_annotation <- readr::read_csv(
  annotation_file,
  show_col_types = FALSE
)

deseq_results <- readr::read_csv(
  annotated_results_file,
  show_col_types = FALSE
)

significant_results <- readr::read_csv(
  significant_results_file,
  show_col_types = FALSE
)

normalized_counts <- readr::read_csv(
  normalized_counts_file,
  show_col_types = FALSE
)

cat(
  "\nGenome annotation rows:",
  nrow(genome_annotation),
  "\n"
)

cat(
  "DESeq2 result rows:",
  nrow(deseq_results),
  "\n"
)

cat(
  "Significant DEG rows:",
  nrow(significant_results),
  "\n"
)

cat(
  "Normalized-count rows:",
  nrow(normalized_counts),
  "\n"
)


# ------------------------------------------------------------
# 6. Validate essential columns
# ------------------------------------------------------------

required_deseq_columns <- c(
  "gene_id",
  "base_mean",
  "log2_fold_change",
  "p_value",
  "adjusted_p_value",
  "significance"
)

missing_deseq_columns <- setdiff(
  required_deseq_columns,
  names(deseq_results)
)

if (length(missing_deseq_columns) > 0) {
  
  stop(
    "Missing DESeq2 columns: ",
    paste(
      missing_deseq_columns,
      collapse = ", "
    )
  )
}

if (!"gene_id" %in%
    names(genome_annotation)) {
  
  stop(
    "Genome annotation does not contain gene_id."
  )
}

if (!"gene_id" %in%
    names(normalized_counts)) {
  
  stop(
    "Normalized-count table does not contain gene_id."
  )
}


# ------------------------------------------------------------
# 7. Helper: add missing character columns
# ------------------------------------------------------------

add_missing_character_column <- function(
    data,
    column_name
) {
  
  if (!column_name %in% names(data)) {
    
    data[[column_name]] <- NA_character_
  }
  
  data
}

annotation_columns <- c(
  "current_locus_tag",
  "gene_symbol",
  "gene_synonym",
  "product",
  "product_short",
  "protein_id",
  "feature_type",
  "ncRNA_class",
  "feature_location",
  "display_label",
  "annotation_status",
  "annotation_source"
)

for (column_name in annotation_columns) {
  
  genome_annotation <-
    add_missing_character_column(
      genome_annotation,
      column_name
    )
}


# ------------------------------------------------------------
# 8. Helper: parse genomic locations
#
# Supported examples:
# 123..456
# complement(123..456)
# join(123..200,300..456)
# ------------------------------------------------------------

parse_genomic_location <- function(location) {
  
  location <- as.character(location)
  
  strand <- dplyr::case_when(
    
    stringr::str_detect(
      location,
      "^complement\\("
    ) ~ "-",
    
    !is.na(location) ~ "+",
    
    TRUE ~ NA_character_
  )
  
  coordinate_values <- stringr::str_extract_all(
    location,
    "[0-9]+"
  )
  
  start_coordinate <- purrr::map_dbl(
    coordinate_values,
    function(values) {
      
      if (length(values) == 0) {
        
        return(NA_real_)
      }
      
      min(
        as.numeric(values),
        na.rm = TRUE
      )
    }
  )
  
  end_coordinate <- purrr::map_dbl(
    coordinate_values,
    function(values) {
      
      if (length(values) == 0) {
        
        return(NA_real_)
      }
      
      max(
        as.numeric(values),
        na.rm = TRUE
      )
    }
  )
  
  tibble::tibble(
    genomic_start = start_coordinate,
    genomic_end = end_coordinate,
    strand = strand
  )
}


# ------------------------------------------------------------
# 9. Prepare strain-specific genome annotation
# ------------------------------------------------------------

location_table <- parse_genomic_location(
  genome_annotation$feature_location
)

annotation_prepared <- genome_annotation |>
  dplyr::mutate(
    
    gene_id =
      as.character(gene_id),
    
    gene_symbol =
      dplyr::na_if(
        stringr::str_squish(
          gene_symbol
        ),
        ""
      ),
    
    product =
      dplyr::na_if(
        stringr::str_squish(
          product
        ),
        ""
      ),
    
    product_short =
      dplyr::na_if(
        stringr::str_squish(
          product_short
        ),
        ""
      ),
    
    feature_location =
      dplyr::na_if(
        stringr::str_squish(
          feature_location
        ),
        ""
      )
  ) |>
  dplyr::bind_cols(
    location_table
  ) |>
  dplyr::mutate(
    
    gene_length_bp =
      dplyr::case_when(
        
        !is.na(genomic_start) &
          !is.na(genomic_end) ~
          genomic_end -
          genomic_start +
          1,
        
        TRUE ~ NA_real_
      ),
    
    is_rna_feature =
      dplyr::case_when(
        
        feature_type %in% c(
          "ncRNA",
          "misc_RNA",
          "tRNA",
          "rRNA",
          "tmRNA"
        ) ~ TRUE,
        
        !is.na(ncRNA_class) ~ TRUE,
        
        TRUE ~ FALSE
      )
  ) |>
  dplyr::distinct(
    gene_id,
    .keep_all = TRUE
  )

cat(
  "\nUnique genome-annotation records:",
  nrow(annotation_prepared),
  "\n"
)


# ------------------------------------------------------------
# 10. Prepare DESeq2 statistics
# ------------------------------------------------------------

deseq_prepared <- deseq_results |>
  dplyr::transmute(
    
    gene_id =
      as.character(gene_id),
    
    base_mean =
      as.numeric(base_mean),
    
    log2_fold_change =
      as.numeric(log2_fold_change),
    
    unshrunk_log2_fold_change =
      if (
        "unshrunk_log2_fold_change" %in%
        names(deseq_results)
      ) {
        
        as.numeric(
          unshrunk_log2_fold_change
        )
        
      } else {
        
        NA_real_
      },
    
    standard_error =
      if (
        "standard_error" %in%
        names(deseq_results)
      ) {
        
        as.numeric(
          standard_error
        )
        
      } else {
        
        NA_real_
      },
    
    test_statistic =
      if (
        "test_statistic" %in%
        names(deseq_results)
      ) {
        
        as.numeric(
          test_statistic
        )
        
      } else {
        
        NA_real_
      },
    
    p_value =
      as.numeric(p_value),
    
    adjusted_p_value =
      as.numeric(adjusted_p_value),
    
    significance =
      as.character(significance)
  ) |>
  dplyr::distinct(
    gene_id,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# 11. Detect normalized-count sample columns
# ------------------------------------------------------------

expression_columns <- setdiff(
  names(normalized_counts),
  "gene_id"
)

control_columns <- expression_columns[
  stringr::str_detect(
    expression_columns,
    stringr::regex(
      "Cont[123]",
      ignore_case = TRUE
    )
  )
]

treated_columns <- expression_columns[
  stringr::str_detect(
    expression_columns,
    stringr::regex(
      "Van[123]",
      ignore_case = TRUE
    )
  )
]

cat(
  "\nControl columns detected:\n",
  paste(
    control_columns,
    collapse = "\n"
  ),
  "\n"
)

cat(
  "\nTreated columns detected:\n",
  paste(
    treated_columns,
    collapse = "\n"
  ),
  "\n"
)

if (length(control_columns) != 3) {
  
  stop(
    "Expected exactly 3 control columns, but detected ",
    length(control_columns),
    ".\nDetected columns:\n",
    paste(
      control_columns,
      collapse = "\n"
    )
  )
}

if (length(treated_columns) != 3) {
  
  stop(
    "Expected exactly 3 treated columns, but detected ",
    length(treated_columns),
    ".\nDetected columns:\n",
    paste(
      treated_columns,
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------
# 12. Build numeric expression matrices
#
# This section fixes the previous object '.' not found error.
# ------------------------------------------------------------

control_matrix <- normalized_counts |>
  dplyr::select(
    dplyr::all_of(
      control_columns
    )
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ suppressWarnings(
        as.numeric(.x)
      )
    )
  ) |>
  as.matrix()

treated_matrix <- normalized_counts |>
  dplyr::select(
    dplyr::all_of(
      treated_columns
    )
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ suppressWarnings(
        as.numeric(.x)
      )
    )
  ) |>
  as.matrix()

if (nrow(control_matrix) !=
    nrow(normalized_counts)) {
  
  stop(
    "Control matrix row number does not match normalized counts."
  )
}

if (nrow(treated_matrix) !=
    nrow(normalized_counts)) {
  
  stop(
    "Treated matrix row number does not match normalized counts."
  )
}


# ------------------------------------------------------------
# 13. Calculate expression summaries safely
# ------------------------------------------------------------

expression_summary <- tibble::tibble(
  
  gene_id =
    as.character(
      normalized_counts$gene_id
    ),
  
  untreated_mean =
    rowMeans(
      control_matrix,
      na.rm = TRUE
    ),
  
  treated_mean =
    rowMeans(
      treated_matrix,
      na.rm = TRUE
    ),
  
  untreated_median =
    apply(
      control_matrix,
      1,
      stats::median,
      na.rm = TRUE
    ),
  
  treated_median =
    apply(
      treated_matrix,
      1,
      stats::median,
      na.rm = TRUE
    )
)

if (nrow(expression_summary) !=
    nrow(normalized_counts)) {
  
  stop(
    "Expression-summary row count does not match normalized counts."
  )
}

if (anyDuplicated(
  expression_summary$gene_id
) > 0) {
  
  stop(
    "Duplicated gene IDs were detected in expression_summary."
  )
}

cat(
  "\nExpression summaries calculated for",
  nrow(expression_summary),
  "genes.\n"
)

print(
  utils::head(
    expression_summary
  )
)


# ------------------------------------------------------------
# 14. Create concise master labels
# ------------------------------------------------------------

make_master_label <- function(
    gene_id,
    gene_symbol,
    product_short,
    product
) {
  
  gene_symbol <- dplyr::na_if(
    stringr::str_squish(
      as.character(gene_symbol)
    ),
    ""
  )
  
  product_short <- dplyr::na_if(
    stringr::str_squish(
      as.character(product_short)
    ),
    ""
  )
  
  product <- dplyr::na_if(
    stringr::str_squish(
      as.character(product)
    ),
    ""
  )
  
  selected_product <- dplyr::coalesce(
    product_short,
    product
  )
  
  selected_product <- dplyr::case_when(
    
    !is.na(selected_product) &
      stringr::str_detect(
        selected_product,
        stringr::regex(
          "^hypothetical protein",
          ignore_case = TRUE
        )
      ) ~ NA_character_,
    
    !is.na(selected_product) &
      stringr::str_detect(
        selected_product,
        stringr::regex(
          "^uncharacterized protein",
          ignore_case = TRUE
        )
      ) ~ NA_character_,
    
    TRUE ~ selected_product
  )
  
  label <- dplyr::coalesce(
    gene_symbol,
    selected_product,
    gene_id
  )
  
  dplyr::if_else(
    
    stringr::str_length(label) >
      45,
    
    paste0(
      stringr::str_sub(
        label,
        1,
        42
      ),
      "..."
    ),
    
    label
  )
}


# ------------------------------------------------------------
# 15. Build master gene knowledgebase
# ------------------------------------------------------------

master_gene_kb <- annotation_prepared |>
  dplyr::full_join(
    deseq_prepared,
    by = "gene_id"
  ) |>
  dplyr::left_join(
    expression_summary,
    by = "gene_id"
  ) |>
  dplyr::mutate(
    
    master_label =
      make_master_label(
        gene_id,
        gene_symbol,
        product_short,
        product
      ),
    
    vancomycin_response =
      dplyr::case_when(
        
        significance ==
          "Upregulated" ~
          "Induced",
        
        significance ==
          "Downregulated" ~
          "Repressed",
        
        significance ==
          "Not significant" ~
          "Not significant",
        
        TRUE ~
          "Not tested"
      ),
    
    is_significant_deg =
      dplyr::case_when(
        
        significance %in% c(
          "Upregulated",
          "Downregulated"
        ) ~ TRUE,
        
        TRUE ~ FALSE
      ),
    
    strong_response =
      dplyr::case_when(
        
        !is.na(adjusted_p_value) &
          adjusted_p_value < 0.01 &
          abs(log2_fold_change) >= 2 ~
          TRUE,
        
        TRUE ~ FALSE
      ),
    
    response_priority =
      dplyr::case_when(
        
        strong_response &
          vancomycin_response ==
          "Induced" ~
          "High-priority induced",
        
        strong_response &
          vancomycin_response ==
          "Repressed" ~
          "High-priority repressed",
        
        is_significant_deg ~
          "Significant response",
        
        significance ==
          "Not significant" ~
          "Background",
        
        TRUE ~
          "Not evaluated"
      ),
    
    # Functional annotation fields
    go_terms =
      NA_character_,
    
    kegg_gene_id =
      NA_character_,
    
    kegg_orthology =
      NA_character_,
    
    kegg_pathways =
      NA_character_,
    
    cog_category =
      NA_character_,
    
    cog_description =
      NA_character_,
    
    operon_id =
      NA_character_,
    
    regulon =
      NA_character_,
    
    transcription_factor =
      NA_character_,
    
    resistance_database_hit =
      NA_character_,
    
    resistance_gene =
      NA_character_,
    
    resistance_mechanism =
      NA_character_,
    
    virulence_database_hit =
      NA_character_,
    
    virulence_factor =
      NA_character_,
    
    essentiality_status =
      NA_character_,
    
    # sRNA-network fields
    is_srna =
      is_rna_feature,
    
    srna_standard_name =
      dplyr::case_when(
        
        is_rna_feature &
          !is.na(gene_symbol) ~
          gene_symbol,
        
        TRUE ~
          NA_character_
      ),
    
    clash_target_status =
      FALSE,
    
    clash_srna_status =
      FALSE,
    
    number_of_clash_srna_regulators =
      0L,
    
    number_of_clash_mrna_targets =
      0L,
    
    clash_srna_regulators =
      NA_character_,
    
    clash_mrna_targets =
      NA_character_,
    
    clash_total_hybrid_count =
      NA_real_,
    
    clash_best_adjusted_p =
      NA_real_,
    
    clash_best_connection_score =
      NA_real_,
    
    clash_evidence_level =
      "No CLASH evidence integrated yet",
    
    predicted_srna_target_status =
      FALSE,
    
    predicted_srna_regulators =
      NA_character_,
    
    srna_network_role =
      "Not yet classified",
    
    transcriptomic_evidence =
      dplyr::case_when(
        
        is_significant_deg ~
          "GSE254530 DESeq2",
        
        significance ==
          "Not significant" ~
          paste(
            "Tested in GSE254530;",
            "not significant"
          ),
        
        TRUE ~
          "No transcriptomic evidence"
      ),
    
    interaction_evidence =
      paste(
        "Pending GSE254532",
        "RNase III-CLASH integration"
      ),
    
    functional_validation_evidence =
      "Pending GSE158830 integration",
    
    evidence_tier =
      dplyr::case_when(
        
        is_significant_deg ~
          "Tier 3: transcriptomic evidence",
        
        TRUE ~
          "Unclassified"
      )
  ) |>
  dplyr::arrange(
    
    dplyr::desc(
      is_significant_deg
    ),
    
    adjusted_p_value,
    
    dplyr::desc(
      abs(log2_fold_change)
    ),
    
    gene_id
  )


# ------------------------------------------------------------
# 16. Check unique gene identifiers
# ------------------------------------------------------------

duplicate_gene_ids <- master_gene_kb |>
  dplyr::group_by(
    gene_id
  ) |>
  dplyr::summarise(
    number_of_rows =
      dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::filter(
    number_of_rows > 1
  )

if (nrow(duplicate_gene_ids) > 0) {
  
  readr::write_csv(
    duplicate_gene_ids,
    file.path(
      audit_folder,
      "duplicate_gene_ids.csv"
    )
  )
  
  stop(
    "Duplicated gene IDs remain in the master knowledgebase.\n",
    "See duplicate_gene_ids.csv."
  )
}


# ------------------------------------------------------------
# 17. Create focused DEG table
# ------------------------------------------------------------

vancomycin_degs <- master_gene_kb |>
  dplyr::filter(
    is_significant_deg
  ) |>
  dplyr::arrange(
    vancomycin_response,
    adjusted_p_value,
    dplyr::desc(
      abs(log2_fold_change)
    )
  )


# ------------------------------------------------------------
# 18. Create high-priority response table
# ------------------------------------------------------------

high_priority_response <- master_gene_kb |>
  dplyr::filter(
    strong_response
  ) |>
  dplyr::arrange(
    adjusted_p_value,
    dplyr::desc(
      abs(log2_fold_change)
    )
  )


# ------------------------------------------------------------
# 19. Create RNA-feature table
# ------------------------------------------------------------

rna_features <- master_gene_kb |>
  dplyr::filter(
    is_rna_feature
  ) |>
  dplyr::arrange(
    feature_type,
    genomic_start
  )


# ------------------------------------------------------------
# 20. Create network node template
# ------------------------------------------------------------

network_node_template <- master_gene_kb |>
  dplyr::transmute(
    
    node_id =
      gene_id,
    
    node_label =
      master_label,
    
    node_type =
      dplyr::case_when(
        
        is_rna_feature ~
          "RNA feature",
        
        TRUE ~
          "mRNA/protein-coding gene"
      ),
    
    gene_symbol,
    
    product,
    
    vancomycin_response,
    
    log2_fold_change,
    
    adjusted_p_value,
    
    is_significant_deg,
    
    strong_response,
    
    srna_network_role,
    
    evidence_tier
  )


# ------------------------------------------------------------
# 21. Create empty network edge template
# ------------------------------------------------------------

network_edge_template <- tibble::tibble(
  
  source_node =
    character(),
  
  target_node =
    character(),
  
  source_type =
    character(),
  
  target_type =
    character(),
  
  interaction_type =
    character(),
  
  treatment_context =
    character(),
  
  hybrid_count =
    numeric(),
  
  raw_p_value =
    numeric(),
  
  adjusted_p_value =
    numeric(),
  
  connection_score =
    numeric(),
  
  experimental_dataset =
    character(),
  
  evidence_type =
    character(),
  
  evidence_tier =
    character(),
  
  notes =
    character()
)


# ------------------------------------------------------------
# 22. Knowledgebase summary
# ------------------------------------------------------------

knowledgebase_summary <- tibble::tibble(
  
  metric = c(
    "Total unique genes/features",
    "Protein-coding or non-RNA features",
    "RNA features",
    "Genes tested by DESeq2",
    "Significant vancomycin-responsive genes",
    "Induced genes",
    "Repressed genes",
    "Strong-response genes",
    "Genes with gene symbols",
    "Genes with product annotation",
    "Unmapped or locus-only genes"
  ),
  
  value = c(
    
    nrow(master_gene_kb),
    
    sum(
      !master_gene_kb$is_rna_feature,
      na.rm = TRUE
    ),
    
    sum(
      master_gene_kb$is_rna_feature,
      na.rm = TRUE
    ),
    
    sum(
      !is.na(
        master_gene_kb$base_mean
      )
    ),
    
    sum(
      master_gene_kb$is_significant_deg,
      na.rm = TRUE
    ),
    
    sum(
      master_gene_kb$
        vancomycin_response ==
        "Induced",
      na.rm = TRUE
    ),
    
    sum(
      master_gene_kb$
        vancomycin_response ==
        "Repressed",
      na.rm = TRUE
    ),
    
    sum(
      master_gene_kb$strong_response,
      na.rm = TRUE
    ),
    
    sum(
      !is.na(
        master_gene_kb$gene_symbol
      )
    ),
    
    sum(
      !is.na(
        master_gene_kb$product
      )
    ),
    
    sum(
      master_gene_kb$
        annotation_status %in%
        c(
          "Locus tag only",
          "Unmapped"
        ),
      na.rm = TRUE
    )
  )
)


# ------------------------------------------------------------
# 23. Create field dictionary
# ------------------------------------------------------------

field_dictionary <- tibble::tribble(
  
  ~field_group,
  ~field_name,
  ~description,
  
  "Identity",
  "gene_id",
  "Original SAA6008 locus identifier",
  
  "Identity",
  "current_locus_tag",
  "Current RefSeq locus tag when available",
  
  "Identity",
  "gene_symbol",
  "Standard gene symbol",
  
  "Identity",
  "product",
  "Full strain-specific product annotation",
  
  "Genome",
  "genomic_start",
  "Minimum genomic coordinate parsed from GenBank",
  
  "Genome",
  "genomic_end",
  "Maximum genomic coordinate parsed from GenBank",
  
  "Genome",
  "strand",
  "Genomic strand inferred from the feature location",
  
  "Transcriptomics",
  "log2_fold_change",
  "Vancomycin-treated versus untreated log2 fold change",
  
  "Transcriptomics",
  "adjusted_p_value",
  "Benjamini-Hochberg-adjusted P value",
  
  "Transcriptomics",
  "vancomycin_response",
  "Induced, repressed, not significant, or not tested",
  
  "sRNA network",
  "clash_target_status",
  "Whether the gene is an experimentally detected CLASH target",
  
  "sRNA network",
  "clash_srna_status",
  "Whether the feature acts as an sRNA node in CLASH",
  
  "sRNA network",
  "clash_srna_regulators",
  "Experimentally detected sRNA regulators",
  
  "sRNA network",
  "clash_mrna_targets",
  "Experimentally detected mRNA targets",
  
  "Evidence",
  "evidence_tier",
  "Integrated confidence tier after data integration"
)


# ------------------------------------------------------------
# 24. Save master knowledgebase
# ------------------------------------------------------------

readr::write_csv(
  master_gene_kb,
  file.path(
    table_folder,
    "JKD6008_master_gene_knowledgebase.csv"
  )
)


# ------------------------------------------------------------
# 25. Save focused tables
# ------------------------------------------------------------

readr::write_csv(
  vancomycin_degs,
  file.path(
    table_folder,
    "JKD6008_vancomycin_responsive_DEGs.csv"
  )
)

readr::write_csv(
  high_priority_response,
  file.path(
    table_folder,
    "JKD6008_high_priority_vancomycin_response.csv"
  )
)

readr::write_csv(
  rna_features,
  file.path(
    table_folder,
    "JKD6008_annotated_RNA_features.csv"
  )
)


# ------------------------------------------------------------
# 26. Save network templates
# ------------------------------------------------------------

readr::write_csv(
  network_node_template,
  file.path(
    table_folder,
    "sRNA_network_node_template.csv"
  )
)

readr::write_csv(
  network_edge_template,
  file.path(
    table_folder,
    "sRNA_network_edge_template.csv"
  )
)


# ------------------------------------------------------------
# 27. Save summary and dictionary
# ------------------------------------------------------------

readr::write_csv(
  knowledgebase_summary,
  file.path(
    table_folder,
    "JKD6008_master_knowledgebase_summary.csv"
  )
)

readr::write_csv(
  field_dictionary,
  file.path(
    audit_folder,
    "master_knowledgebase_field_dictionary.csv"
  )
)


# ------------------------------------------------------------
# 28. Save input audit
# ------------------------------------------------------------

input_paths <- c(
  annotation_file,
  annotated_results_file,
  significant_results_file,
  normalized_counts_file
)

input_audit <- tibble::tibble(
  
  input_type = c(
    "Genome annotation",
    "Annotated DESeq2 results",
    "Significant DESeq2 results",
    "Normalized counts"
  ),
  
  input_file =
    input_paths,
  
  exists =
    file.exists(
      input_paths
    ),
  
  size_MB =
    round(
      file.info(
        input_paths
      )$size /
        1024^2,
      3
    )
)

readr::write_csv(
  input_audit,
  file.path(
    audit_folder,
    "Script12_input_audit.csv"
  )
)


# ------------------------------------------------------------
# 29. Save session information
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
# 30. Verify output files
# ------------------------------------------------------------

expected_outputs <- c(
  
  file.path(
    table_folder,
    "JKD6008_master_gene_knowledgebase.csv"
  ),
  
  file.path(
    table_folder,
    "JKD6008_vancomycin_responsive_DEGs.csv"
  ),
  
  file.path(
    table_folder,
    "JKD6008_high_priority_vancomycin_response.csv"
  ),
  
  file.path(
    table_folder,
    "JKD6008_annotated_RNA_features.csv"
  ),
  
  file.path(
    table_folder,
    "sRNA_network_node_template.csv"
  ),
  
  file.path(
    table_folder,
    "sRNA_network_edge_template.csv"
  ),
  
  file.path(
    table_folder,
    "JKD6008_master_knowledgebase_summary.csv"
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
    "Script12_output_verification.csv"
  )
)


# ------------------------------------------------------------
# 31. Final console report
# ------------------------------------------------------------

cat("\n============================================\n")
cat("SCRIPT 12 COMPLETED SUCCESSFULLY\n")
cat("============================================\n\n")

print(
  knowledgebase_summary,
  n = Inf
)

cat(
  "\nOutput verification:\n"
)

print(
  output_verification,
  n = Inf
)

cat(
  "\nMaster knowledgebase:\n",
  file.path(
    table_folder,
    "JKD6008_master_gene_knowledgebase.csv"
  ),
  "\n"
)

cat(
  "\nNetwork node template:\n",
  file.path(
    table_folder,
    "sRNA_network_node_template.csv"
  ),
  "\n"
)

cat(
  "\nNetwork edge template:\n",
  file.path(
    table_folder,
    "sRNA_network_edge_template.csv"
  ),
  "\n"
)

cat(
  "\nNext step:\n",
  "Parse and standardize the GSE254532 ",
  "RNase III-CLASH interaction table.\n"
)