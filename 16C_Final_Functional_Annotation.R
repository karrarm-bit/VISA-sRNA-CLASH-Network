# ============================================================
# SCRIPT 16C FINAL
# Database-Grounded Functional Annotation
# Evidence-Layered sRNA-mRNA CLASH Network
# Staphylococcus aureus JKD6008
#
# INPUT:
# Script 16B identifier-resolved network
#
# PRINCIPLES:
# 1. Preserve all 211 CLASH-supported interactions.
# 2. Use resolved JKD6008 master-KB annotations.
# 3. GO / KEGG / COG / regulon / resistance / virulence
#    are annotation layers, NOT interaction-validation layers.
# 4. Target DEG remains independent from CLASH evidence.
# 5. No integrated biological-priority score.
# 6. No Priority A/B/C ranking.
# 7. Unresolved identifiers remain unresolved.
# ============================================================


# ============================================================
# 01. PROJECT PATHS
# ============================================================

project_folder <- "D:/Bac-sRNA"

input_file <- file.path(
  project_folder,
  "16B_identifier_resolution",
  "Evidence_Layered_Network_ID_Resolved.csv"
)

output_folder <- file.path(
  project_folder,
  "16C_final_functional_annotation"
)

table_folder <- file.path(
  output_folder,
  "tables"
)

figure_folder <- file.path(
  output_folder,
  "figures"
)

audit_folder <- file.path(
  output_folder,
  "audit"
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

dir.create(
  audit_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 02. PACKAGES
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "tibble",
  "ggplot2",
  "scales"
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
library(ggplot2)
library(scales)


# ============================================================
# 03. HELPER FUNCTIONS
# ============================================================

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  stringr::str_squish(x)
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
  
  z <- stringr::str_to_lower(
    stringr::str_squish(
      as.character(x)
    )
  )
  
  z %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}


add_missing <- function(
    dat,
    nm,
    default = NA_character_
) {
  
  if (!nm %in% names(dat)) {
    dat[[nm]] <- rep(
      default,
      nrow(dat)
    )
  }
  
  dat
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


has_text <- function(x) {
  
  !is.na(x) &
    stringr::str_squish(
      as.character(x)
    ) != ""
}


# ============================================================
# 04. READ RESOLVED NETWORK
# ============================================================

if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Script 16B resolved network not found:\n",
      input_file
    )
  )
}


dat <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 16C FINAL FUNCTIONAL ANNOTATION\n"
)

cat(
  "============================================\n\n"
)

cat(
  "Input interactions:",
  nrow(dat),
  "\n"
)


# ============================================================
# 05. ENSURE REQUIRED COLUMNS
# ============================================================

required_defaults <- list(
  
  source_node = "",
  source_label = "",
  source_raw_names = "",
  source_RNA_class = "",
  
  target_node = "",
  target_label = "",
  target_raw_names = "",
  
  identifier_match_status = "Unmatched",
  identifier_match_method = "No master KB match",
  
  resolved_gene_symbol = "",
  resolved_product = "",
  
  gene_id = "",
  current_locus_tag = "",
  gene_symbol = "",
  gene_synonym = "",
  protein_id = "",
  display_label = "",
  
  annotation_status = "",
  annotation_source = "",
  
  go_terms = "",
  kegg_gene_id = "",
  kegg_orthology = "",
  kegg_pathways = "",
  cog_category = "",
  cog_description = "",
  
  operon_id = "",
  regulon = "",
  transcription_factor = "",
  
  resistance_database_hit = "",
  resistance_gene = "",
  resistance_mechanism = "",
  
  virulence_database_hit = "",
  virulence_factor = "",
  
  essentiality_status = "",
  
  clash_evidence_score = NA_real_,
  clash_support_category = "",
  
  number_of_supporting_rows = NA_real_,
  total_hybrid_count = NA_real_,
  maximum_number_of_experiments = NA_real_,
  best_adjusted_p_value = NA_real_,
  best_connection_score = NA_real_,
  
  target_log2_fold_change = NA_real_,
  target_deseq_adjusted_p = NA_real_,
  target_is_significant_deg = FALSE,
  target_strong_response = FALSE,
  
  target_transcriptomic_layer = "",
  benchmark_annotation = "",
  evidence_profile = ""
)


for (nm in names(required_defaults)) {
  
  dat <- add_missing(
    dat,
    nm,
    required_defaults[[nm]]
  )
}


# ============================================================
# 06. STANDARDIZE TYPES
# ============================================================

character_columns <- c(
  "source_node",
  "source_label",
  "source_raw_names",
  "source_RNA_class",
  "target_node",
  "target_label",
  "target_raw_names",
  "identifier_match_status",
  "identifier_match_method",
  "resolved_gene_symbol",
  "resolved_product",
  "gene_id",
  "current_locus_tag",
  "gene_symbol",
  "gene_synonym",
  "protein_id",
  "display_label",
  "annotation_status",
  "annotation_source",
  "go_terms",
  "kegg_gene_id",
  "kegg_orthology",
  "kegg_pathways",
  "cog_category",
  "cog_description",
  "operon_id",
  "regulon",
  "transcription_factor",
  "resistance_database_hit",
  "resistance_gene",
  "resistance_mechanism",
  "virulence_database_hit",
  "virulence_factor",
  "essentiality_status",
  "clash_support_category",
  "target_transcriptomic_layer",
  "benchmark_annotation",
  "evidence_profile"
)


for (nm in character_columns) {
  
  dat[[nm]] <- safe_chr(
    dat[[nm]]
  )
}


numeric_columns <- c(
  "clash_evidence_score",
  "number_of_supporting_rows",
  "total_hybrid_count",
  "maximum_number_of_experiments",
  "best_adjusted_p_value",
  "best_connection_score",
  "target_log2_fold_change",
  "target_deseq_adjusted_p"
)


for (nm in numeric_columns) {
  
  dat[[nm]] <- safe_num(
    dat[[nm]]
  )
}


dat$target_is_significant_deg <-
  safe_logical(
    dat$target_is_significant_deg
  )


dat$target_strong_response <-
  safe_logical(
    dat$target_strong_response
  )


# ============================================================
# 07. FINAL TARGET IDENTITY
# ============================================================

dat <- dat |>
  
  dplyr::mutate(
    
    target_gene_final =
      first_nonblank(
        .data$gene_symbol,
        .data$resolved_gene_symbol,
        .data$target_label,
        .data$current_locus_tag,
        .data$gene_id,
        .data$protein_id,
        .data$target_node
      ),
    
    target_locus_final =
      first_nonblank(
        .data$current_locus_tag,
        .data$gene_id,
        .data$target_node
      ),
    
    target_protein_final =
      first_nonblank(
        .data$protein_id,
        .data$target_label
      ),
    
    target_product_final =
      first_nonblank(
        .data$resolved_product,
        .data$display_label,
        .data$target_label,
        .data$target_node
      ),
    
    KB_match =
      stringr::str_detect(
        .data$identifier_match_status,
        "^Matched"
      )
  )


# ============================================================
# 08. ANNOTATION AVAILABILITY FLAGS
# ============================================================

dat <- dat |>
  
  dplyr::mutate(
    
    has_GO =
      has_text(
        .data$go_terms
      ),
    
    has_KEGG =
      has_text(
        .data$kegg_gene_id
      ) |
      has_text(
        .data$kegg_orthology
      ) |
      has_text(
        .data$kegg_pathways
      ),
    
    has_COG =
      has_text(
        .data$cog_category
      ) |
      has_text(
        .data$cog_description
      ),
    
    has_regulon =
      has_text(
        .data$regulon
      ) |
      has_text(
        .data$transcription_factor
      ),
    
    has_resistance_annotation =
      has_text(
        .data$resistance_database_hit
      ) |
      has_text(
        .data$resistance_gene
      ) |
      has_text(
        .data$resistance_mechanism
      ),
    
    has_virulence_annotation =
      has_text(
        .data$virulence_database_hit
      ) |
      has_text(
        .data$virulence_factor
      )
  )


# ============================================================
# 09. BUILD DATABASE ANNOTATION TEXT
#
# Used only to normalize broad functional systems.
# It does NOT create interaction evidence.
# ============================================================

dat <- dat |>
  
  dplyr::mutate(
    
    database_annotation_text =
      stringr::str_to_lower(
        
        paste(
          .data$target_product_final,
          .data$go_terms,
          .data$kegg_orthology,
          .data$kegg_pathways,
          .data$cog_category,
          .data$cog_description,
          sep = " | "
        )
      )
  )


# ============================================================
# 10. BROAD FUNCTIONAL SYSTEM NORMALIZATION
#
# This maps existing product/GO/KEGG/COG descriptions into
# manuscript-friendly broad categories.
#
# It does NOT imply resistance, virulence or regulatory effect.
# ============================================================

dat <- dat |>
  
  dplyr::mutate(
    
    functional_system =
      dplyr::case_when(
        
        !.data$KB_match ~
          "Unresolved identifier",
        
        stringr::str_detect(
          .data$database_annotation_text,
          paste0(
            "peptidoglycan|",
            "cell wall|",
            "cell envelope|",
            "teichoic|",
            "penicillin.binding|",
            "mur[a-z0-9]"
          )
        ) ~
          "Cell wall and envelope",
        
        stringr::str_detect(
          .data$database_annotation_text,
          paste0(
            "ribosom|",
            "translation|",
            "trna|",
            "aminoacyl.trna|",
            "translation factor"
          )
        ) ~
          "Translation and protein synthesis",
        
        stringr::str_detect(
          .data$database_annotation_text,
          paste0(
            "dna replication|",
            "dna repair|",
            "ribonucleotide|",
            "recombination|",
            "chromosome|",
            "dna polymerase"
          )
        ) ~
          "DNA synthesis, replication and repair",
        
        stringr::str_detect(
          .data$database_annotation_text,
          paste0(
            "purine|",
            "pyrimidine|",
            "nucleotide|",
            "inosine monophosphate|",
            "imp dehydrogenase|",
            "guanine|",
            "adenine"
          )
        ) ~
          "Nucleotide metabolism",
        
        stringr::str_detect(
          .data$database_annotation_text,
          paste0(
            "amino acid|",
            "glutamate|",
            "glutamine|",
            "arginine|",
            "histidine|",
            "leucine|",
            "isoleucine|",
            "valine|",
            "nitrogen metabolism"
          )
        ) ~
          "Amino-acid and nitrogen metabolism",
        
        stringr::str_detect(
          .data$database_annotation_text,
          paste0(
            "transporter|",
            "transport|",
            "permease|",
            "abc-type|",
            "abc transporter|",
            "substrate.binding|",
            "uptake"
          )
        ) ~
          "Transport and membrane processes",
        
        stringr::str_detect(
          .data$database_annotation_text,
          paste0(
            "respirat|",
            "oxidase|",
            "cytochrome|",
            "electron transport|",
            "quinol|",
            "redox"
          )
        ) ~
          "Respiration and redox metabolism",
        
        stringr::str_detect(
          .data$database_annotation_text,
          paste0(
            "glycolysis|",
            "carbon metabolism|",
            "carbohydrate|",
            "dihydroxyacetone|",
            "sugar|",
            "central metabolism"
          )
        ) ~
          "Carbon and central metabolism",
        
        stringr::str_detect(
          .data$database_annotation_text,
          paste0(
            "stress response|",
            "heat shock|",
            "oxidative stress|",
            "chaperone|",
            "detoxification"
          )
        ) ~
          "Cellular stress response",
        
        .data$has_virulence_annotation ~
          "Virulence and host interaction",
        
        stringr::str_detect(
          .data$database_annotation_text,
          paste0(
            "virulence|",
            "adhesin|",
            "surface protein|",
            "hemolysin|",
            "toxin|",
            "immunoglobulin.binding"
          )
        ) ~
          "Virulence and host interaction",
        
        stringr::str_detect(
          .data$database_annotation_text,
          paste0(
            "transcription regulator|",
            "transcriptional regulator|",
            "regulatory protein|",
            "two.component|",
            "response regulator|",
            "sigma factor"
          )
        ) ~
          "Transcription and signal regulation",
        
        stringr::str_detect(
          .data$database_annotation_text,
          paste0(
            "protein folding|",
            "protease|",
            "peptidase|",
            "protein quality|",
            "protein turnover"
          )
        ) ~
          "Protein processing and turnover",
        
        TRUE ~
          "Other annotated function"
      )
  )


# ============================================================
# 11. FUNCTIONAL ANNOTATION BASIS
# ============================================================

dat <- dat |>
  
  dplyr::mutate(
    
    functional_annotation_basis =
      dplyr::case_when(
        
        !.data$KB_match ~
          "Unresolved identifier",
        
        .data$has_GO &
          .data$has_KEGG &
          .data$has_COG ~
          "Product + GO + KEGG + COG",
        
        .data$has_GO &
          .data$has_KEGG ~
          "Product + GO + KEGG",
        
        .data$has_GO &
          .data$has_COG ~
          "Product + GO + COG",
        
        .data$has_KEGG &
          .data$has_COG ~
          "Product + KEGG + COG",
        
        .data$has_GO ~
          "Product + GO",
        
        .data$has_KEGG ~
          "Product + KEGG",
        
        .data$has_COG ~
          "Product + COG",
        
        TRUE ~
          "Master-KB product annotation"
      )
  )


# ============================================================
# 12. RESISTANCE EVIDENCE
#
# No resistance inference from generic metabolism/stress.
# ============================================================

dat <- dat |>
  
  dplyr::mutate(
    
    resistance_evidence =
      dplyr::case_when(
        
        !.data$KB_match ~
          "Identifier unresolved",
        
        .data$has_resistance_annotation ~
          "Explicit resistance annotation in master KB",
        
        TRUE ~
          "No explicit resistance annotation"
      ),
    
    resistance_annotation =
      dplyr::case_when(
        
        .data$has_resistance_annotation ~
          first_nonblank(
            .data$resistance_gene,
            .data$resistance_mechanism,
            .data$resistance_database_hit
          ),
        
        TRUE ~
          ""
      )
  )


# ============================================================
# 13. VIRULENCE EVIDENCE
# ============================================================

dat <- dat |>
  
  dplyr::mutate(
    
    virulence_evidence =
      dplyr::case_when(
        
        !.data$KB_match ~
          "Identifier unresolved",
        
        .data$has_virulence_annotation ~
          "Explicit virulence annotation in master KB",
        
        TRUE ~
          "No explicit virulence annotation"
      ),
    
    virulence_annotation =
      dplyr::case_when(
        
        .data$has_virulence_annotation ~
          first_nonblank(
            .data$virulence_factor,
            .data$virulence_database_hit
          ),
        
        TRUE ~
          ""
      )
  )


# ============================================================
# 14. REGULON EVIDENCE
# ============================================================

dat <- dat |>
  
  dplyr::mutate(
    
    regulon_evidence =
      dplyr::case_when(
        
        !.data$KB_match ~
          "Identifier unresolved",
        
        .data$has_regulon ~
          "Mapped regulon/transcription-factor annotation",
        
        TRUE ~
          "No mapped regulon annotation"
      ),
    
    regulon_annotation =
      dplyr::case_when(
        
        .data$has_regulon ~
          first_nonblank(
            .data$regulon,
            .data$transcription_factor
          ),
        
        TRUE ~
          ""
      )
  )


# ============================================================
# 15. TRANSCRIPTOMIC STATUS
# ============================================================

dat <- dat |>
  
  dplyr::mutate(
    
    target_transcriptomic_status =
      dplyr::case_when(
        
        .data$target_strong_response ~
          "Strong target DEG",
        
        .data$target_is_significant_deg ~
          "Significant target DEG",
        
        !is.na(
          .data$target_log2_fold_change
        ) ~
          "Measured but below DEG threshold",
        
        TRUE ~
          "No mapped transcriptomic evidence"
      ),
    
    target_response_direction =
      dplyr::case_when(
        
        .data$target_is_significant_deg &
          .data$target_log2_fold_change > 0 ~
          "Induced",
        
        .data$target_is_significant_deg &
          .data$target_log2_fold_change < 0 ~
          "Repressed",
        
        !is.na(
          .data$target_log2_fold_change
        ) ~
          "Measured non-DEG",
        
        TRUE ~
          "Not mapped"
      )
  )


# ============================================================
# 16. CONSERVATIVE SOURCE RNA TERMINOLOGY
# ============================================================

dat <- dat |>
  
  dplyr::mutate(
    
    source_RNA_class_final =
      dplyr::case_when(
        
        stringr::str_detect(
          .data$source_RNA_class,
          stringr::regex(
            "3UTR",
            ignore_case = TRUE
          )
        ) ~
          "Putative 3′UTR-associated RNA",
        
        stringr::str_detect(
          .data$source_RNA_class,
          stringr::regex(
            "5UTR",
            ignore_case = TRUE
          )
        ) ~
          "Putative 5′UTR-associated RNA",
        
        stringr::str_detect(
          .data$source_RNA_class,
          stringr::regex(
            "Named sRNA",
            ignore_case = TRUE
          )
        ) ~
          "Named sRNA",
        
        stringr::str_detect(
          .data$source_RNA_class,
          stringr::regex(
            "Intergenic",
            ignore_case = TRUE
          )
        ) ~
          "Intergenic RNA candidate",
        
        .data$source_RNA_class != "" ~
          .data$source_RNA_class,
        
        TRUE ~
          "Regulatory RNA candidate"
      )
  )


# ============================================================
# 17. DATABASE ANNOTATION COMPLETENESS
# ============================================================

dat <- dat |>
  
  dplyr::mutate(
    
    database_annotation_layers =
      .data$has_GO +
      .data$has_KEGG +
      .data$has_COG +
      .data$has_regulon +
      .data$has_resistance_annotation +
      .data$has_virulence_annotation,
    
    database_annotation_completeness =
      dplyr::case_when(
        
        !.data$KB_match ~
          "Identifier unresolved",
        
        .data$database_annotation_layers >= 4 ~
          "Extensively annotated",
        
        .data$database_annotation_layers >= 2 ~
          "Multiply annotated",
        
        .data$database_annotation_layers == 1 ~
          "Single structured annotation layer",
        
        TRUE ~
          "Product annotation only"
      )
  )


# ============================================================
# 18. FULL FINAL ANNOTATED NETWORK
# ============================================================

final_network <- dat |>
  
  dplyr::arrange(
    dplyr::desc(
      .data$clash_evidence_score
    ),
    dplyr::desc(
      .data$total_hybrid_count
    )
  )


# ============================================================
# 19. MANUSCRIPT-READY INTERACTION TABLE
# ============================================================

manuscript_table <- final_network |>
  
  dplyr::transmute(
    
    regulatory_RNA =
      .data$source_label,
    
    RNA_class =
      .data$source_RNA_class_final,
    
    target_gene =
      .data$target_gene_final,
    
    current_locus_tag =
      .data$current_locus_tag,
    
    protein_id =
      .data$protein_id,
    
    target_product =
      .data$target_product_final,
    
    functional_system =
      .data$functional_system,
    
    functional_annotation_basis =
      .data$functional_annotation_basis,
    
    GO_terms =
      .data$go_terms,
    
    KEGG_orthology =
      .data$kegg_orthology,
    
    KEGG_pathways =
      .data$kegg_pathways,
    
    COG_category =
      .data$cog_category,
    
    COG_description =
      .data$cog_description,
    
    regulon =
      .data$regulon_annotation,
    
    resistance_evidence =
      .data$resistance_evidence,
    
    resistance_annotation =
      .data$resistance_annotation,
    
    virulence_evidence =
      .data$virulence_evidence,
    
    virulence_annotation =
      .data$virulence_annotation,
    
    CLASH_support =
      .data$clash_support_category,
    
    CLASH_evidence_score =
      .data$clash_evidence_score,
    
    supporting_rows =
      .data$number_of_supporting_rows,
    
    hybrid_count =
      .data$total_hybrid_count,
    
    experiments =
      .data$maximum_number_of_experiments,
    
    CLASH_adjusted_P =
      .data$best_adjusted_p_value,
    
    connection_score =
      .data$best_connection_score,
    
    target_transcriptomic_status =
      .data$target_transcriptomic_status,
    
    target_response_direction =
      .data$target_response_direction,
    
    target_log2FC =
      .data$target_log2_fold_change,
    
    target_adjusted_P =
      .data$target_deseq_adjusted_p,
    
    identifier_match_status =
      .data$identifier_match_status,
    
    benchmark_annotation =
      .data$benchmark_annotation,
    
    evidence_profile =
      .data$evidence_profile
  )


# ============================================================
# 20. FUNCTIONAL SYSTEM SUMMARY
# ============================================================

functional_summary <- final_network |>
  
  dplyr::group_by(
    .data$functional_system
  ) |>
  
  dplyr::summarise(
    
    interactions =
      dplyr::n(),
    
    regulatory_RNAs =
      dplyr::n_distinct(
        .data$source_node
      ),
    
    unique_targets =
      dplyr::n_distinct(
        .data$target_node
      ),
    
    significant_target_DEG =
      sum(
        .data$target_is_significant_deg,
        na.rm = TRUE
      ),
    
    strong_target_DEG =
      sum(
        .data$target_strong_response,
        na.rm = TRUE
      ),
    
    higher_CLASH =
      sum(
        .data$clash_support_category ==
          "Higher CLASH support",
        na.rm = TRUE
      ),
    
    intermediate_CLASH =
      sum(
        .data$clash_support_category ==
          "Intermediate CLASH support",
        na.rm = TRUE
      ),
    
    limited_CLASH =
      sum(
        .data$clash_support_category ==
          "Limited CLASH support",
        na.rm = TRUE
      ),
    
    minimal_CLASH =
      sum(
        .data$clash_support_category ==
          "Minimal CLASH support",
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  ) |>
  
  dplyr::arrange(
    dplyr::desc(
      .data$interactions
    )
  )


# ============================================================
# 21. DATABASE ANNOTATION SUMMARY
# ============================================================

database_summary <- tibble::tibble(
  
  annotation_layer =
    c(
      "Identifier matched to master KB",
      "GO annotation",
      "KEGG annotation",
      "COG annotation",
      "Regulon / transcription factor",
      "Explicit resistance annotation",
      "Explicit virulence annotation"
    ),
  
  interactions =
    c(
      sum(
        final_network$KB_match,
        na.rm = TRUE
      ),
      sum(
        final_network$has_GO,
        na.rm = TRUE
      ),
      sum(
        final_network$has_KEGG,
        na.rm = TRUE
      ),
      sum(
        final_network$has_COG,
        na.rm = TRUE
      ),
      sum(
        final_network$has_regulon,
        na.rm = TRUE
      ),
      sum(
        final_network$has_resistance_annotation,
        na.rm = TRUE
      ),
      sum(
        final_network$has_virulence_annotation,
        na.rm = TRUE
      )
    )
) |>
  
  dplyr::mutate(
    
    percent_of_211 =
      round(
        100 *
          .data$interactions /
          nrow(final_network),
        2
      )
  )


# ============================================================
# 22. UNIQUE TARGET DATABASE SUMMARY
# ============================================================

unique_target_annotation <- final_network |>
  
  dplyr::group_by(
    .data$target_node
  ) |>
  
  dplyr::summarise(
    
    target_gene =
      dplyr::first(
        .data$target_gene_final
      ),
    
    product =
      dplyr::first(
        .data$target_product_final
      ),
    
    KB_match =
      any(
        .data$KB_match
      ),
    
    GO =
      any(
        .data$has_GO
      ),
    
    KEGG =
      any(
        .data$has_KEGG
      ),
    
    COG =
      any(
        .data$has_COG
      ),
    
    regulon =
      any(
        .data$has_regulon
      ),
    
    resistance =
      any(
        .data$has_resistance_annotation
      ),
    
    virulence =
      any(
        .data$has_virulence_annotation
      ),
    
    functional_system =
      dplyr::first(
        .data$functional_system
      ),
    
    .groups =
      "drop"
  )


# ============================================================
# 23. RESPONSIVE-TARGET SUBSET
# ============================================================

responsive_subset <- final_network |>
  
  dplyr::filter(
    .data$target_is_significant_deg
  )


# ============================================================
# 24. CRITICAL BENCHMARK / nrdF SUBSET
# ============================================================

critical_subset <- final_network |>
  
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
            .data$target_gene_final,
            .data$target_label,
            .data$target_node
          )
        ),
        
        "\\bnrdf\\b"
      )
  ) |>
  
  dplyr::select(
    
    "source_label",
    "source_RNA_class_final",
    "target_gene_final",
    "current_locus_tag",
    "protein_id",
    "target_product_final",
    "identifier_match_status",
    "functional_system",
    "functional_annotation_basis",
    "go_terms",
    "kegg_orthology",
    "kegg_pathways",
    "cog_category",
    "cog_description",
    "regulon_annotation",
    "resistance_evidence",
    "virulence_evidence",
    "clash_evidence_score",
    "clash_support_category",
    "total_hybrid_count",
    "maximum_number_of_experiments",
    "best_connection_score",
    "best_adjusted_p_value",
    "target_transcriptomic_status",
    "target_log2_fold_change",
    "target_deseq_adjusted_p",
    "benchmark_annotation",
    "evidence_profile"
  )


# ============================================================
# 25. nrdF SUBSET
# ============================================================

nrdf_subset <- final_network |>
  
  dplyr::filter(
    
    stringr::str_detect(
      
      stringr::str_to_lower(
        paste(
          .data$target_gene_final,
          .data$target_label
        )
      ),
      
      "\\bnrdf\\b"
    )
  ) |>
  
  dplyr::select(
    
    "source_label",
    "source_RNA_class_final",
    "target_gene_final",
    "current_locus_tag",
    "protein_id",
    "target_product_final",
    "functional_system",
    "functional_annotation_basis",
    "go_terms",
    "kegg_orthology",
    "kegg_pathways",
    "cog_category",
    "cog_description",
    "clash_evidence_score",
    "clash_support_category",
    "number_of_supporting_rows",
    "total_hybrid_count",
    "maximum_number_of_experiments",
    "best_connection_score",
    "best_adjusted_p_value",
    "target_log2_fold_change",
    "target_deseq_adjusted_p",
    "target_transcriptomic_status",
    "evidence_profile"
  )


# ============================================================
# 26. UNRESOLVED IDENTIFIERS
# ============================================================

unresolved_targets <- final_network |>
  
  dplyr::filter(
    !.data$KB_match
  ) |>
  
  dplyr::select(
    
    "source_label",
    "target_node",
    "target_label",
    "resolved_gene_symbol",
    "identifier_match_status",
    "identifier_match_method",
    "clash_evidence_score",
    "clash_support_category",
    "total_hybrid_count"
  )


# ============================================================
# 27. FIGURE DATA:
# UNIQUE TARGETS BY FUNCTIONAL SYSTEM
# ============================================================

figA_data <- unique_target_annotation |>
  
  dplyr::count(
    .data$functional_system,
    name = "unique_targets"
  ) |>
  
  dplyr::arrange(
    .data$unique_targets
  ) |>
  
  dplyr::mutate(
    
    functional_system =
      factor(
        .data$functional_system,
        levels =
          .data$functional_system
      )
  )


# ============================================================
# 28. FIGURE 3A
# ============================================================

figA <- ggplot2::ggplot(
  
  figA_data,
  
  ggplot2::aes(
    x = .data$unique_targets,
    y = .data$functional_system
  )
) +
  
  ggplot2::geom_col(
    width = 0.7
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label = .data$unique_targets
    ),
    
    hjust = -0.2,
    size = 3.2
  ) +
  
  ggplot2::scale_x_continuous(
    
    expand =
      ggplot2::expansion(
        mult = c(
          0,
          0.15
        )
      )
  ) +
  
  ggplot2::labs(
    
    title =
      "Functional context of CLASH-linked mRNA targets",
    
    subtitle =
      paste(
        "Broad categories summarize master-KB product,",
        "GO, KEGG and COG annotations"
      ),
    
    x =
      "Unique mRNA targets",
    
    y =
      "Functional system"
  ) +
  
  ggplot2::theme_classic(
    base_size = 10
  ) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    axis.title =
      ggplot2::element_text(
        face = "bold"
      )
  )


# ============================================================
# 29. FIGURE 3B:
# ALL NETWORK VS RESPONSIVE TARGETS
# ============================================================

all_unique <- final_network |>
  
  dplyr::distinct(
    .data$target_node,
    .data$functional_system
  ) |>
  
  dplyr::count(
    .data$functional_system,
    name = "n"
  ) |>
  
  dplyr::mutate(
    scope = "All CLASH-linked targets"
  )


responsive_unique <- responsive_subset |>
  
  dplyr::distinct(
    .data$target_node,
    .data$functional_system
  ) |>
  
  dplyr::count(
    .data$functional_system,
    name = "n"
  ) |>
  
  dplyr::mutate(
    scope = "Significant target-DEG subset"
  )


figB_data <- dplyr::bind_rows(
  all_unique,
  responsive_unique
) |>
  
  dplyr::group_by(
    .data$scope
  ) |>
  
  dplyr::mutate(
    
    proportion =
      .data$n /
      sum(
        .data$n
      )
  ) |>
  
  dplyr::ungroup()


figB <- ggplot2::ggplot(
  
  figB_data,
  
  ggplot2::aes(
    x = .data$functional_system,
    y = .data$proportion,
    fill = .data$scope
  )
) +
  
  ggplot2::geom_col(
    position = "dodge"
  ) +
  
  ggplot2::scale_y_continuous(
    labels =
      scales::percent_format(
        accuracy = 1
      )
  ) +
  
  ggplot2::labs(
    
    title =
      "Functional context across independent evidence layers",
    
    subtitle =
      paste(
        "Target differential expression is shown as",
        "a subset and was not required for CLASH inclusion"
      ),
    
    x =
      "Functional system",
    
    y =
      "Proportion of unique targets",
    
    fill =
      NULL
  ) +
  
  ggplot2::theme_classic(
    base_size = 10
  ) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    axis.text.x =
      ggplot2::element_text(
        angle = 40,
        hjust = 1
      ),
    
    legend.position =
      "top"
  )


# ============================================================
# 30. SAVE FIGURES
# ============================================================

print(figA)
print(figB)


ggplot2::ggsave(
  file.path(
    figure_folder,
    "Figure3A_Functional_Context_Unique_Targets.pdf"
  ),
  figA,
  width = 8,
  height = 5.8,
  device = grDevices::cairo_pdf,
  bg = "white"
)


ggplot2::ggsave(
  file.path(
    figure_folder,
    "Figure3A_Functional_Context_Unique_Targets_600dpi.png"
  ),
  figA,
  width = 8,
  height = 5.8,
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  file.path(
    figure_folder,
    "Figure3B_All_vs_DEG_Functional_Context.pdf"
  ),
  figB,
  width = 10,
  height = 6,
  device = grDevices::cairo_pdf,
  bg = "white"
)


ggplot2::ggsave(
  file.path(
    figure_folder,
    "Figure3B_All_vs_DEG_Functional_Context_600dpi.png"
  ),
  figB,
  width = 10,
  height = 6,
  dpi = 600,
  bg = "white"
)


# ============================================================
# 31. WRITE TABLES
# ============================================================

readr::write_csv(
  final_network,
  file.path(
    table_folder,
    "FINAL_Evidence_Layered_Functionally_Annotated_Network.csv"
  )
)


readr::write_csv(
  manuscript_table,
  file.path(
    table_folder,
    "FINAL_Manuscript_Interaction_Annotation_Table.csv"
  )
)


readr::write_csv(
  functional_summary,
  file.path(
    table_folder,
    "FINAL_Functional_System_Summary.csv"
  )
)


readr::write_csv(
  database_summary,
  file.path(
    table_folder,
    "FINAL_Database_Annotation_Summary.csv"
  )
)


readr::write_csv(
  unique_target_annotation,
  file.path(
    table_folder,
    "FINAL_Unique_Target_Annotation_Table.csv"
  )
)


readr::write_csv(
  responsive_subset,
  file.path(
    table_folder,
    "FINAL_Responsive_Target_Subset.csv"
  )
)


readr::write_csv(
  critical_subset,
  file.path(
    table_folder,
    "FINAL_Benchmark_nrdF_Context.csv"
  )
)


readr::write_csv(
  nrdf_subset,
  file.path(
    table_folder,
    "FINAL_nrdF_Candidate_Annotation.csv"
  )
)


readr::write_csv(
  unresolved_targets,
  file.path(
    audit_folder,
    "FINAL_Unresolved_Target_Audit.csv"
  )
)


# ============================================================
# 32. SANITY COUNTS
# ============================================================

total_edges <-
  nrow(
    final_network
  )


unique_rnas <-
  dplyr::n_distinct(
    final_network$source_node
  )


unique_targets <-
  dplyr::n_distinct(
    final_network$target_node
  )


matched_edges <-
  sum(
    final_network$KB_match,
    na.rm = TRUE
  )


matched_unique_targets <-
  unique_target_annotation |>
  dplyr::filter(
    .data$KB_match
  ) |>
  nrow()


responsive_edges <-
  sum(
    final_network$target_is_significant_deg,
    na.rm = TRUE
  )


responsive_targets <-
  final_network |>
  dplyr::filter(
    .data$target_is_significant_deg
  ) |>
  dplyr::summarise(
    n =
      dplyr::n_distinct(
        .data$target_node
      )
  ) |>
  dplyr::pull(
    .data$n
  )


nrdf_edges <-
  nrow(
    nrdf_subset
  )


critical_edges <-
  nrow(
    critical_subset
  )


# ============================================================
# 33. FINAL SUMMARY TABLE
# ============================================================

final_summary <- tibble::tibble(
  
  metric =
    c(
      "Total CLASH-supported interactions",
      "Unique regulatory RNAs",
      "Unique mRNA targets",
      "Master-KB matched interactions",
      "Master-KB unmatched interactions",
      "Master-KB matched unique targets",
      "Master-KB unmatched unique targets",
      "Significant target-DEG interactions",
      "Unique significant target-DEG genes",
      "nrdF candidate interactions",
      "Benchmark/nrdF context interactions",
      "Interactions with GO annotation",
      "Interactions with KEGG annotation",
      "Interactions with COG annotation",
      "Interactions with regulon annotation",
      "Interactions with explicit resistance annotation",
      "Interactions with explicit virulence annotation"
    ),
  
  value =
    c(
      total_edges,
      unique_rnas,
      unique_targets,
      matched_edges,
      total_edges - matched_edges,
      matched_unique_targets,
      unique_targets - matched_unique_targets,
      responsive_edges,
      responsive_targets,
      nrdf_edges,
      critical_edges,
      sum(
        final_network$has_GO,
        na.rm = TRUE
      ),
      sum(
        final_network$has_KEGG,
        na.rm = TRUE
      ),
      sum(
        final_network$has_COG,
        na.rm = TRUE
      ),
      sum(
        final_network$has_regulon,
        na.rm = TRUE
      ),
      sum(
        final_network$has_resistance_annotation,
        na.rm = TRUE
      ),
      sum(
        final_network$has_virulence_annotation,
        na.rm = TRUE
      )
    )
)


readr::write_csv(
  final_summary,
  file.path(
    audit_folder,
    "FINAL_Script16C_Summary.csv"
  )
)


# ============================================================
# 34. INTERPRETATION NOTES
# ============================================================

interpretation_notes <- c(
  
  "SCRIPT 16C FINAL - INTERPRETATION NOTES",
  
  "",
  
  paste(
    "1. The complete 211-edge CLASH-supported network",
    "is retained."
  ),
  
  paste(
    "2. Master-KB identifier resolution is inherited",
    "from Script 16B."
  ),
  
  paste(
    "3. Functional systems summarize existing product,",
    "GO, KEGG and COG information."
  ),
  
  paste(
    "4. Functional-system labels are contextual",
    "annotations and do not validate RNA regulation."
  ),
  
  paste(
    "5. Resistance annotations are reported only when",
    "explicitly represented in the master KB."
  ),
  
  paste(
    "6. Virulence annotations are reported only when",
    "explicitly represented in the master KB."
  ),
  
  paste(
    "7. Regulon information is an independent contextual",
    "annotation layer."
  ),
  
  paste(
    "8. Target transcriptomic response is independent",
    "from CLASH interaction evidence."
  ),
  
  paste(
    "9. Absence of target differential expression does",
    "not exclude translational or other",
    "post-transcriptional regulation."
  ),
  
  paste(
    "10. No biological-priority score or Priority A/B/C",
    "classification is calculated."
  ),
  
  paste(
    "11. Unresolved targets remain unresolved rather",
    "than being assigned speculative annotations."
  ),
  
  paste(
    "12. Putative 3UTR-associated RNA terminology does",
    "not establish transcript processing or release."
  ),
  
  paste(
    "13. The two nrdF interactions remain candidate",
    "interactions supported by CLASH and target",
    "transcriptomic response, not proof of direct",
    "functional regulation."
  )
)


writeLines(
  interpretation_notes,
  file.path(
    audit_folder,
    "FINAL_Interpretation_Notes.txt"
  )
)


# ============================================================
# 35. OUTPUT VERIFICATION
# ============================================================

expected_outputs <- c(
  
  file.path(
    table_folder,
    "FINAL_Evidence_Layered_Functionally_Annotated_Network.csv"
  ),
  
  file.path(
    table_folder,
    "FINAL_Manuscript_Interaction_Annotation_Table.csv"
  ),
  
  file.path(
    table_folder,
    "FINAL_Functional_System_Summary.csv"
  ),
  
  file.path(
    table_folder,
    "FINAL_Database_Annotation_Summary.csv"
  ),
  
  file.path(
    table_folder,
    "FINAL_Benchmark_nrdF_Context.csv"
  ),
  
  file.path(
    table_folder,
    "FINAL_nrdF_Candidate_Annotation.csv"
  ),
  
  file.path(
    audit_folder,
    "FINAL_Script16C_Summary.csv"
  ),
  
  file.path(
    figure_folder,
    "Figure3A_Functional_Context_Unique_Targets_600dpi.png"
  ),
  
  file.path(
    figure_folder,
    "Figure3B_All_vs_DEG_Functional_Context_600dpi.png"
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
    "FINAL_Output_Verification.csv"
  )
)


# ============================================================
# 36. SESSION INFO
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    audit_folder,
    "sessionInfo.txt"
  )
)


# ============================================================
# 37. SANITY WARNINGS
# ============================================================

if (total_edges != 211) {
  
  warning(
    paste0(
      "Expected 211 CLASH-supported interactions; observed ",
      total_edges,
      "."
    )
  )
}


if (unique_targets != 182) {
  
  warning(
    paste0(
      "Expected 182 unique targets; observed ",
      unique_targets,
      "."
    )
  )
}


if (matched_edges != 193) {
  
  warning(
    paste0(
      "Script 16B previously resolved 193 interactions; ",
      "Script 16C currently sees ",
      matched_edges,
      "."
    )
  )
}


if (responsive_edges != 10) {
  
  warning(
    paste0(
      "Expected 10 significant target-DEG interactions; observed ",
      responsive_edges,
      "."
    )
  )
}


if (responsive_targets != 7) {
  
  warning(
    paste0(
      "Expected 7 unique significant target-DEG genes; observed ",
      responsive_targets,
      "."
    )
  )
}


if (nrdf_edges != 2) {
  
  warning(
    paste0(
      "Expected 2 nrdF candidate interactions; observed ",
      nrdf_edges,
      "."
    )
  )
}


if (critical_edges != 6) {
  
  warning(
    paste0(
      "Expected 6 benchmark/nrdF context interactions; observed ",
      critical_edges,
      "."
    )
  )
}


# ============================================================
# 38. FINAL CONSOLE REPORT
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 16C COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================\n\n"
)


cat(
  "FINAL ANALYSIS SUMMARY\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  final_summary,
  n = Inf
)


cat(
  "\nDATABASE ANNOTATION SUMMARY\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  database_summary,
  n = Inf
)


cat(
  "\nFUNCTIONAL SYSTEM SUMMARY\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  functional_summary,
  n = Inf,
  width = Inf
)


cat(
  "\nCRITICAL BENCHMARK / nrdF AUDIT\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  critical_subset,
  n = Inf,
  width = Inf
)


cat(
  "\nnrdF CANDIDATE ANNOTATION\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  nrdf_subset,
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
  "INTERPRETATION\n"
)

cat(
  "============================================\n"
)

cat(
  paste(
    "\nCLASH defines the interaction network.",
    "\nTarget differential expression is an independent evidence layer.",
    "\nGO/KEGG/COG/regulon annotations provide biological context only.",
    "\nResistance and virulence annotations are not inferred from generic metabolism.",
    "\nNo integrated biological-priority score was calculated.",
    "\nNo Priority A/B/C classification was calculated.",
    "\nUnresolved identifiers remain explicitly unresolved.",
    "\nThe nrdF interactions remain candidates requiring independent functional validation.\n"
  )
)


cat(
  "\nMain output folder:\n",
  output_folder,
  "\n"
)


cat(
  "\n============================================\n"
)

cat(
  "END OF SCRIPT 16C FINAL\n"
)

cat(
  "============================================\n"
)
list.files(
  "D:/Bac-sRNA",
  pattern = "dRNA|Term|GSE158830|bed|wig|bigwig|bw|coverage",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
untar(
  "D:/Bac-sRNA/07_processed_data/downloaded_files/GSE158830/GSE158830_RAW.tar",
  list = TRUE
)