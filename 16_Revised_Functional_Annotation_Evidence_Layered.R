# ============================================================
# SCRIPT 16 REVISED
# Evidence-Layered Functional Annotation
# Vancomycin-responsive sRNA / mRNA CLASH Network
# Staphylococcus aureus JKD6008
#
# IMPORTANT PRINCIPLES
# 1. CLASH defines the interaction network.
# 2. Target transcriptomic response is an independent layer.
# 3. Functional annotation is descriptive/contextual.
# 4. Functional annotation is NOT interaction validation.
# 5. No integrated biological-priority score.
# 6. No Priority A / B / C classification.
# 7. Putative UTR-associated RNA terminology is retained
#    until independent processing evidence is demonstrated.
# ============================================================


# ============================================================
# 01. PROJECT PATHS
# ============================================================

project_folder <- "D:/Bac-sRNA"

script14_folder <- file.path(
  project_folder,
  "14_revised_evidence_layered_network"
)

script14_table_folder <- file.path(
  script14_folder,
  "tables"
)


# ============================================================
# 02. INPUT FILES
# ============================================================

network_file <- file.path(
  script14_table_folder,
  "Evidence_Layered_CLASH_Network_ALL.csv"
)

responsive_file <- file.path(
  script14_table_folder,
  "Vancomycin_Responsive_Target_Subset.csv"
)

nrdf_file <- file.path(
  script14_table_folder,
  "nrdF_Candidate_Interactions.csv"
)


# ============================================================
# 03. MASTER KNOWLEDGEBASE CANDIDATES
# ============================================================

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
    paste0(
      "No master knowledgebase was found.\n",
      "Expected one of:\n",
      paste(
        master_candidates,
        collapse = "\n"
      )
    )
  )
}

master_kb_file <- existing_master_files[1]


# ============================================================
# 04. OUTPUT FOLDERS
# ============================================================

output_folder <- file.path(
  project_folder,
  "16_revised_functional_annotation"
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

dir.create(
  audit_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 05. REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "tidyr",
  "tibble",
  "ggplot2",
  "scales"
)


missing_packages <- setdiff(
  required_packages,
  rownames(
    installed.packages()
  )
)


if (length(missing_packages) > 0) {
  
  install.packages(
    missing_packages,
    dependencies = TRUE
  )
}


library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)


# ============================================================
# 06. HELPER FUNCTIONS
# ============================================================

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
    default
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
    
    take <- (
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
# 07. VERIFY INPUT FILES
# ============================================================

if (!file.exists(network_file)) {
  
  stop(
    paste0(
      "Missing main network file:\n",
      network_file
    )
  )
}


if (!file.exists(responsive_file)) {
  
  stop(
    paste0(
      "Missing responsive-target file:\n",
      responsive_file
    )
  )
}


if (!file.exists(nrdf_file)) {
  
  stop(
    paste0(
      "Missing nrdF candidate file:\n",
      nrdf_file
    )
  )
}


# ============================================================
# 08. READ INPUT DATA
# ============================================================

network <- readr::read_csv(
  network_file,
  show_col_types = FALSE
)


responsive <- readr::read_csv(
  responsive_file,
  show_col_types = FALSE
)


nrdf <- readr::read_csv(
  nrdf_file,
  show_col_types = FALSE
)


master_kb <- readr::read_csv(
  master_kb_file,
  show_col_types = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 16 REVISED\n"
)

cat(
  "============================================\n\n"
)

cat(
  "Main CLASH network rows:",
  nrow(network),
  "\n"
)

cat(
  "Responsive-target subset rows:",
  nrow(responsive),
  "\n"
)

cat(
  "nrdF candidate rows:",
  nrow(nrdf),
  "\n"
)

cat(
  "Master knowledgebase:\n",
  master_kb_file,
  "\n\n"
)


# ============================================================
# 09. ENSURE REQUIRED NETWORK COLUMNS EXIST
# ============================================================

network_defaults <- list(
  
  source_node =
    NA_character_,
  
  source_label =
    NA_character_,
  
  source_raw_names =
    NA_character_,
  
  source_RNA_class =
    NA_character_,
  
  target_node =
    NA_character_,
  
  target_label =
    NA_character_,
  
  target_raw_names =
    NA_character_,
  
  number_of_supporting_rows =
    NA_real_,
  
  total_hybrid_count =
    NA_real_,
  
  maximum_number_of_experiments =
    NA_real_,
  
  best_adjusted_p_value =
    NA_real_,
  
  best_connection_score =
    NA_real_,
  
  clash_evidence_score =
    NA_real_,
  
  clash_support_category =
    NA_character_,
  
  target_log2_fold_change =
    NA_real_,
  
  target_deseq_adjusted_p =
    NA_real_,
  
  target_is_significant_deg =
    FALSE,
  
  target_strong_response =
    FALSE,
  
  target_transcriptomic_layer =
    NA_character_,
  
  benchmark_annotation =
    NA_character_,
  
  evidence_profile =
    NA_character_
)


for (nm in names(network_defaults)) {
  
  network <- add_missing(
    network,
    nm,
    network_defaults[[nm]]
  )
}


# ============================================================
# 10. STANDARDIZE NETWORK COLUMN TYPES
# ============================================================

network <- network |>
  
  dplyr::mutate(
    
    source_node =
      safe_chr(
        .data$source_node
      ),
    
    source_label =
      safe_chr(
        .data$source_label
      ),
    
    source_raw_names =
      safe_chr(
        .data$source_raw_names
      ),
    
    source_RNA_class =
      safe_chr(
        .data$source_RNA_class
      ),
    
    target_node =
      safe_chr(
        .data$target_node
      ),
    
    target_label =
      safe_chr(
        .data$target_label
      ),
    
    target_raw_names =
      safe_chr(
        .data$target_raw_names
      ),
    
    clash_support_category =
      safe_chr(
        .data$clash_support_category
      ),
    
    target_transcriptomic_layer =
      safe_chr(
        .data$target_transcriptomic_layer
      ),
    
    benchmark_annotation =
      safe_chr(
        .data$benchmark_annotation
      ),
    
    evidence_profile =
      safe_chr(
        .data$evidence_profile
      ),
    
    number_of_supporting_rows =
      safe_num(
        .data$number_of_supporting_rows
      ),
    
    total_hybrid_count =
      safe_num(
        .data$total_hybrid_count
      ),
    
    maximum_number_of_experiments =
      safe_num(
        .data$maximum_number_of_experiments
      ),
    
    best_adjusted_p_value =
      safe_num(
        .data$best_adjusted_p_value
      ),
    
    best_connection_score =
      safe_num(
        .data$best_connection_score
      ),
    
    clash_evidence_score =
      safe_num(
        .data$clash_evidence_score
      ),
    
    target_log2_fold_change =
      safe_num(
        .data$target_log2_fold_change
      ),
    
    target_deseq_adjusted_p =
      safe_num(
        .data$target_deseq_adjusted_p
      ),
    
    target_is_significant_deg =
      safe_logical(
        .data$target_is_significant_deg
      ),
    
    target_strong_response =
      safe_logical(
        .data$target_strong_response
      )
  )


# ============================================================
# 11. CONSERVATIVE RNA CLASSIFICATION
# ============================================================

network <- network |>
  
  dplyr::mutate(
    
    source_RNA_class_revised =
      dplyr::case_when(
        
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
            "Intergenic",
            ignore_case = TRUE
          )
        ) ~
          "Intergenic RNA candidate",
        
        stringr::str_detect(
          .data$source_RNA_class,
          stringr::regex(
            "Numbered",
            ignore_case = TRUE
          )
        ) ~
          "Numbered sRNA candidate",
        
        TRUE ~
          "Other regulatory RNA candidate"
      )
  )


# ============================================================
# 12. PREPARE MASTER KNOWLEDGEBASE
# ============================================================

master_defaults <- list(
  
  gene_id =
    NA_character_,
  
  current_locus_tag =
    NA_character_,
  
  gene_symbol =
    NA_character_,
  
  gene_synonym =
    NA_character_,
  
  master_label =
    NA_character_,
  
  product =
    NA_character_,
  
  product_short =
    NA_character_,
  
  resistance_gene =
    NA_character_,
  
  resistance_mechanism =
    NA_character_,
  
  virulence_factor =
    NA_character_,
  
  regulon =
    NA_character_,
  
  go_terms =
    NA_character_,
  
  pathway =
    NA_character_,
  
  functional_category =
    NA_character_
)


for (nm in names(master_defaults)) {
  
  master_kb <- add_missing(
    master_kb,
    nm,
    master_defaults[[nm]]
  )
}


# ============================================================
# 13. BUILD MASTER ANNOTATION TABLE
# ============================================================

master_annotation <- master_kb |>
  
  dplyr::transmute(
    
    target_node =
      safe_chr(
        .data$gene_id
      ),
    
    current_locus_tag =
      safe_chr(
        .data$current_locus_tag
      ),
    
    gene_symbol_master =
      safe_chr(
        .data$gene_symbol
      ),
    
    gene_synonym_master =
      safe_chr(
        .data$gene_synonym
      ),
    
    master_label =
      safe_chr(
        .data$master_label
      ),
    
    product_master =
      first_nonblank(
        .data$product_short,
        .data$product
      ),
    
    existing_resistance_gene =
      safe_chr(
        .data$resistance_gene
      ),
    
    existing_resistance_mechanism =
      safe_chr(
        .data$resistance_mechanism
      ),
    
    existing_virulence_factor =
      safe_chr(
        .data$virulence_factor
      ),
    
    existing_regulon =
      safe_chr(
        .data$regulon
      ),
    
    existing_go_terms =
      safe_chr(
        .data$go_terms
      ),
    
    existing_pathway =
      safe_chr(
        .data$pathway
      ),
    
    existing_functional_category =
      safe_chr(
        .data$functional_category
      )
  ) |>
  
  dplyr::filter(
    .data$target_node != ""
  ) |>
  
  dplyr::distinct(
    .data$target_node,
    .keep_all = TRUE
  )


# ============================================================
# 14. MANUAL DESCRIPTIVE DICTIONARY
#
# IMPORTANT:
# These descriptions provide functional context only.
# They are NOT used to score or rank interactions.
# ============================================================

manual_dictionary <- tibble::tribble(
  
  ~gene_symbol,
  ~manual_primary_function,
  ~manual_functional_system,
  ~manual_note,
  
  "spa",
  "Immunoglobulin-binding surface protein A",
  "Virulence and host interaction",
  paste(
    "Known surface protein.",
    "Functional context is distinct from evidence",
    "for the specific RNA interaction."
  ),
  
  "guaB",
  "Inosine monophosphate dehydrogenase",
  "Nucleotide metabolism",
  paste(
    "Purine-biosynthesis enzyme.",
    "No direct resistance role is assigned here."
  ),
  
  "nrdF",
  "Ribonucleotide-diphosphate reductase beta subunit",
  "DNA synthesis and repair",
  paste(
    "Deoxyribonucleotide synthesis.",
    "The RNA-nrdF interactions remain candidates."
  ),
  
  "gltB",
  "Glutamate synthase large subunit",
  "Amino-acid and nitrogen metabolism",
  paste(
    "Metabolic annotation only.",
    "No direct resistance role is assigned here."
  ),
  
  "modA",
  "Molybdate ABC transporter substrate-binding protein",
  "Transport and membrane processes",
  paste(
    "Transport annotation only.",
    "No direct resistance role is assigned here."
  ),
  
  "dhaK",
  "Dihydroxyacetone kinase subunit",
  "Carbon and central metabolism",
  paste(
    "Carbon-metabolism annotation only.",
    "No direct resistance role is assigned here."
  ),
  
  "qoxB",
  "Quinol oxidase subunit",
  "Respiration and redox metabolism",
  paste(
    "Respiratory annotation only.",
    "No direct resistance role is assigned here."
  )
) |>
  
  dplyr::mutate(
    
    gene_symbol_norm =
      stringr::str_to_lower(
        .data$gene_symbol
      )
  )


# ============================================================
# 15. JOIN NETWORK TO MASTER KNOWLEDGEBASE
# ============================================================

annotated_network <- network |>
  
  dplyr::left_join(
    master_annotation,
    by = "target_node"
  ) |>
  
  dplyr::mutate(
    
    final_gene_symbol =
      first_nonblank(
        
        .data$gene_symbol_master,
        
        .data$target_label,
        
        .data$master_label,
        
        .data$target_raw_names,
        
        .data$target_node
      ),
    
    final_product =
      first_nonblank(
        
        .data$product_master,
        
        .data$target_label,
        
        rep(
          "Uncharacterized product",
          dplyr::n()
        )
      ),
    
    gene_symbol_norm =
      stringr::str_to_lower(
        .data$final_gene_symbol
      )
  ) |>
  
  dplyr::left_join(
    manual_dictionary,
    by = "gene_symbol_norm"
  )


# ============================================================
# 16. BUILD TEXT FOR FUNCTIONAL INFERENCE
# ============================================================

annotated_network <- annotated_network |>
  
  dplyr::mutate(
    
    annotation_text =
      stringr::str_to_lower(
        
        paste(
          
          .data$final_gene_symbol,
          
          .data$final_product,
          
          .data$existing_functional_category,
          
          .data$existing_pathway,
          
          .data$existing_go_terms,
          
          sep = " "
        )
      )
  )


# ============================================================
# 17. KEYWORD-BASED FUNCTIONAL SYSTEM INFERENCE
#
# This is a descriptive fallback only.
# ============================================================

annotated_network <- annotated_network |>
  
  dplyr::mutate(
    
    inferred_functional_system =
      dplyr::case_when(
        
        stringr::str_detect(
          .data$annotation_text,
          paste0(
            "cell wall|",
            "peptidoglycan|",
            "penicillin.?binding|",
            "\\bmur[a-z]\\b|",
            "\\bpbp|",
            "teichoic"
          )
        ) ~
          "Cell wall and envelope",
        
        stringr::str_detect(
          .data$annotation_text,
          paste0(
            "transporter|",
            "transport|",
            "permease|",
            "\\babc\\b|",
            "efflux|",
            "uptake"
          )
        ) ~
          "Transport and membrane processes",
        
        stringr::str_detect(
          .data$annotation_text,
          paste0(
            "virulence|",
            "adhesin|",
            "surface protein|",
            "toxin|",
            "hemolysin|",
            "protein a"
          )
        ) ~
          "Virulence and host interaction",
        
        stringr::str_detect(
          .data$annotation_text,
          paste0(
            "ribonucleotide|",
            "dna |",
            "replication|",
            "repair|",
            "polymerase"
          )
        ) ~
          "DNA synthesis and repair",
        
        stringr::str_detect(
          .data$annotation_text,
          paste0(
            "ribosom|",
            "translation|",
            "trna|",
            "aminoacyl"
          )
        ) ~
          "Translation and protein synthesis",
        
        stringr::str_detect(
          .data$annotation_text,
          paste0(
            "purine|",
            "pyrimidine|",
            "nucleotide|",
            "guanine|",
            "adenine"
          )
        ) ~
          "Nucleotide metabolism",
        
        stringr::str_detect(
          .data$annotation_text,
          paste0(
            "glutamate|",
            "amino acid|",
            "nitrogen|",
            "arginine|",
            "histidine|",
            "leucine|",
            "glycine"
          )
        ) ~
          "Amino-acid and nitrogen metabolism",
        
        stringr::str_detect(
          .data$annotation_text,
          paste0(
            "glycol|",
            "carbon|",
            "sugar|",
            "dihydroxyacetone|",
            "dehydrogenase"
          )
        ) ~
          "Carbon and central metabolism",
        
        stringr::str_detect(
          .data$annotation_text,
          paste0(
            "oxidase|",
            "respiration|",
            "electron transport|",
            "cytochrome|",
            "quinol"
          )
        ) ~
          "Respiration and redox metabolism",
        
        stringr::str_detect(
          .data$annotation_text,
          paste0(
            "stress|",
            "chaperone|",
            "heat shock|",
            "oxidative|",
            "detox"
          )
        ) ~
          "Cellular stress response",
        
        TRUE ~
          "Other or unclassified function"
      )
  )


# ============================================================
# 18. FINAL FUNCTIONAL SYSTEM
# ============================================================

annotated_network <- annotated_network |>
  
  dplyr::mutate(
    
    final_functional_system =
      dplyr::case_when(
        
        !is.na(
          .data$manual_functional_system
        ) &
          .data$manual_functional_system != "" ~
          
          .data$manual_functional_system,
        
        .data$existing_functional_category != "" ~
          
          .data$existing_functional_category,
        
        TRUE ~
          
          .data$inferred_functional_system
      ),
    
    final_primary_function =
      dplyr::case_when(
        
        !is.na(
          .data$manual_primary_function
        ) &
          .data$manual_primary_function != "" ~
          
          .data$manual_primary_function,
        
        .data$final_product != "" ~
          
          .data$final_product,
        
        TRUE ~
          
          "Uncharacterized function"
      )
  )


# ============================================================
# 19. FUNCTIONAL ANNOTATION PROVENANCE
# ============================================================

annotated_network <- annotated_network |>
  
  dplyr::mutate(
    
    functional_annotation_source =
      dplyr::case_when(
        
        !is.na(
          .data$manual_primary_function
        ) &
          .data$manual_primary_function != "" ~
          
          "Manual descriptive annotation",
        
        .data$existing_functional_category != "" |
          .data$existing_pathway != "" |
          .data$existing_go_terms != "" ~
          
          "Existing master annotation",
        
        .data$inferred_functional_system !=
          "Other or unclassified function" ~
          
          "Keyword-based functional inference",
        
        TRUE ~
          
          "Unresolved"
      ),
    
    annotation_evidence_level =
      dplyr::case_when(
        
        .data$existing_functional_category != "" |
          .data$existing_pathway != "" |
          .data$existing_go_terms != "" ~
          
          "Master/database-derived annotation",
        
        !is.na(
          .data$manual_primary_function
        ) &
          .data$manual_primary_function != "" ~
          
          "Manual descriptive annotation",
        
        .data$inferred_functional_system !=
          "Other or unclassified function" ~
          
          "Keyword-based inference",
        
        TRUE ~
          
          "Unresolved"
      )
  )


# ============================================================
# 20. RESISTANCE CONTEXT
#
# IMPORTANT:
# Explicit annotations are kept separate from inferred context.
# ============================================================

annotated_network <- annotated_network |>
  
  dplyr::mutate(
    
    resistance_master_supported =
      
      .data$existing_resistance_gene != "" |
      
      .data$existing_resistance_mechanism != "",
    
    resistance_context =
      dplyr::case_when(
        
        .data$resistance_master_supported ~
          
          first_nonblank(
            
            .data$existing_resistance_mechanism,
            
            .data$existing_resistance_gene
          ),
        
        .data$final_functional_system ==
          "Cell wall and envelope" ~
          
          paste(
            "Functional context potentially relevant",
            "to glycopeptide response;",
            "not a resistance annotation"
          ),
        
        .data$final_functional_system %in%
          c(
            "Transport and membrane processes",
            "Respiration and redox metabolism",
            "Cellular stress response"
          ) ~
          
          paste(
            "Stress/adaptation-related functional context;",
            "not a direct resistance annotation"
          ),
        
        TRUE ~
          
          "No direct resistance annotation assigned"
      ),
    
    resistance_annotation_source =
      dplyr::case_when(
        
        .data$resistance_master_supported ~
          
          "Existing master annotation",
        
        .data$resistance_context !=
          "No direct resistance annotation assigned" ~
          
          "Functional-context inference",
        
        TRUE ~
          
          "None"
      )
  )


# ============================================================
# 21. VIRULENCE CONTEXT
# ============================================================

annotated_network <- annotated_network |>
  
  dplyr::mutate(
    
    virulence_master_supported =
      
      .data$existing_virulence_factor != "",
    
    virulence_context =
      dplyr::case_when(
        
        .data$virulence_master_supported ~
          
          .data$existing_virulence_factor,
        
        .data$final_functional_system ==
          "Virulence and host interaction" ~
          
          paste(
            "Virulence-related functional context;",
            "interaction-specific effect not established"
          ),
        
        TRUE ~
          
          "No direct virulence annotation assigned"
      ),
    
    virulence_annotation_source =
      dplyr::case_when(
        
        .data$virulence_master_supported ~
          
          "Existing master annotation",
        
        .data$final_functional_system ==
          "Virulence and host interaction" ~
          
          "Functional-context inference",
        
        TRUE ~
          
          "None"
      )
  )


# ============================================================
# 22. REGULON CONTEXT
# ============================================================

annotated_network <- annotated_network |>
  
  dplyr::mutate(
    
    regulon_context =
      dplyr::case_when(
        
        .data$existing_regulon != "" ~
          
          .data$existing_regulon,
        
        TRUE ~
          
          "No mapped regulon annotation"
      ),
    
    regulon_annotation_source =
      dplyr::if_else(
        
        .data$existing_regulon != "",
        
        "Existing master annotation",
        
        "None"
      )
  )


# ============================================================
# 23. TARGET TRANSCRIPTOMIC LAYER
# ============================================================

annotated_network <- annotated_network |>
  
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
          
          "Measured, below DEG threshold",
        
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
# 24. MANUSCRIPT-SAFE FUNCTIONAL TABLE
# ============================================================

manuscript_table <- annotated_network |>
  
  dplyr::transmute(
    
    regulatory_RNA =
      .data$source_label,
    
    RNA_class =
      .data$source_RNA_class_revised,
    
    target_gene =
      .data$final_gene_symbol,
    
    target_locus =
      .data$target_node,
    
    annotated_function =
      .data$final_primary_function,
    
    functional_system =
      .data$final_functional_system,
    
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
    
    target_log2_fold_change =
      .data$target_log2_fold_change,
    
    target_adjusted_P =
      .data$target_deseq_adjusted_p,
    
    resistance_context =
      .data$resistance_context,
    
    resistance_annotation_source =
      .data$resistance_annotation_source,
    
    virulence_context =
      .data$virulence_context,
    
    virulence_annotation_source =
      .data$virulence_annotation_source,
    
    regulon_context =
      .data$regulon_context,
    
    annotation_evidence_level =
      .data$annotation_evidence_level,
    
    functional_annotation_source =
      .data$functional_annotation_source,
    
    benchmark_annotation =
      .data$benchmark_annotation,
    
    evidence_profile =
      .data$evidence_profile,
    
    interpretation_caution =
      paste(
        "CLASH supports RNA-RNA contact;",
        "functional direction and causal regulation",
        "require independent validation."
      )
  ) |>
  
  dplyr::arrange(
    
    dplyr::desc(
      .data$CLASH_evidence_score
    ),
    
    dplyr::desc(
      .data$hybrid_count
    )
  )


# ============================================================
# 25. DEG-ASSOCIATED SUBSET
# ============================================================

responsive_functional_subset <- manuscript_table |>
  
  dplyr::filter(
    
    .data$target_transcriptomic_status %in%
      
      c(
        "Strong target DEG",
        "Significant target DEG"
      )
  )


# ============================================================
# 26. BENCHMARK + nrdF CONTEXT
# ============================================================

critical_subset <- annotated_network |>
  
  dplyr::filter(
    
    stringr::str_detect(
      
      stringr::str_to_lower(
        
        paste(
          
          .data$source_label,
          
          .data$source_raw_names,
          
          .data$source_node
        )
      ),
      
      "rsaoi|sau[-_ ]?6477|spra2|rsaj"
    ) |
      
      stringr::str_detect(
        
        stringr::str_to_lower(
          
          paste(
            
            .data$final_gene_symbol,
            
            .data$target_label,
            
            .data$target_node
          )
        ),
        
        "\\bnrdf\\b"
      )
  ) |>
  
  dplyr::transmute(
    
    regulatory_RNA =
      .data$source_label,
    
    RNA_class =
      .data$source_RNA_class_revised,
    
    target_gene =
      .data$final_gene_symbol,
    
    functional_system =
      .data$final_functional_system,
    
    CLASH_support =
      .data$clash_support_category,
    
    CLASH_evidence_score =
      .data$clash_evidence_score,
    
    hybrid_count =
      .data$total_hybrid_count,
    
    experiments =
      .data$maximum_number_of_experiments,
    
    connection_score =
      .data$best_connection_score,
    
    CLASH_adjusted_P =
      .data$best_adjusted_p_value,
    
    target_transcriptomic_status =
      .data$target_transcriptomic_status,
    
    target_log2_fold_change =
      .data$target_log2_fold_change,
    
    target_adjusted_P =
      .data$target_deseq_adjusted_p,
    
    functional_annotation_source =
      .data$functional_annotation_source,
    
    resistance_context =
      .data$resistance_context,
    
    virulence_context =
      .data$virulence_context,
    
    benchmark_annotation =
      .data$benchmark_annotation,
    
    evidence_profile =
      .data$evidence_profile
  )


# ============================================================
# 27. FUNCTIONAL SYSTEM SUMMARY
# ============================================================

functional_summary <- annotated_network |>
  
  dplyr::group_by(
    .data$final_functional_system
  ) |>
  
  dplyr::summarise(
    
    number_of_interactions =
      dplyr::n(),
    
    number_of_regulatory_RNAs =
      dplyr::n_distinct(
        .data$source_node
      ),
    
    number_of_target_genes =
      dplyr::n_distinct(
        .data$target_node
      ),
    
    significant_DEG_interactions =
      sum(
        .data$target_is_significant_deg,
        na.rm = TRUE
      ),
    
    strong_DEG_interactions =
      sum(
        .data$target_strong_response,
        na.rm = TRUE
      ),
    
    higher_CLASH_support =
      sum(
        .data$clash_support_category ==
          "Higher CLASH support",
        na.rm = TRUE
      ),
    
    intermediate_CLASH_support =
      sum(
        .data$clash_support_category ==
          "Intermediate CLASH support",
        na.rm = TRUE
      ),
    
    limited_CLASH_support =
      sum(
        .data$clash_support_category ==
          "Limited CLASH support",
        na.rm = TRUE
      ),
    
    minimal_CLASH_support =
      sum(
        .data$clash_support_category ==
          "Minimal CLASH support",
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  
  dplyr::arrange(
    
    dplyr::desc(
      .data$number_of_interactions
    )
  )


# ============================================================
# 28. ANNOTATION SOURCE SUMMARY
# ============================================================

annotation_source_summary <- annotated_network |>
  
  dplyr::count(
    
    .data$annotation_evidence_level,
    
    name =
      "number_of_interactions"
  ) |>
  
  dplyr::arrange(
    
    dplyr::desc(
      .data$number_of_interactions
    )
  )


# ============================================================
# 29. SCRIPT 16 SUMMARY
# ============================================================

nrdf_count <- sum(
  
  stringr::str_detect(
    
    stringr::str_to_lower(
      
      paste(
        
        annotated_network$final_gene_symbol,
        
        annotated_network$target_label
      )
    ),
    
    "\\bnrdf\\b"
  ),
  
  na.rm = TRUE
)


rsaoi_count <- sum(
  
  stringr::str_detect(
    
    stringr::str_to_lower(
      
      paste(
        
        annotated_network$source_label,
        
        annotated_network$source_raw_names
      )
    ),
    
    "rsaoi|sau[-_ ]?6477"
  ),
  
  na.rm = TRUE
)


spra2_count <- sum(
  
  stringr::str_detect(
    
    stringr::str_to_lower(
      
      paste(
        
        annotated_network$source_label,
        
        annotated_network$source_raw_names
      )
    ),
    
    "spra2|rsaj"
  ),
  
  na.rm = TRUE
)


script16_summary <- tibble::tibble(
  
  metric = c(
    
    "Total CLASH-supported interactions",
    
    "Unique regulatory RNAs",
    
    "Unique mRNA targets",
    
    paste(
      "Higher/intermediate",
      "CLASH-support interactions"
    ),
    
    "Significant target-DEG interactions",
    
    "Unique significant target-DEG genes",
    
    paste(
      "Measured target but below",
      "DEG threshold interactions"
    ),
    
    paste(
      "Interactions without mapped",
      "target transcriptomic evidence"
    ),
    
    "nrdF candidate interactions",
    
    "RsaOI interactions",
    
    "SprA2/RsaJ-associated interactions",
    
    paste(
      "Master/database-derived",
      "functional annotations"
    ),
    
    paste(
      "Manual descriptive",
      "functional annotations"
    ),
    
    paste(
      "Keyword-inferred",
      "functional annotations"
    ),
    
    "Unresolved functional annotations"
  ),
  
  value = c(
    
    nrow(
      annotated_network
    ),
    
    dplyr::n_distinct(
      annotated_network$source_node
    ),
    
    dplyr::n_distinct(
      annotated_network$target_node
    ),
    
    sum(
      
      annotated_network$
        clash_support_category %in%
        
        c(
          "Higher CLASH support",
          "Intermediate CLASH support"
        ),
      
      na.rm = TRUE
    ),
    
    sum(
      
      annotated_network$
        target_is_significant_deg,
      
      na.rm = TRUE
    ),
    
    dplyr::n_distinct(
      
      annotated_network$
        target_node[
          annotated_network$
            target_is_significant_deg
        ]
    ),
    
    sum(
      
      !annotated_network$
        target_is_significant_deg &
        
        !is.na(
          annotated_network$
            target_log2_fold_change
        ),
      
      na.rm = TRUE
    ),
    
    sum(
      
      is.na(
        annotated_network$
          target_log2_fold_change
      ),
      
      na.rm = TRUE
    ),
    
    nrdf_count,
    
    rsaoi_count,
    
    spra2_count,
    
    sum(
      
      annotated_network$
        annotation_evidence_level ==
        "Master/database-derived annotation",
      
      na.rm = TRUE
    ),
    
    sum(
      
      annotated_network$
        annotation_evidence_level ==
        "Manual descriptive annotation",
      
      na.rm = TRUE
    ),
    
    sum(
      
      annotated_network$
        annotation_evidence_level ==
        "Keyword-based inference",
      
      na.rm = TRUE
    ),
    
    sum(
      
      annotated_network$
        annotation_evidence_level ==
        "Unresolved",
      
      na.rm = TRUE
    )
  )
)


# ============================================================
# 30. MISSING / UNCLASSIFIED ANNOTATION AUDIT
# ============================================================

missing_annotation_audit <- annotated_network |>
  
  dplyr::filter(
    
    .data$final_functional_system ==
      "Other or unclassified function" |
      
      .data$annotation_evidence_level ==
      "Unresolved"
  ) |>
  
  dplyr::select(
    
    .data$source_label,
    
    .data$target_node,
    
    .data$final_gene_symbol,
    
    .data$final_product,
    
    .data$final_functional_system,
    
    .data$annotation_evidence_level
  )


# ============================================================
# 31. FIGURE 3A DATA
# ============================================================

functional_plot_data <- functional_summary |>
  
  dplyr::mutate(
    
    display_category =
      stringr::str_wrap(
        .data$final_functional_system,
        width = 34
      ),
    
    display_category =
      stats::reorder(
        
        .data$display_category,
        
        .data$number_of_interactions
      )
  )


# ============================================================
# 32. FIGURE 3A
# Functional composition of complete CLASH network
# ============================================================

functional_plot <- ggplot2::ggplot(
  
  functional_plot_data,
  
  ggplot2::aes(
    
    x =
      .data$number_of_interactions,
    
    y =
      .data$display_category
  )
) +
  
  ggplot2::geom_col(
    width = 0.70
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      
      label =
        .data$number_of_interactions
    ),
    
    hjust = -0.20,
    
    size = 3.2
  ) +
  
  ggplot2::scale_x_continuous(
    
    expand =
      ggplot2::expansion(
        
        mult =
          c(
            0,
            0.18
          )
      )
  ) +
  
  ggplot2::labs(
    
    x =
      paste(
        "Number of CLASH-supported",
        "RNA-mRNA interactions"
      ),
    
    y =
      "Functional system",
    
    title =
      paste(
        "Functional composition of the",
        "CLASH-supported network"
      ),
    
    subtitle =
      paste(
        "Functional categories are descriptive",
        "and do not establish regulatory effect"
      )
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
      ),
    
    plot.margin =
      ggplot2::margin(
        12,
        35,
        12,
        12
      )
  )


# ============================================================
# 33. FIGURE 3B DATA
# ============================================================

evidence_plot_data <- annotated_network |>
  
  dplyr::count(
    
    .data$clash_support_category,
    
    .data$target_transcriptomic_status,
    
    name = "n"
  )


# ============================================================
# 34. FIGURE 3B
# CLASH evidence versus transcriptomic context
# ============================================================

evidence_plot <- ggplot2::ggplot(
  
  evidence_plot_data,
  
  ggplot2::aes(
    
    x =
      .data$clash_support_category,
    
    y =
      .data$n,
    
    fill =
      .data$target_transcriptomic_status
  )
) +
  
  ggplot2::geom_col() +
  
  ggplot2::labs(
    
    x =
      "Descriptive CLASH-support category",
    
    y =
      "Number of interactions",
    
    fill =
      "Target transcriptomic layer",
    
    title =
      paste(
        "CLASH evidence and target",
        "transcriptomic context"
      ),
    
    subtitle =
      paste(
        "Target differential expression is",
        "an orthogonal evidence layer,",
        "not an inclusion criterion"
      )
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
      ),
    
    axis.text.x =
      ggplot2::element_text(
        angle = 25,
        hjust = 1
      ),
    
    legend.position =
      "right"
  )


# ============================================================
# 35. DISPLAY FIGURES
# ============================================================

print(
  functional_plot
)

print(
  evidence_plot
)


# ============================================================
# 36. SAVE FIGURE 3A
# ============================================================

ggplot2::ggsave(
  
  filename =
    file.path(
      figure_folder,
      "Figure3A_Functional_Composition.pdf"
    ),
  
  plot =
    functional_plot,
  
  width =
    8,
  
  height =
    5.5,
  
  device =
    grDevices::cairo_pdf,
  
  bg =
    "white"
)


ggplot2::ggsave(
  
  filename =
    file.path(
      figure_folder,
      "Figure3A_Functional_Composition_600dpi.png"
    ),
  
  plot =
    functional_plot,
  
  width =
    8,
  
  height =
    5.5,
  
  dpi =
    600,
  
  bg =
    "white"
)


ggplot2::ggsave(
  
  filename =
    file.path(
      figure_folder,
      "Figure3A_Functional_Composition_600dpi.tiff"
    ),
  
  plot =
    functional_plot,
  
  width =
    8,
  
  height =
    5.5,
  
  dpi =
    600,
  
  compression =
    "lzw",
  
  bg =
    "white"
)


# ============================================================
# 37. SAVE FIGURE 3B
# ============================================================

ggplot2::ggsave(
  
  filename =
    file.path(
      figure_folder,
      "Figure3B_CLASH_vs_Transcriptomic_Context.pdf"
    ),
  
  plot =
    evidence_plot,
  
  width =
    9,
  
  height =
    5.8,
  
  device =
    grDevices::cairo_pdf,
  
  bg =
    "white"
)


ggplot2::ggsave(
  
  filename =
    file.path(
      figure_folder,
      "Figure3B_CLASH_vs_Transcriptomic_Context_600dpi.png"
    ),
  
  plot =
    evidence_plot,
  
  width =
    9,
  
  height =
    5.8,
  
  dpi =
    600,
  
  bg =
    "white"
)


ggplot2::ggsave(
  
  filename =
    file.path(
      figure_folder,
      "Figure3B_CLASH_vs_Transcriptomic_Context_600dpi.tiff"
    ),
  
  plot =
    evidence_plot,
  
  width =
    9,
  
  height =
    5.8,
  
  dpi =
    600,
  
  compression =
    "lzw",
  
  bg =
    "white"
)


# ============================================================
# 38. SAVE MAIN ANNOTATED NETWORK
# ============================================================

readr::write_csv(
  
  annotated_network,
  
  file.path(
    
    table_folder,
    
    "Evidence_Layered_Network_Functionally_Annotated.csv"
  )
)


# ============================================================
# 39. SAVE MANUSCRIPT TABLE
# ============================================================

readr::write_csv(
  
  manuscript_table,
  
  file.path(
    
    table_folder,
    
    "Manuscript_Functional_Annotation_Table.csv"
  )
)


# ============================================================
# 40. SAVE DEG-ASSOCIATED SUBSET
# ============================================================

readr::write_csv(
  
  responsive_functional_subset,
  
  file.path(
    
    table_folder,
    
    "DEG_Associated_Functional_Subset.csv"
  )
)


# ============================================================
# 41. SAVE BENCHMARK + nrdF TABLE
# ============================================================

readr::write_csv(
  
  critical_subset,
  
  file.path(
    
    table_folder,
    
    "Benchmark_and_nrdF_Functional_Context.csv"
  )
)


# ============================================================
# 42. SAVE FUNCTIONAL SUMMARY
# ============================================================

readr::write_csv(
  
  functional_summary,
  
  file.path(
    
    table_folder,
    
    "Functional_System_Summary.csv"
  )
)


# ============================================================
# 43. SAVE ANNOTATION AUDIT
# ============================================================

readr::write_csv(
  
  annotation_source_summary,
  
  file.path(
    
    audit_folder,
    
    "Annotation_Evidence_Source_Summary.csv"
  )
)


readr::write_csv(
  
  missing_annotation_audit,
  
  file.path(
    
    audit_folder,
    
    "Missing_or_Unclassified_Functional_Annotations.csv"
  )
)


readr::write_csv(
  
  script16_summary,
  
  file.path(
    
    audit_folder,
    
    "Script16_Revised_Summary.csv"
  )
)


# ============================================================
# 44. METHODOLOGICAL NOTES
# ============================================================

methodological_notes <- c(
  
  "Script 16 Revised - methodological notes",
  
  "",
  
  paste(
    "1. The complete CLASH-supported sRNA-mRNA network",
    "is annotated; target DEG is not an inclusion criterion."
  ),
  
  paste(
    "2. CLASH evidence, target transcriptomic response,",
    "and functional annotation are retained as",
    "separate evidence layers."
  ),
  
  paste(
    "3. No integrated biological-priority score",
    "or Priority A/B/C classification is calculated."
  ),
  
  paste(
    "4. Keyword-based functional assignments are",
    "descriptive inference, not database validation",
    "or experimental validation."
  ),
  
  paste(
    "5. Resistance and virulence labels from the",
    "master knowledgebase are kept separate from",
    "functional-context inference."
  ),
  
  paste(
    "6. Metabolic or stress context is not labeled",
    "as direct antimicrobial resistance without",
    "an explicit pre-existing annotation."
  ),
  
  paste(
    "7. Putative 3UTR-associated RNA is used instead",
    "of 3UTR-derived RNA until independent processing",
    "or boundary evidence is established."
  ),
  
  paste(
    "8. CLASH supports RNA-RNA contact but does not",
    "establish activation, repression, translational",
    "effect, or causality."
  ),
  
  paste(
    "9. Target mRNA differential abundance is",
    "orthogonal context; lack of DEG does not exclude",
    "post-transcriptional regulation."
  ),
  
  paste(
    "10. Ribo-seq without matched vancomycin exposure",
    "should not be used to claim vancomycin-induced",
    "translational regulation."
  )
)


writeLines(
  
  methodological_notes,
  
  con =
    file.path(
      
      audit_folder,
      
      "Script16_Revised_Methodological_Notes.txt"
    )
)


# ============================================================
# 45. SAVE SESSION INFORMATION
# ============================================================

writeLines(
  
  capture.output(
    sessionInfo()
  ),
  
  con =
    file.path(
      
      audit_folder,
      
      "sessionInfo.txt"
    )
)


# ============================================================
# 46. EXPECTED OUTPUT FILES
# ============================================================

expected_outputs <- c(
  
  file.path(
    
    table_folder,
    
    "Evidence_Layered_Network_Functionally_Annotated.csv"
  ),
  
  file.path(
    
    table_folder,
    
    "Manuscript_Functional_Annotation_Table.csv"
  ),
  
  file.path(
    
    table_folder,
    
    "DEG_Associated_Functional_Subset.csv"
  ),
  
  file.path(
    
    table_folder,
    
    "Benchmark_and_nrdF_Functional_Context.csv"
  ),
  
  file.path(
    
    table_folder,
    
    "Functional_System_Summary.csv"
  ),
  
  file.path(
    
    audit_folder,
    
    "Script16_Revised_Summary.csv"
  ),
  
  file.path(
    
    figure_folder,
    
    "Figure3A_Functional_Composition_600dpi.png"
  ),
  
  file.path(
    
    figure_folder,
    
    "Figure3B_CLASH_vs_Transcriptomic_Context_600dpi.png"
  )
)


# ============================================================
# 47. OUTPUT VERIFICATION
# ============================================================

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
    
    "Script16_Revised_Output_Verification.csv"
  )
)


# ============================================================
# 48. BASIC SANITY CHECKS
# ============================================================

if (nrow(annotated_network) != nrow(network)) {
  
  warning(
    paste(
      "The number of annotated rows differs",
      "from the input CLASH network."
    )
  )
}


if (nrow(network) != 211) {
  
  warning(
    paste0(
      "Expected approximately 211 collapsed CLASH edges, ",
      "but found ",
      nrow(network),
      ". Review Script 14 outputs."
    )
  )
}


if (
  sum(
    annotated_network$
    target_is_significant_deg,
    na.rm = TRUE
  ) != 10
) {
  
  warning(
    paste0(
      "Expected approximately 10 interactions with ",
      "significant target DEG, but found ",
      sum(
        annotated_network$
          target_is_significant_deg,
        na.rm = TRUE
      ),
      "."
    )
  )
}


if (nrdf_count != 2) {
  
  warning(
    paste0(
      "Expected 2 nrdF candidate interactions, but found ",
      nrdf_count,
      "."
    )
  )
}


# ============================================================
# 49. FINAL CONSOLE REPORT
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 16 REVISED COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================\n\n"
)


cat(
  "FUNCTIONAL ANNOTATION SUMMARY\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  script16_summary,
  n = Inf
)


cat(
  "\nCRITICAL INTERACTION AUDIT\n"
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
  "\nANNOTATION EVIDENCE SUMMARY\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  annotation_source_summary,
  n = Inf
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
  "IMPORTANT INTERPRETATION\n"
)

cat(
  "============================================\n"
)

cat(
  paste(
    "\nThe complete CLASH-supported network is retained.",
    "\nTarget differential expression is an independent layer.",
    "\nFunctional annotation does not increase CLASH evidence strength.",
    "\nNo integrated biological-priority score was calculated.",
    "\nNo Priority A/B/C categories were calculated.",
    "\nPutative 3UTR-associated RNAs remain candidate RNAs",
    "pending independent processing/boundary evidence.\n"
  )
)


cat(
  "\nMain output folder:\n",
  output_folder,
  "\n"
)


cat(
  "\nMain manuscript-ready table:\n",
  file.path(
    table_folder,
    "Manuscript_Functional_Annotation_Table.csv"
  ),
  "\n"
)


cat(
  "\nCritical benchmark/nrdF table:\n",
  file.path(
    table_folder,
    "Benchmark_and_nrdF_Functional_Context.csv"
  ),
  "\n"
)


cat(
  "\nSummary audit:\n",
  file.path(
    audit_folder,
    "Script16_Revised_Summary.csv"
  ),
  "\n"
)


cat(
  "\n============================================\n"
)

cat(
  "END OF SCRIPT 16 REVISED\n"
)

cat(
  "============================================\n"
)