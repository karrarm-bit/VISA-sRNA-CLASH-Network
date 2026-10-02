# ============================================================
# BacRegRNA Project
# Script 18
#
# Sequence extraction and structural-interaction preparation
# for the vancomycin-focused sRNA-mRNA regulatory network
#
# Main tasks:
# 1. Read validated interactions from Script 17.
# 2. Extract sRNA genomic sequences.
# 3. Extract target-mRNA regions around translation starts.
# 4. Generate pairwise and combined FASTA files.
# 5. Prepare IntaRNA commands.
# 6. Optionally run IntaRNA when the executable is available.
#
# Important:
# Structural predictions are supporting computational evidence.
# They do not independently validate functional regulation.
# ============================================================


# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

script17_table_folder <- file.path(
  project_folder,
  "17_network_validation",
  "tables"
)

validated_network_file <- file.path(
  script17_table_folder,
  "Vancomycin_sRNA_network_validated_and_enriched.csv"
)

master_kb_file <- file.path(
  project_folder,
  "13_CLASH_network",
  "tables",
  "JKD6008_master_gene_knowledgebase_CLASH_updated.csv"
)

reference_folder <- file.path(
  project_folder,
  "18_reference_genome"
)

output_folder <- file.path(
  project_folder,
  "18_structural_interaction"
)

sequence_folder <- file.path(
  output_folder,
  "sequences"
)

srna_folder <- file.path(
  sequence_folder,
  "sRNA"
)

target_folder <- file.path(
  sequence_folder,
  "target_mRNA"
)

pair_folder <- file.path(
  sequence_folder,
  "interaction_pairs"
)

intarna_folder <- file.path(
  output_folder,
  "IntaRNA"
)

intarna_result_folder <- file.path(
  intarna_folder,
  "results"
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
  reference_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  srna_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  target_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  pair_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  intarna_result_folder,
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
# 2. Analysis settings
# ------------------------------------------------------------

# Target region around the predicted translation start:
# 150 nt upstream and 100 nt downstream.

upstream_nt <- 150L
downstream_nt <- 100L

# Minimum acceptable sequence lengths.

minimum_srna_length <- 20L
minimum_target_length <- 50L

# Run IntaRNA automatically only when its executable exists.

run_IntaRNA_if_available <- TRUE


# ------------------------------------------------------------
# 3. Install and load packages
# ------------------------------------------------------------

cran_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "tibble",
  "purrr"
)

installed_packages <- rownames(
  installed.packages()
)

missing_cran_packages <- setdiff(
  cran_packages,
  installed_packages
)

if (length(missing_cran_packages) > 0) {
  
  install.packages(
    missing_cran_packages,
    dependencies = TRUE
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  
  install.packages("BiocManager")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  
  BiocManager::install(
    "Biostrings",
    ask = FALSE,
    update = FALSE
  )
}

library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(purrr)
library(Biostrings)


# ============================================================
# HELPER FUNCTIONS
# ============================================================


# ------------------------------------------------------------
# 4. Safe character conversion
# ------------------------------------------------------------

safe_character <- function(x) {
  
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
          
          as.character(value[1])
        },
        FUN.VALUE = character(1)
      )
    )
  }
  
  as.character(x)
}


# ------------------------------------------------------------
# 5. Safe numeric conversion
# ------------------------------------------------------------

safe_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


# ------------------------------------------------------------
# 6. Safe logical conversion
# ------------------------------------------------------------

safe_logical <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  value <- stringr::str_to_lower(
    stringr::str_squish(
      safe_character(x)
    )
  )
  
  dplyr::case_when(
    
    value %in% c(
      "true",
      "t",
      "1",
      "yes",
      "y"
    ) ~ TRUE,
    
    value %in% c(
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
# 7. Clean column names
# ------------------------------------------------------------

clean_column_names <- function(x) {
  
  cleaned <- safe_character(x)
  
  cleaned <- stringr::str_to_lower(cleaned)
  
  cleaned <- stringr::str_replace_all(
    cleaned,
    "[^a-z0-9]+",
    "_"
  )
  
  cleaned <- stringr::str_replace_all(
    cleaned,
    "^_+|_+$",
    ""
  )
  
  make.unique(
    cleaned,
    sep = "_"
  )
}


# ------------------------------------------------------------
# 8. Add missing column
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
# 9. Normalize identifier
# ------------------------------------------------------------

normalize_identifier <- function(x) {
  
  x <- safe_character(x)
  
  x <- stringr::str_squish(x)
  
  x <- stringr::str_replace_all(
    x,
    "[^A-Za-z0-9]+",
    "_"
  )
  
  x <- stringr::str_replace_all(
    x,
    "^_+|_+$",
    ""
  )
  
  stringr::str_to_lower(x)
}


# ------------------------------------------------------------
# 10. Create safe filename
# ------------------------------------------------------------

safe_filename <- function(x) {
  
  x <- safe_character(x)
  
  x <- stringr::str_replace_all(
    x,
    "[^A-Za-z0-9_-]+",
    "_"
  )
  
  x <- stringr::str_replace_all(
    x,
    "_+",
    "_"
  )
  
  x <- stringr::str_replace_all(
    x,
    "^_|_$",
    ""
  )
  
  x
}


# ------------------------------------------------------------
# 11. Extract sRNA coordinates from source_node
#
# Expected source-node ending:
# 631699-631784:-
# ------------------------------------------------------------

extract_srna_start <- function(x) {
  
  coordinate_text <- stringr::str_extract(
    safe_character(x),
    "[0-9]+-[0-9]+(?=:[+-]$)"
  )
  
  safe_numeric(
    stringr::str_extract(
      coordinate_text,
      "^[0-9]+"
    )
  )
}


extract_srna_end <- function(x) {
  
  coordinate_text <- stringr::str_extract(
    safe_character(x),
    "[0-9]+-[0-9]+(?=:[+-]$)"
  )
  
  safe_numeric(
    stringr::str_extract(
      coordinate_text,
      "[0-9]+$"
    )
  )
}


extract_srna_strand <- function(x) {
  
  stringr::str_extract(
    safe_character(x),
    "[+-]$"
  )
}


# ------------------------------------------------------------
# 12. Clamp genomic coordinates
# ------------------------------------------------------------

clamp_coordinate <- function(
    position,
    lower_limit,
    upper_limit
) {
  
  pmax(
    lower_limit,
    pmin(
      position,
      upper_limit
    )
  )
}


# ------------------------------------------------------------
# 13. Extract sequence from genome
# ------------------------------------------------------------

extract_genomic_sequence <- function(
    genome_sequence,
    start_position,
    end_position,
    strand = "+"
) {
  
  genome_length <- length(
    genome_sequence
  )
  
  if (
    is.na(start_position) ||
    is.na(end_position)
  ) {
    
    return(
      Biostrings::DNAString("")
    )
  }
  
  start_position <- as.integer(
    round(start_position)
  )
  
  end_position <- as.integer(
    round(end_position)
  )
  
  start_position <- clamp_coordinate(
    start_position,
    1L,
    genome_length
  )
  
  end_position <- clamp_coordinate(
    end_position,
    1L,
    genome_length
  )
  
  interval_start <- min(
    start_position,
    end_position
  )
  
  interval_end <- max(
    start_position,
    end_position
  )
  
  sequence_result <- Biostrings::subseq(
    genome_sequence,
    start = interval_start,
    end = interval_end
  )
  
  if (
    !is.na(strand) &&
    strand == "-"
  ) {
    
    sequence_result <- Biostrings::reverseComplement(
      sequence_result
    )
  }
  
  sequence_result
}


# ------------------------------------------------------------
# 14. Extract mRNA interaction region
#
# Plus strand:
# gene start - upstream to gene start + downstream
#
# Minus strand:
# gene end - downstream to gene end + upstream,
# followed by reverse complementation.
# ------------------------------------------------------------

extract_target_window <- function(
    genome_sequence,
    gene_start,
    gene_end,
    strand,
    upstream_length = 150L,
    downstream_length = 100L
) {
  
  if (
    is.na(gene_start) ||
    is.na(gene_end) ||
    is.na(strand)
  ) {
    
    return(
      list(
        sequence = Biostrings::DNAString(""),
        window_start = NA_integer_,
        window_end = NA_integer_,
        translation_start = NA_integer_
      )
    )
  }
  
  gene_start <- as.integer(
    round(gene_start)
  )
  
  gene_end <- as.integer(
    round(gene_end)
  )
  
  if (strand == "+") {
    
    translation_start <- min(
      gene_start,
      gene_end
    )
    
    window_start <- translation_start -
      upstream_length
    
    window_end <- translation_start +
      downstream_length
    
  } else {
    
    translation_start <- max(
      gene_start,
      gene_end
    )
    
    window_start <- translation_start -
      downstream_length
    
    window_end <- translation_start +
      upstream_length
  }
  
  target_sequence <- extract_genomic_sequence(
    genome_sequence = genome_sequence,
    start_position = window_start,
    end_position = window_end,
    strand = strand
  )
  
  list(
    sequence = target_sequence,
    window_start = max(
      1L,
      window_start
    ),
    window_end = min(
      length(genome_sequence),
      window_end
    ),
    translation_start = translation_start
  )
}


# ------------------------------------------------------------
# 15. Sequence composition
# ------------------------------------------------------------

calculate_gc_percent <- function(sequence_value) {
  
  sequence_text <- toupper(
    as.character(sequence_value)
  )
  
  if (
    is.na(sequence_text) ||
    nchar(sequence_text) == 0
  ) {
    
    return(NA_real_)
  }
  
  bases <- strsplit(
    sequence_text,
    split = ""
  )[[1]]
  
  valid_bases <- bases[
    bases %in% c(
      "A",
      "C",
      "G",
      "T"
    )
  ]
  
  if (length(valid_bases) == 0) {
    return(NA_real_)
  }
  
  100 *
    sum(
      valid_bases %in% c(
        "G",
        "C"
      )
    ) /
    length(valid_bases)
}


# ------------------------------------------------------------
# 16. Write one FASTA record
# ------------------------------------------------------------

write_single_fasta <- function(
    sequence_value,
    sequence_name,
    output_file
) {
  
  fasta_object <- Biostrings::DNAStringSet(
    sequence_value
  )
  
  names(fasta_object) <- sequence_name
  
  Biostrings::writeXStringSet(
    fasta_object,
    filepath = output_file,
    format = "fasta"
  )
}


# ============================================================
# VERIFY INPUT FILES
# ============================================================


# ------------------------------------------------------------
# 17. Check Script 17 and master files
# ------------------------------------------------------------

required_files <- c(
  validated_network_file,
  master_kb_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    "Required files were not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------
# 18. Locate genome FASTA file
# ------------------------------------------------------------

genome_candidates <- list.files(
  reference_folder,
  pattern = "\\.(fna|fa|fasta)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(genome_candidates) == 0) {
  
  stop(
    paste0(
      "\nNo genome FASTA file was found.\n\n",
      "Place the JKD6008 genome file inside:\n",
      reference_folder,
      "\n\nAccepted extensions: .fna, .fa or .fasta\n"
    )
  )
}

genome_file <- genome_candidates[1]

cat(
  "\nGenome file selected:\n",
  genome_file,
  "\n"
)


# ============================================================
# READ DATA
# ============================================================


# ------------------------------------------------------------
# 19. Read validated network
# ------------------------------------------------------------

validated_network <- readr::read_csv(
  validated_network_file,
  show_col_types = FALSE
)

master_kb <- readr::read_csv(
  master_kb_file,
  show_col_types = FALSE
)

names(validated_network) <- clean_column_names(
  names(validated_network)
)

names(master_kb) <- clean_column_names(
  names(master_kb)
)

if (nrow(validated_network) == 0) {
  
  stop(
    "The validated interaction table contains zero rows."
  )
}


# ------------------------------------------------------------
# 20. Read genome
# ------------------------------------------------------------

genome_set <- Biostrings::readDNAStringSet(
  genome_file,
  format = "fasta"
)

if (length(genome_set) == 0) {
  
  stop(
    "No sequences were read from the genome FASTA file."
  )
}

genome_lengths <- Biostrings::width(
  genome_set
)

primary_sequence_index <- which.max(
  genome_lengths
)

genome_sequence <- genome_set[
  primary_sequence_index
][[1]]

genome_name <- names(
  genome_set
)[primary_sequence_index]

genome_length <- length(
  genome_sequence
)

cat(
  "\nGenome sequence used:",
  genome_name,
  "\nGenome length:",
  genome_length,
  "nt\n"
)


# ============================================================
# STANDARDIZE COLUMNS
# ============================================================


# ------------------------------------------------------------
# 21. Required network columns
# ------------------------------------------------------------

required_network_columns <- list(
  
  source_node = NA_character_,
  target_node = NA_character_,
  validated_gene_id = NA_character_,
  
  curated_source_label = NA_character_,
  final_gene_symbol = NA_character_,
  
  total_confidence_score = NA_real_,
  integrated_biological_priority_score = NA_real_,
  validation_adjusted_priority_score = NA_real_,
  
  biological_priority = NA_character_,
  validation_priority = NA_character_,
  
  identifier_validation_status = NA_character_,
  requires_manual_review = FALSE
)

for (
  column_name in names(
    required_network_columns
  )
) {
  
  validated_network <- add_missing_column(
    validated_network,
    column_name,
    required_network_columns[[column_name]]
  )
}


# ------------------------------------------------------------
# 22. Required master columns
# ------------------------------------------------------------

required_master_columns <- list(
  
  gene_id = NA_character_,
  current_locus_tag = NA_character_,
  gene_symbol = NA_character_,
  genomic_start = NA_real_,
  genomic_end = NA_real_,
  strand = NA_character_
)

for (
  column_name in names(
    required_master_columns
  )
) {
  
  master_kb <- add_missing_column(
    master_kb,
    column_name,
    required_master_columns[[column_name]]
  )
}


# ------------------------------------------------------------
# 23. Standardize master coordinates
# ------------------------------------------------------------

master_coordinates <- master_kb |>
  dplyr::mutate(
    
    gene_id =
      safe_character(
        .data$gene_id
      ),
    
    gene_id_normalized =
      normalize_identifier(
        .data$gene_id
      ),
    
    genomic_start =
      safe_numeric(
        .data$genomic_start
      ),
    
    genomic_end =
      safe_numeric(
        .data$genomic_end
      ),
    
    strand =
      safe_character(
        .data$strand
      )
  ) |>
  dplyr::transmute(
    
    validated_gene_id_normalized =
      .data$gene_id_normalized,
    
    master_gene_id =
      .data$gene_id,
    
    target_locus_tag =
      safe_character(
        .data$current_locus_tag
      ),
    
    target_gene_symbol_master =
      safe_character(
        .data$gene_symbol
      ),
    
    target_genomic_start =
      .data$genomic_start,
    
    target_genomic_end =
      .data$genomic_end,
    
    target_strand =
      .data$strand
  ) |>
  dplyr::filter(
    !is.na(
      .data$validated_gene_id_normalized
    )
  ) |>
  dplyr::distinct(
    .data$validated_gene_id_normalized,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# 24. Join interactions to coordinates
# ------------------------------------------------------------

interaction_sequences <- validated_network |>
  dplyr::mutate(
    
    interaction_id =
      paste0(
        "INTPAIR_",
        stringr::str_pad(
          dplyr::row_number(),
          width = 3,
          side = "left",
          pad = "0"
        )
      ),
    
    source_node =
      safe_character(
        .data$source_node
      ),
    
    target_node =
      safe_character(
        .data$target_node
      ),
    
    validated_gene_id =
      dplyr::coalesce(
        safe_character(
          .data$validated_gene_id
        ),
        safe_character(
          .data$target_node
        )
      ),
    
    validated_gene_id_normalized =
      normalize_identifier(
        .data$validated_gene_id
      ),
    
    srna_start =
      extract_srna_start(
        .data$source_node
      ),
    
    srna_end =
      extract_srna_end(
        .data$source_node
      ),
    
    srna_strand =
      extract_srna_strand(
        .data$source_node
      )
  ) |>
  dplyr::left_join(
    master_coordinates,
    by = "validated_gene_id_normalized"
  )


# ============================================================
# EXTRACT SEQUENCES
# ============================================================


# ------------------------------------------------------------
# 25. Extract all sRNA and target sequences
# ------------------------------------------------------------

sequence_results <- vector(
  mode = "list",
  length = nrow(
    interaction_sequences
  )
)

for (
  row_index in seq_len(
    nrow(interaction_sequences)
  )
) {
  
  current_row <- interaction_sequences[
    row_index,
  ]
  
  srna_sequence <- extract_genomic_sequence(
    
    genome_sequence =
      genome_sequence,
    
    start_position =
      current_row$srna_start,
    
    end_position =
      current_row$srna_end,
    
    strand =
      current_row$srna_strand
  )
  
  target_result <- extract_target_window(
    
    genome_sequence =
      genome_sequence,
    
    gene_start =
      current_row$target_genomic_start,
    
    gene_end =
      current_row$target_genomic_end,
    
    strand =
      current_row$target_strand,
    
    upstream_length =
      upstream_nt,
    
    downstream_length =
      downstream_nt
  )
  
  sequence_results[[row_index]] <- list(
    
    srna_sequence =
      as.character(
        srna_sequence
      ),
    
    target_sequence =
      as.character(
        target_result$sequence
      ),
    
    target_window_start =
      target_result$window_start,
    
    target_window_end =
      target_result$window_end,
    
    target_translation_start =
      target_result$translation_start
  )
}


# ------------------------------------------------------------
# 26. Add extracted sequences
# ------------------------------------------------------------

interaction_sequences$srna_sequence <-
  vapply(
    sequence_results,
    function(x) {
      x$srna_sequence
    },
    FUN.VALUE = character(1)
  )

interaction_sequences$target_sequence <-
  vapply(
    sequence_results,
    function(x) {
      x$target_sequence
    },
    FUN.VALUE = character(1)
  )

interaction_sequences$target_window_start <-
  vapply(
    sequence_results,
    function(x) {
      
      if (is.null(x$target_window_start)) {
        return(NA_real_)
      }
      
      as.numeric(
        x$target_window_start
      )
    },
    FUN.VALUE = numeric(1)
  )

interaction_sequences$target_window_end <-
  vapply(
    sequence_results,
    function(x) {
      
      if (is.null(x$target_window_end)) {
        return(NA_real_)
      }
      
      as.numeric(
        x$target_window_end
      )
    },
    FUN.VALUE = numeric(1)
  )

interaction_sequences$target_translation_start <-
  vapply(
    sequence_results,
    function(x) {
      
      if (is.null(x$target_translation_start)) {
        return(NA_real_)
      }
      
      as.numeric(
        x$target_translation_start
      )
    },
    FUN.VALUE = numeric(1)
  )


# ------------------------------------------------------------
# 27. Calculate sequence metrics
# ------------------------------------------------------------

interaction_sequences <- interaction_sequences |>
  dplyr::mutate(
    
    srna_length =
      nchar(
        .data$srna_sequence
      ),
    
    target_sequence_length =
      nchar(
        .data$target_sequence
      ),
    
    srna_GC_percent =
      vapply(
        .data$srna_sequence,
        calculate_gc_percent,
        FUN.VALUE = numeric(1)
      ),
    
    target_GC_percent =
      vapply(
        .data$target_sequence,
        calculate_gc_percent,
        FUN.VALUE = numeric(1)
      ),
    
    srna_coordinates_available =
      !is.na(
        .data$srna_start
      ) &
      !is.na(
        .data$srna_end
      ) &
      .data$srna_strand %in%
      c(
        "+",
        "-"
      ),
    
    target_coordinates_available =
      !is.na(
        .data$target_genomic_start
      ) &
      !is.na(
        .data$target_genomic_end
      ) &
      .data$target_strand %in%
      c(
        "+",
        "-"
      ),
    
    srna_sequence_valid =
      .data$srna_length >=
      minimum_srna_length &
      stringr::str_detect(
        .data$srna_sequence,
        "^[ACGTNacgtn]+$"
      ),
    
    target_sequence_valid =
      .data$target_sequence_length >=
      minimum_target_length &
      stringr::str_detect(
        .data$target_sequence,
        "^[ACGTNacgtn]+$"
      ),
    
    interaction_sequence_status =
      dplyr::case_when(
        
        .data$srna_sequence_valid &
          .data$target_sequence_valid ~
          "Ready for structural prediction",
        
        !.data$srna_coordinates_available ~
          "sRNA coordinates unresolved",
        
        !.data$target_coordinates_available ~
          "Target coordinates unresolved",
        
        !.data$srna_sequence_valid ~
          "Invalid or short sRNA sequence",
        
        !.data$target_sequence_valid ~
          "Invalid or short target sequence",
        
        TRUE ~
          "Sequence review required"
      )
  )


# ============================================================
# WRITE FASTA FILES
# ============================================================


# ------------------------------------------------------------
# 28. Write individual and pairwise FASTA files
# ------------------------------------------------------------

fasta_manifest <- vector(
  mode = "list",
  length = nrow(
    interaction_sequences
  )
)

for (
  row_index in seq_len(
    nrow(interaction_sequences)
  )
) {
  
  current_row <- interaction_sequences[
    row_index,
  ]
  
  source_label_safe <- safe_filename(
    current_row$curated_source_label
  )
  
  target_label_safe <- safe_filename(
    current_row$final_gene_symbol
  )
  
  if (
    is.na(source_label_safe) ||
    source_label_safe == ""
  ) {
    
    source_label_safe <-
      current_row$interaction_id
  }
  
  if (
    is.na(target_label_safe) ||
    target_label_safe == ""
  ) {
    
    target_label_safe <-
      safe_filename(
        current_row$validated_gene_id
      )
  }
  
  pair_prefix <- paste0(
    current_row$interaction_id,
    "_",
    source_label_safe,
    "_to_",
    target_label_safe
  )
  
  srna_fasta_file <- file.path(
    srna_folder,
    paste0(
      pair_prefix,
      "_sRNA.fa"
    )
  )
  
  target_fasta_file <- file.path(
    target_folder,
    paste0(
      pair_prefix,
      "_target.fa"
    )
  )
  
  pair_fasta_file <- file.path(
    pair_folder,
    paste0(
      pair_prefix,
      "_pair.fa"
    )
  )
  
  if (
    current_row$interaction_sequence_status ==
    "Ready for structural prediction"
  ) {
    
    srna_name <- paste0(
      current_row$interaction_id,
      "|sRNA|",
      current_row$curated_source_label,
      "|",
      current_row$srna_start,
      "-",
      current_row$srna_end,
      "|",
      current_row$srna_strand
    )
    
    target_name <- paste0(
      current_row$interaction_id,
      "|target|",
      current_row$final_gene_symbol,
      "|",
      current_row$target_window_start,
      "-",
      current_row$target_window_end,
      "|",
      current_row$target_strand
    )
    
    write_single_fasta(
      sequence_value =
        current_row$srna_sequence,
      
      sequence_name =
        srna_name,
      
      output_file =
        srna_fasta_file
    )
    
    write_single_fasta(
      sequence_value =
        current_row$target_sequence,
      
      sequence_name =
        target_name,
      
      output_file =
        target_fasta_file
    )
    
    pair_set <- Biostrings::DNAStringSet(
      c(
        current_row$srna_sequence,
        current_row$target_sequence
      )
    )
    
    names(pair_set) <- c(
      srna_name,
      target_name
    )
    
    Biostrings::writeXStringSet(
      pair_set,
      filepath = pair_fasta_file,
      format = "fasta"
    )
    
  } else {
    
    srna_fasta_file <- NA_character_
    target_fasta_file <- NA_character_
    pair_fasta_file <- NA_character_
  }
  
  fasta_manifest[[row_index]] <- tibble::tibble(
    
    interaction_id =
      current_row$interaction_id,
    
    regulatory_RNA =
      current_row$curated_source_label,
    
    target_gene =
      current_row$final_gene_symbol,
    
    status =
      current_row$interaction_sequence_status,
    
    srna_fasta =
      srna_fasta_file,
    
    target_fasta =
      target_fasta_file,
    
    pair_fasta =
      pair_fasta_file
  )
}

fasta_manifest <- dplyr::bind_rows(
  fasta_manifest
)


# ------------------------------------------------------------
# 29. Write combined sRNA FASTA
# ------------------------------------------------------------

ready_sequences <- interaction_sequences |>
  dplyr::filter(
    .data$interaction_sequence_status ==
      "Ready for structural prediction"
  )

if (nrow(ready_sequences) > 0) {
  
  combined_srna_set <- Biostrings::DNAStringSet(
    ready_sequences$srna_sequence
  )
  
  names(combined_srna_set) <- paste0(
    ready_sequences$interaction_id,
    "|",
    ready_sequences$curated_source_label
  )
  
  Biostrings::writeXStringSet(
    
    combined_srna_set,
    
    filepath = file.path(
      sequence_folder,
      "All_interaction_sRNAs.fa"
    ),
    
    format = "fasta"
  )
  
  
  combined_target_set <- Biostrings::DNAStringSet(
    ready_sequences$target_sequence
  )
  
  names(combined_target_set) <- paste0(
    ready_sequences$interaction_id,
    "|",
    ready_sequences$final_gene_symbol
  )
  
  Biostrings::writeXStringSet(
    
    combined_target_set,
    
    filepath = file.path(
      sequence_folder,
      "All_target_interaction_windows.fa"
    ),
    
    format = "fasta"
  )
}


# ============================================================
# PREPARE INTARNA COMMANDS
# ============================================================


# ------------------------------------------------------------
# 30. Build IntaRNA command table
# ------------------------------------------------------------

intarna_command_table <- fasta_manifest |>
  dplyr::filter(
    .data$status ==
      "Ready for structural prediction"
  ) |>
  dplyr::mutate(
    
    result_file =
      file.path(
        intarna_result_folder,
        paste0(
          .data$interaction_id,
          "_IntaRNA.csv"
        )
      ),
    
    interaction_output_file =
      file.path(
        intarna_result_folder,
        paste0(
          .data$interaction_id,
          "_interaction.txt"
        )
      ),
    
    command =
      paste(
        
        "IntaRNA",
        
        paste0(
          "--query=\"",
          .data$srna_fasta,
          "\""
        ),
        
        paste0(
          "--target=\"",
          .data$target_fasta,
          "\""
        ),
        
        "--outMode=C",
        
        "--outCsvCols=id1,id2,start1,end1,start2,end2,E",
        
        paste0(
          "--out=\"",
          .data$result_file,
          "\""
        )
      )
  )


# ------------------------------------------------------------
# 31. Save IntaRNA shell commands
# ------------------------------------------------------------

if (nrow(intarna_command_table) > 0) {
  
  writeLines(
    
    c(
      "#!/usr/bin/env bash",
      "",
      "# IntaRNA commands generated by Script 18",
      "",
      intarna_command_table$command
    ),
    
    con = file.path(
      intarna_folder,
      "run_all_IntaRNA.sh"
    )
  )
  
  writeLines(
    
    intarna_command_table$command,
    
    con = file.path(
      intarna_folder,
      "IntaRNA_commands.txt"
    )
  )
}


# ============================================================
# OPTIONAL INTARNA EXECUTION
# ============================================================


# ------------------------------------------------------------
# 32. Check IntaRNA availability
# ------------------------------------------------------------

intarna_executable <- Sys.which(
  "IntaRNA"
)

intarna_available <- (
  nchar(
    intarna_executable
  ) > 0
)

cat(
  "\nIntaRNA executable detected:",
  intarna_available,
  "\n"
)

if (intarna_available) {
  
  cat(
    "IntaRNA path:",
    intarna_executable,
    "\n"
  )
}


# ------------------------------------------------------------
# 33. Run IntaRNA if available
# ------------------------------------------------------------

intarna_run_report <- tibble::tibble()

if (
  run_IntaRNA_if_available &&
  intarna_available &&
  nrow(intarna_command_table) > 0
) {
  
  intarna_run_list <- vector(
    mode = "list",
    length = nrow(
      intarna_command_table
    )
  )
  
  for (
    row_index in seq_len(
      nrow(intarna_command_table)
    )
  ) {
    
    current_command <-
      intarna_command_table[
        row_index,
      ]
    
    command_arguments <- c(
      
      paste0(
        "--query=",
        current_command$srna_fasta
      ),
      
      paste0(
        "--target=",
        current_command$target_fasta
      ),
      
      "--outMode=C",
      
      "--outCsvCols=id1,id2,start1,end1,start2,end2,E",
      
      paste0(
        "--out=",
        current_command$result_file
      )
    )
    
    command_status <- tryCatch(
      
      {
        
        system2(
          command =
            intarna_executable,
          
          args =
            command_arguments,
          
          stdout = TRUE,
          
          stderr = TRUE
        )
        
        if (
          file.exists(
            current_command$result_file
          )
        ) {
          "Completed"
        } else {
          "Command executed but result file was not created"
        }
      },
      
      error = function(e) {
        
        paste0(
          "Failed: ",
          conditionMessage(e)
        )
      }
    )
    
    intarna_run_list[[row_index]] <-
      tibble::tibble(
        
        interaction_id =
          current_command$interaction_id,
        
        regulatory_RNA =
          current_command$regulatory_RNA,
        
        target_gene =
          current_command$target_gene,
        
        result_file =
          current_command$result_file,
        
        run_status =
          paste(
            command_status,
            collapse = " | "
          )
      )
  }
  
  intarna_run_report <- dplyr::bind_rows(
    intarna_run_list
  )
  
} else {
  
  intarna_run_report <- intarna_command_table |>
    dplyr::transmute(
      
      interaction_id =
        .data$interaction_id,
      
      regulatory_RNA =
        .data$regulatory_RNA,
      
      target_gene =
        .data$target_gene,
      
      result_file =
        .data$result_file,
      
      run_status =
        dplyr::case_when(
          
          !intarna_available ~
            "IntaRNA not installed; FASTA files and commands prepared",
          
          !run_IntaRNA_if_available ~
            "Automatic IntaRNA execution disabled",
          
          TRUE ~
            "No interaction available for execution"
        )
    )
}


# ============================================================
# IMPORT AVAILABLE INTARNA RESULTS
# ============================================================


# ------------------------------------------------------------
# 34. Read completed IntaRNA result files
# ------------------------------------------------------------

intarna_result_list <- list()

if (nrow(intarna_command_table) > 0) {
  
  for (
    row_index in seq_len(
      nrow(intarna_command_table)
    )
  ) {
    
    current_row <- intarna_command_table[
      row_index,
    ]
    
    if (
      file.exists(
        current_row$result_file
      ) &&
      file.info(
        current_row$result_file
      )$size > 0
    ) {
      
      result_data <- tryCatch(
        
        readr::read_delim(
          
          current_row$result_file,
          
          delim = ";",
          
          show_col_types = FALSE,
          
          progress = FALSE
        ),
        
        error = function(e) {
          
          tryCatch(
            
            readr::read_csv(
              current_row$result_file,
              show_col_types = FALSE
            ),
            
            error = function(e2) {
              NULL
            }
          )
        }
      )
      
      if (
        !is.null(result_data) &&
        nrow(result_data) > 0
      ) {
        
        result_data$interaction_id <-
          current_row$interaction_id
        
        result_data$regulatory_RNA <-
          current_row$regulatory_RNA
        
        result_data$target_gene <-
          current_row$target_gene
        
        intarna_result_list[
          [
            length(
              intarna_result_list
            ) + 1
          ]
        ] <- list(
          result_data
        )
      }
    }
  }
}

if (length(intarna_result_list) > 0) {
  
  combined_intarna_results <-
    dplyr::bind_rows(
      intarna_result_list
    )
  
} else {
  
  combined_intarna_results <-
    tibble::tibble(
      
      interaction_id = character(),
      regulatory_RNA = character(),
      target_gene = character(),
      id1 = character(),
      id2 = character(),
      start1 = numeric(),
      end1 = numeric(),
      start2 = numeric(),
      end2 = numeric(),
      E = numeric()
    )
}


# ============================================================
# SAVE TABLES AND AUDIT FILES
# ============================================================


# ------------------------------------------------------------
# 35. Save sequence metadata
# ------------------------------------------------------------

readr::write_csv(
  
  interaction_sequences,
  
  file.path(
    table_folder,
    "Interaction_sequence_metadata.csv"
  )
)

readr::write_csv(
  
  fasta_manifest,
  
  file.path(
    table_folder,
    "Interaction_FASTA_manifest.csv"
  )
)

readr::write_csv(
  
  intarna_command_table,
  
  file.path(
    table_folder,
    "IntaRNA_command_manifest.csv"
  )
)

readr::write_csv(
  
  intarna_run_report,
  
  file.path(
    table_folder,
    "IntaRNA_run_report.csv"
  )
)

readr::write_csv(
  
  combined_intarna_results,
  
  file.path(
    table_folder,
    "Combined_IntaRNA_predictions.csv"
  )
)


# ------------------------------------------------------------
# 36. Sequence extraction audit
# ------------------------------------------------------------

sequence_audit <- interaction_sequences |>
  dplyr::transmute(
    
    interaction_id =
      .data$interaction_id,
    
    regulatory_RNA =
      .data$curated_source_label,
    
    target_gene =
      .data$final_gene_symbol,
    
    srna_coordinates =
      paste0(
        .data$srna_start,
        "-",
        .data$srna_end,
        ":",
        .data$srna_strand
      ),
    
    srna_length =
      .data$srna_length,
    
    srna_GC_percent =
      round(
        .data$srna_GC_percent,
        2
      ),
    
    target_coordinates =
      paste0(
        .data$target_window_start,
        "-",
        .data$target_window_end,
        ":",
        .data$target_strand
      ),
    
    target_sequence_length =
      .data$target_sequence_length,
    
    target_GC_percent =
      round(
        .data$target_GC_percent,
        2
      ),
    
    status =
      .data$interaction_sequence_status
  )

readr::write_csv(
  
  sequence_audit,
  
  file.path(
    audit_folder,
    "Sequence_extraction_audit.csv"
  )
)


# ------------------------------------------------------------
# 37. Script summary
# ------------------------------------------------------------

script18_summary <- tibble::tibble(
  
  metric = c(
    
    "Validated interactions read",
    
    "Interactions with sRNA coordinates",
    
    "Interactions with target coordinates",
    
    "Interactions ready for structural prediction",
    
    "Interactions requiring sequence review",
    
    "Individual sRNA FASTA files",
    
    "Individual target FASTA files",
    
    "IntaRNA executable available",
    
    "IntaRNA predictions imported"
  ),
  
  value = c(
    
    nrow(
      interaction_sequences
    ),
    
    sum(
      interaction_sequences$
        srna_coordinates_available,
      na.rm = TRUE
    ),
    
    sum(
      interaction_sequences$
        target_coordinates_available,
      na.rm = TRUE
    ),
    
    sum(
      interaction_sequences$
        interaction_sequence_status ==
        "Ready for structural prediction",
      na.rm = TRUE
    ),
    
    sum(
      interaction_sequences$
        interaction_sequence_status !=
        "Ready for structural prediction",
      na.rm = TRUE
    ),
    
    sum(
      !is.na(
        fasta_manifest$srna_fasta
      )
    ),
    
    sum(
      !is.na(
        fasta_manifest$target_fasta
      )
    ),
    
    as.integer(
      intarna_available
    ),
    
    dplyr::n_distinct(
      combined_intarna_results$
        interaction_id
    )
  )
)

readr::write_csv(
  
  script18_summary,
  
  file.path(
    audit_folder,
    "Script18_summary.csv"
  )
)


# ------------------------------------------------------------
# 38. Methodological notes
# ------------------------------------------------------------

methodological_notes <- c(
  
  "Script 18 methodological notes",
  
  "",
  
  paste0(
    "1. Target interaction windows were extracted from ",
    upstream_nt,
    " nt upstream to ",
    downstream_nt,
    " nt downstream of the predicted translation start."
  ),
  
  "2. Minus-strand sequences were reverse-complemented to preserve transcript orientation.",
  
  "3. sRNA coordinates were parsed from the source-node identifiers generated during CLASH processing.",
  
  "4. Target coordinates were obtained from the JKD6008 master gene knowledgebase.",
  
  "5. The longest FASTA sequence was treated as the primary chromosome or principal genome sequence.",
  
  "6. IntaRNA structural predictions are supporting computational evidence and do not establish functional regulation.",
  
  "7. Absence of the IntaRNA executable does not invalidate sequence extraction; FASTA files and executable commands are retained for later analysis.",
  
  "8. The extracted mRNA windows emphasize the translation-initiation region, which is commonly relevant to bacterial sRNA regulation.",
  
  "9. Genomic-coordinate conventions were treated as one-based inclusive coordinates.",
  
  "10. Sequence records failing minimum-length or nucleotide-content checks were excluded from automatic structural prediction."
)

writeLines(
  
  methodological_notes,
  
  con = file.path(
    audit_folder,
    "Script18_methodological_notes.txt"
  )
)


# ------------------------------------------------------------
# 39. Save session information
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
# 40. Verify important outputs
# ------------------------------------------------------------

expected_outputs <- c(
  
  file.path(
    table_folder,
    "Interaction_sequence_metadata.csv"
  ),
  
  file.path(
    table_folder,
    "Interaction_FASTA_manifest.csv"
  ),
  
  file.path(
    table_folder,
    "IntaRNA_command_manifest.csv"
  ),
  
  file.path(
    table_folder,
    "IntaRNA_run_report.csv"
  ),
  
  file.path(
    table_folder,
    "Combined_IntaRNA_predictions.csv"
  ),
  
  file.path(
    audit_folder,
    "Sequence_extraction_audit.csv"
  ),
  
  file.path(
    audit_folder,
    "Script18_summary.csv"
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
      )$size / 1024,
      2
    )
)

readr::write_csv(
  
  output_verification,
  
  file.path(
    audit_folder,
    "Script18_output_verification.csv"
  )
)


# ------------------------------------------------------------
# 41. Final console report
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 18 COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================\n\n"
)

cat(
  "Genome used:\n",
  genome_file,
  "\n\n"
)

cat(
  "Sequence extraction summary:\n"
)

print(
  script18_summary,
  n = Inf
)

cat(
  "\nInteraction sequence status:\n"
)

print(
  sequence_audit,
  n = Inf
)

cat(
  "\nFASTA files saved inside:\n",
  sequence_folder,
  "\n\n"
)

cat(
  "IntaRNA command file:\n",
  file.path(
    intarna_folder,
    "IntaRNA_commands.txt"
  ),
  "\n\n"
)

if (!intarna_available) {
  
  cat(
    "NOTE:\n",
    "IntaRNA was not detected on this computer.\n",
    "The sequence files and commands were prepared successfully.\n",
    "Structural prediction can be run later using IntaRNA ",
    "through Linux, WSL, Conda or another supported environment.\n\n"
  )
}

cat(
  "Combined prediction table:\n",
  file.path(
    table_folder,
    "Combined_IntaRNA_predictions.csv"
  ),
  "\n\n"
)

cat(
  "Next step:\n",
  "Inspect sequence extraction, then run or import IntaRNA ",
  "predictions and visualize the strongest sRNA-mRNA ",
  "binding regions.\n"
)
# ============================================================
# Script 18 — Fix 2
# Robust coordinate rescue without assuming columns exist
# ============================================================


# ------------------------------------------------------------
# 1. Ensure required columns exist before rescue
# ------------------------------------------------------------

required_interaction_columns <- list(
  
  master_gene_id = NA_character_,
  target_locus_tag = NA_character_,
  target_gene_symbol_master = NA_character_,
  
  target_genomic_start = NA_real_,
  target_genomic_end = NA_real_,
  target_strand = NA_character_,
  
  target_symbol_normalized = NA_character_
)

for (column_name in names(required_interaction_columns)) {
  
  if (!column_name %in% names(interaction_sequences)) {
    
    interaction_sequences[[column_name]] <- rep(
      required_interaction_columns[[column_name]],
      nrow(interaction_sequences)
    )
  }
}


# ------------------------------------------------------------
# 2. Ensure normalized target symbol exists
# ------------------------------------------------------------

interaction_sequences$target_symbol_normalized <-
  normalize_identifier(
    interaction_sequences$final_gene_symbol
  )


# ------------------------------------------------------------
# 3. Build symbol-based coordinate lookup
# ------------------------------------------------------------

symbol_coordinate_lookup <- master_coordinates |>
  dplyr::filter(
    !is.na(.data$symbol_normalized),
    .data$symbol_normalized != ""
  ) |>
  dplyr::distinct(
    .data$symbol_normalized,
    .keep_all = TRUE
  ) |>
  dplyr::transmute(
    
    target_symbol_normalized =
      .data$symbol_normalized,
    
    rescue_master_gene_id =
      .data$master_gene_id,
    
    rescue_locus_tag =
      .data$target_locus_tag,
    
    rescue_gene_symbol =
      .data$target_gene_symbol_master,
    
    rescue_genomic_start =
      .data$target_genomic_start,
    
    rescue_genomic_end =
      .data$target_genomic_end,
    
    rescue_strand =
      .data$target_strand
  )


# ------------------------------------------------------------
# 4. Join rescue coordinates
# ------------------------------------------------------------

interaction_sequences <- interaction_sequences |>
  dplyr::left_join(
    symbol_coordinate_lookup,
    by = "target_symbol_normalized"
  )


# ------------------------------------------------------------
# 5. Fill missing values safely
# ------------------------------------------------------------

interaction_sequences <- interaction_sequences |>
  dplyr::mutate(
    
    master_gene_id =
      dplyr::coalesce(
        as.character(.data$master_gene_id),
        as.character(.data$rescue_master_gene_id)
      ),
    
    target_locus_tag =
      dplyr::coalesce(
        as.character(.data$target_locus_tag),
        as.character(.data$rescue_locus_tag)
      ),
    
    target_gene_symbol_master =
      dplyr::coalesce(
        as.character(.data$target_gene_symbol_master),
        as.character(.data$rescue_gene_symbol),
        as.character(.data$final_gene_symbol)
      ),
    
    target_genomic_start =
      dplyr::coalesce(
        suppressWarnings(
          as.numeric(.data$target_genomic_start)
        ),
        suppressWarnings(
          as.numeric(.data$rescue_genomic_start)
        )
      ),
    
    target_genomic_end =
      dplyr::coalesce(
        suppressWarnings(
          as.numeric(.data$target_genomic_end)
        ),
        suppressWarnings(
          as.numeric(.data$rescue_genomic_end)
        )
      ),
    
    target_strand =
      dplyr::coalesce(
        as.character(.data$target_strand),
        as.character(.data$rescue_strand)
      )
  ) |>
  dplyr::select(
    -dplyr::starts_with("rescue_")
  )


# ------------------------------------------------------------
# 6. Standardize strand once more
# ------------------------------------------------------------

interaction_sequences <- interaction_sequences |>
  dplyr::mutate(
    
    target_strand =
      dplyr::case_when(
        
        stringr::str_to_lower(
          stringr::str_squish(
            as.character(.data$target_strand)
          )
        ) %in% c(
          "+",
          "plus",
          "forward",
          "1",
          "positive"
        ) ~ "+",
        
        stringr::str_to_lower(
          stringr::str_squish(
            as.character(.data$target_strand)
          )
        ) %in% c(
          "-",
          "minus",
          "reverse",
          "-1",
          "negative"
        ) ~ "-",
        
        TRUE ~ NA_character_
      )
  )


# ------------------------------------------------------------
# 7. Check recovered coordinates
# ------------------------------------------------------------

coordinate_check <- interaction_sequences |>
  dplyr::transmute(
    
    interaction_id =
      .data$interaction_id,
    
    regulatory_RNA =
      .data$curated_source_label,
    
    target_gene =
      .data$final_gene_symbol,
    
    target_gene_id =
      .data$validated_gene_id,
    
    gene_start =
      .data$target_genomic_start,
    
    gene_end =
      .data$target_genomic_end,
    
    strand =
      .data$target_strand,
    
    coordinate_status =
      dplyr::case_when(
        
        !is.na(.data$target_genomic_start) &
          !is.na(.data$target_genomic_end) &
          .data$target_strand %in% c("+", "-") ~
          
          "Coordinates available",
        
        TRUE ~
          "Coordinates missing"
      )
  )

print(
  coordinate_check,
  n = Inf
)

cat(
  "\nTargets with coordinates:",
  sum(
    coordinate_check$coordinate_status ==
      "Coordinates available"
  ),
  "of",
  nrow(coordinate_check),
  "\n"
)
# ============================================================
# Script 18 continuation
# Continue after:
# Targets with coordinates: 10 of 10
# ============================================================


# ------------------------------------------------------------
# 1. Verify required objects
# ------------------------------------------------------------

required_objects <- c(
  "interaction_sequences",
  "genome_sequence",
  "srna_folder",
  "target_folder",
  "pair_folder",
  "sequence_folder",
  "intarna_folder",
  "intarna_result_folder",
  "table_folder",
  "audit_folder",
  "upstream_nt",
  "downstream_nt",
  "minimum_srna_length",
  "minimum_target_length",
  "run_IntaRNA_if_available"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    FUN.VALUE = logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0) {
  
  stop(
    "Required objects are missing:\n",
    paste(
      missing_objects,
      collapse = "\n"
    ),
    "\nRun Script 18 up to the coordinate-recovery section first."
  )
}


# ------------------------------------------------------------
# 2. Ensure helper functions exist
# ------------------------------------------------------------

required_functions <- c(
  "extract_genomic_sequence",
  "extract_target_window",
  "calculate_gc_percent",
  "safe_filename",
  "write_single_fasta"
)

missing_functions <- required_functions[
  !vapply(
    required_functions,
    exists,
    FUN.VALUE = logical(1),
    inherits = TRUE
  )
]

if (length(missing_functions) > 0) {
  
  stop(
    "Required helper functions are missing:\n",
    paste(
      missing_functions,
      collapse = "\n"
    )
  )
}


# ============================================================
# EXTRACT SEQUENCES
# ============================================================


# ------------------------------------------------------------
# 3. Extract all sRNA and target sequences
# ------------------------------------------------------------

sequence_results <- vector(
  mode = "list",
  length = nrow(interaction_sequences)
)

for (
  row_index in seq_len(
    nrow(interaction_sequences)
  )
) {
  
  current_row <- interaction_sequences[
    row_index,
  ]
  
  srna_sequence <- extract_genomic_sequence(
    
    genome_sequence =
      genome_sequence,
    
    start_position =
      current_row$srna_start,
    
    end_position =
      current_row$srna_end,
    
    strand =
      current_row$srna_strand
  )
  
  target_result <- extract_target_window(
    
    genome_sequence =
      genome_sequence,
    
    gene_start =
      current_row$target_genomic_start,
    
    gene_end =
      current_row$target_genomic_end,
    
    strand =
      current_row$target_strand,
    
    upstream_length =
      upstream_nt,
    
    downstream_length =
      downstream_nt
  )
  
  sequence_results[[row_index]] <- list(
    
    srna_sequence =
      as.character(
        srna_sequence
      ),
    
    target_sequence =
      as.character(
        target_result$sequence
      ),
    
    target_window_start =
      target_result$window_start,
    
    target_window_end =
      target_result$window_end,
    
    target_translation_start =
      target_result$translation_start
  )
}


# ------------------------------------------------------------
# 4. Add extracted sequences
# ------------------------------------------------------------

interaction_sequences$srna_sequence <-
  vapply(
    sequence_results,
    function(x) {
      as.character(x$srna_sequence)
    },
    FUN.VALUE = character(1)
  )

interaction_sequences$target_sequence <-
  vapply(
    sequence_results,
    function(x) {
      as.character(x$target_sequence)
    },
    FUN.VALUE = character(1)
  )

interaction_sequences$target_window_start <-
  vapply(
    sequence_results,
    function(x) {
      
      if (
        is.null(x$target_window_start) ||
        length(x$target_window_start) == 0 ||
        is.na(x$target_window_start)
      ) {
        return(NA_real_)
      }
      
      as.numeric(
        x$target_window_start
      )
    },
    FUN.VALUE = numeric(1)
  )

interaction_sequences$target_window_end <-
  vapply(
    sequence_results,
    function(x) {
      
      if (
        is.null(x$target_window_end) ||
        length(x$target_window_end) == 0 ||
        is.na(x$target_window_end)
      ) {
        return(NA_real_)
      }
      
      as.numeric(
        x$target_window_end
      )
    },
    FUN.VALUE = numeric(1)
  )

interaction_sequences$target_translation_start <-
  vapply(
    sequence_results,
    function(x) {
      
      if (
        is.null(x$target_translation_start) ||
        length(x$target_translation_start) == 0 ||
        is.na(x$target_translation_start)
      ) {
        return(NA_real_)
      }
      
      as.numeric(
        x$target_translation_start
      )
    },
    FUN.VALUE = numeric(1)
  )


# ------------------------------------------------------------
# 5. Calculate sequence metrics
# ------------------------------------------------------------

interaction_sequences <- interaction_sequences |>
  dplyr::mutate(
    
    srna_length =
      nchar(
        .data$srna_sequence
      ),
    
    target_sequence_length =
      nchar(
        .data$target_sequence
      ),
    
    srna_GC_percent =
      vapply(
        .data$srna_sequence,
        calculate_gc_percent,
        FUN.VALUE = numeric(1)
      ),
    
    target_GC_percent =
      vapply(
        .data$target_sequence,
        calculate_gc_percent,
        FUN.VALUE = numeric(1)
      ),
    
    srna_coordinates_available =
      !is.na(
        .data$srna_start
      ) &
      !is.na(
        .data$srna_end
      ) &
      .data$srna_strand %in%
      c(
        "+",
        "-"
      ),
    
    target_coordinates_available =
      !is.na(
        .data$target_genomic_start
      ) &
      !is.na(
        .data$target_genomic_end
      ) &
      .data$target_strand %in%
      c(
        "+",
        "-"
      ),
    
    srna_sequence_valid =
      !is.na(
        .data$srna_sequence
      ) &
      .data$srna_length >=
      minimum_srna_length &
      stringr::str_detect(
        .data$srna_sequence,
        "^[ACGTNacgtn]+$"
      ),
    
    target_sequence_valid =
      !is.na(
        .data$target_sequence
      ) &
      .data$target_sequence_length >=
      minimum_target_length &
      stringr::str_detect(
        .data$target_sequence,
        "^[ACGTNacgtn]+$"
      ),
    
    interaction_sequence_status =
      dplyr::case_when(
        
        .data$srna_sequence_valid &
          .data$target_sequence_valid ~
          "Ready for structural prediction",
        
        !.data$srna_coordinates_available ~
          "sRNA coordinates unresolved",
        
        !.data$target_coordinates_available ~
          "Target coordinates unresolved",
        
        !.data$srna_sequence_valid ~
          "Invalid or short sRNA sequence",
        
        !.data$target_sequence_valid ~
          "Invalid or short target sequence",
        
        TRUE ~
          "Sequence review required"
      )
  )


# ------------------------------------------------------------
# 6. Print sequence extraction status
# ------------------------------------------------------------

sequence_status_check <- interaction_sequences |>
  dplyr::transmute(
    
    interaction_id =
      .data$interaction_id,
    
    regulatory_RNA =
      .data$curated_source_label,
    
    target_gene =
      .data$final_gene_symbol,
    
    srna_length =
      .data$srna_length,
    
    target_length =
      .data$target_sequence_length,
    
    srna_GC =
      round(
        .data$srna_GC_percent,
        2
      ),
    
    target_GC =
      round(
        .data$target_GC_percent,
        2
      ),
    
    status =
      .data$interaction_sequence_status
  )

print(
  sequence_status_check,
  n = Inf
)

cat(
  "\nReady interactions:",
  sum(
    interaction_sequences$
      interaction_sequence_status ==
      "Ready for structural prediction",
    na.rm = TRUE
  ),
  "of",
  nrow(interaction_sequences),
  "\n"
)


# ============================================================
# WRITE FASTA FILES
# ============================================================


# ------------------------------------------------------------
# 7. Write individual and pair FASTA files
# ------------------------------------------------------------

fasta_manifest <- vector(
  mode = "list",
  length = nrow(interaction_sequences)
)

for (
  row_index in seq_len(
    nrow(interaction_sequences)
  )
) {
  
  current_row <- interaction_sequences[
    row_index,
  ]
  
  source_label_safe <- safe_filename(
    current_row$curated_source_label
  )
  
  target_label_safe <- safe_filename(
    current_row$final_gene_symbol
  )
  
  if (
    is.na(source_label_safe) ||
    source_label_safe == ""
  ) {
    
    source_label_safe <-
      current_row$interaction_id
  }
  
  if (
    is.na(target_label_safe) ||
    target_label_safe == ""
  ) {
    
    target_label_safe <-
      safe_filename(
        current_row$validated_gene_id
      )
  }
  
  pair_prefix <- paste0(
    current_row$interaction_id,
    "_",
    source_label_safe,
    "_to_",
    target_label_safe
  )
  
  srna_fasta_file <- file.path(
    srna_folder,
    paste0(
      pair_prefix,
      "_sRNA.fa"
    )
  )
  
  target_fasta_file <- file.path(
    target_folder,
    paste0(
      pair_prefix,
      "_target.fa"
    )
  )
  
  pair_fasta_file <- file.path(
    pair_folder,
    paste0(
      pair_prefix,
      "_pair.fa"
    )
  )
  
  if (
    current_row$interaction_sequence_status ==
    "Ready for structural prediction"
  ) {
    
    srna_name <- paste0(
      current_row$interaction_id,
      "|sRNA|",
      current_row$curated_source_label,
      "|",
      current_row$srna_start,
      "-",
      current_row$srna_end,
      "|",
      current_row$srna_strand
    )
    
    target_name <- paste0(
      current_row$interaction_id,
      "|target|",
      current_row$final_gene_symbol,
      "|",
      current_row$target_window_start,
      "-",
      current_row$target_window_end,
      "|",
      current_row$target_strand
    )
    
    write_single_fasta(
      
      sequence_value =
        current_row$srna_sequence,
      
      sequence_name =
        srna_name,
      
      output_file =
        srna_fasta_file
    )
    
    write_single_fasta(
      
      sequence_value =
        current_row$target_sequence,
      
      sequence_name =
        target_name,
      
      output_file =
        target_fasta_file
    )
    
    pair_set <- Biostrings::DNAStringSet(
      c(
        current_row$srna_sequence,
        current_row$target_sequence
      )
    )
    
    names(pair_set) <- c(
      srna_name,
      target_name
    )
    
    Biostrings::writeXStringSet(
      
      pair_set,
      
      filepath =
        pair_fasta_file,
      
      format =
        "fasta"
    )
    
  } else {
    
    srna_fasta_file <- NA_character_
    target_fasta_file <- NA_character_
    pair_fasta_file <- NA_character_
  }
  
  fasta_manifest[[row_index]] <- tibble::tibble(
    
    interaction_id =
      current_row$interaction_id,
    
    regulatory_RNA =
      current_row$curated_source_label,
    
    target_gene =
      current_row$final_gene_symbol,
    
    status =
      current_row$interaction_sequence_status,
    
    srna_fasta =
      srna_fasta_file,
    
    target_fasta =
      target_fasta_file,
    
    pair_fasta =
      pair_fasta_file
  )
}

fasta_manifest <- dplyr::bind_rows(
  fasta_manifest
)


# ------------------------------------------------------------
# 8. Write combined FASTA files
# ------------------------------------------------------------

ready_sequences <- interaction_sequences |>
  dplyr::filter(
    .data$interaction_sequence_status ==
      "Ready for structural prediction"
  )

if (nrow(ready_sequences) > 0) {
  
  combined_srna_set <- Biostrings::DNAStringSet(
    ready_sequences$srna_sequence
  )
  
  names(combined_srna_set) <- paste0(
    ready_sequences$interaction_id,
    "|",
    ready_sequences$curated_source_label
  )
  
  Biostrings::writeXStringSet(
    
    combined_srna_set,
    
    filepath = file.path(
      sequence_folder,
      "All_interaction_sRNAs.fa"
    ),
    
    format = "fasta"
  )
  
  combined_target_set <- Biostrings::DNAStringSet(
    ready_sequences$target_sequence
  )
  
  names(combined_target_set) <- paste0(
    ready_sequences$interaction_id,
    "|",
    ready_sequences$final_gene_symbol
  )
  
  Biostrings::writeXStringSet(
    
    combined_target_set,
    
    filepath = file.path(
      sequence_folder,
      "All_target_interaction_windows.fa"
    ),
    
    format = "fasta"
  )
}


# ============================================================
# PREPARE INTARNA COMMANDS
# ============================================================


# ------------------------------------------------------------
# 9. Build IntaRNA command table
# ------------------------------------------------------------

intarna_command_table <- fasta_manifest |>
  dplyr::filter(
    .data$status ==
      "Ready for structural prediction"
  ) |>
  dplyr::mutate(
    
    result_file =
      file.path(
        intarna_result_folder,
        paste0(
          .data$interaction_id,
          "_IntaRNA.csv"
        )
      ),
    
    command =
      paste(
        
        "IntaRNA",
        
        paste0(
          "--query=\"",
          .data$srna_fasta,
          "\""
        ),
        
        paste0(
          "--target=\"",
          .data$target_fasta,
          "\""
        ),
        
        "--outMode=C",
        
        "--outCsvCols=id1,id2,start1,end1,start2,end2,E",
        
        paste0(
          "--out=\"",
          .data$result_file,
          "\""
        )
      )
  )


# ------------------------------------------------------------
# 10. Save IntaRNA commands
# ------------------------------------------------------------

if (nrow(intarna_command_table) > 0) {
  
  writeLines(
    
    c(
      "#!/usr/bin/env bash",
      "",
      "# IntaRNA commands generated by Script 18",
      "",
      intarna_command_table$command
    ),
    
    con = file.path(
      intarna_folder,
      "run_all_IntaRNA.sh"
    )
  )
  
  writeLines(
    
    intarna_command_table$command,
    
    con = file.path(
      intarna_folder,
      "IntaRNA_commands.txt"
    )
  )
}


# ============================================================
# OPTIONAL INTARNA EXECUTION
# ============================================================


# ------------------------------------------------------------
# 11. Detect IntaRNA
# ------------------------------------------------------------

intarna_executable <- Sys.which(
  "IntaRNA"
)

intarna_available <- (
  nchar(
    intarna_executable
  ) > 0
)

cat(
  "\nIntaRNA executable detected:",
  intarna_available,
  "\n"
)


# ------------------------------------------------------------
# 12. Run IntaRNA when available
# ------------------------------------------------------------

intarna_run_report <- tibble::tibble()

if (
  run_IntaRNA_if_available &&
  intarna_available &&
  nrow(intarna_command_table) > 0
) {
  
  intarna_run_list <- vector(
    mode = "list",
    length = nrow(
      intarna_command_table
    )
  )
  
  for (
    row_index in seq_len(
      nrow(intarna_command_table)
    )
  ) {
    
    current_command <-
      intarna_command_table[
        row_index,
      ]
    
    command_arguments <- c(
      
      paste0(
        "--query=",
        current_command$srna_fasta
      ),
      
      paste0(
        "--target=",
        current_command$target_fasta
      ),
      
      "--outMode=C",
      
      "--outCsvCols=id1,id2,start1,end1,start2,end2,E",
      
      paste0(
        "--out=",
        current_command$result_file
      )
    )
    
    command_output <- tryCatch(
      
      system2(
        command =
          intarna_executable,
        
        args =
          command_arguments,
        
        stdout = TRUE,
        
        stderr = TRUE
      ),
      
      error = function(e) {
        
        paste0(
          "ERROR: ",
          conditionMessage(e)
        )
      }
    )
    
    result_created <- file.exists(
      current_command$result_file
    )
    
    intarna_run_list[[row_index]] <-
      tibble::tibble(
        
        interaction_id =
          current_command$interaction_id,
        
        regulatory_RNA =
          current_command$regulatory_RNA,
        
        target_gene =
          current_command$target_gene,
        
        result_file =
          current_command$result_file,
        
        result_created =
          result_created,
        
        command_output =
          paste(
            command_output,
            collapse = " | "
          )
      )
  }
  
  intarna_run_report <- dplyr::bind_rows(
    intarna_run_list
  )
  
} else {
  
  intarna_run_report <- intarna_command_table |>
    dplyr::transmute(
      
      interaction_id =
        .data$interaction_id,
      
      regulatory_RNA =
        .data$regulatory_RNA,
      
      target_gene =
        .data$target_gene,
      
      result_file =
        .data$result_file,
      
      result_created =
        FALSE,
      
      command_output =
        dplyr::case_when(
          
          !intarna_available ~
            "IntaRNA not installed; FASTA files and commands prepared",
          
          !run_IntaRNA_if_available ~
            "Automatic IntaRNA execution disabled",
          
          TRUE ~
            "No interaction available for execution"
        )
    )
}


# ============================================================
# IMPORT AVAILABLE INTARNA RESULTS
# ============================================================


# ------------------------------------------------------------
# 13. Import completed IntaRNA outputs
# ------------------------------------------------------------

intarna_result_list <- list()

if (nrow(intarna_command_table) > 0) {
  
  for (
    row_index in seq_len(
      nrow(intarna_command_table)
    )
  ) {
    
    current_row <- intarna_command_table[
      row_index,
    ]
    
    if (
      file.exists(
        current_row$result_file
      ) &&
      file.info(
        current_row$result_file
      )$size > 0
    ) {
      
      result_data <- tryCatch(
        
        readr::read_delim(
          
          current_row$result_file,
          
          delim = ";",
          
          show_col_types = FALSE,
          
          progress = FALSE
        ),
        
        error = function(e) {
          
          tryCatch(
            
            readr::read_csv(
              
              current_row$result_file,
              
              show_col_types = FALSE
            ),
            
            error = function(e2) {
              NULL
            }
          )
        }
      )
      
      if (
        !is.null(result_data) &&
        nrow(result_data) > 0
      ) {
        
        result_data$interaction_id <-
          current_row$interaction_id
        
        result_data$regulatory_RNA <-
          current_row$regulatory_RNA
        
        result_data$target_gene <-
          current_row$target_gene
        
        intarna_result_list[
          [
            length(
              intarna_result_list
            ) + 1
          ]
        ] <- list(
          result_data
        )
      }
    }
  }
}

if (length(intarna_result_list) > 0) {
  
  combined_intarna_results <-
    dplyr::bind_rows(
      intarna_result_list
    )
  
} else {
  
  combined_intarna_results <-
    tibble::tibble(
      
      interaction_id = character(),
      regulatory_RNA = character(),
      target_gene = character(),
      
      id1 = character(),
      id2 = character(),
      
      start1 = numeric(),
      end1 = numeric(),
      
      start2 = numeric(),
      end2 = numeric(),
      
      E = numeric()
    )
}


# ============================================================
# SAVE TABLES
# ============================================================


# ------------------------------------------------------------
# 14. Save sequence metadata
# ------------------------------------------------------------

readr::write_csv(
  
  interaction_sequences,
  
  file.path(
    table_folder,
    "Interaction_sequence_metadata.csv"
  )
)

readr::write_csv(
  
  fasta_manifest,
  
  file.path(
    table_folder,
    "Interaction_FASTA_manifest.csv"
  )
)

readr::write_csv(
  
  intarna_command_table,
  
  file.path(
    table_folder,
    "IntaRNA_command_manifest.csv"
  )
)

readr::write_csv(
  
  intarna_run_report,
  
  file.path(
    table_folder,
    "IntaRNA_run_report.csv"
  )
)

readr::write_csv(
  
  combined_intarna_results,
  
  file.path(
    table_folder,
    "Combined_IntaRNA_predictions.csv"
  )
)


# ------------------------------------------------------------
# 15. Save sequence audit
# ------------------------------------------------------------

sequence_audit <- interaction_sequences |>
  dplyr::transmute(
    
    interaction_id =
      .data$interaction_id,
    
    regulatory_RNA =
      .data$curated_source_label,
    
    target_gene =
      .data$final_gene_symbol,
    
    srna_coordinates =
      paste0(
        .data$srna_start,
        "-",
        .data$srna_end,
        ":",
        .data$srna_strand
      ),
    
    srna_length =
      .data$srna_length,
    
    srna_GC_percent =
      round(
        .data$srna_GC_percent,
        2
      ),
    
    target_coordinates =
      paste0(
        .data$target_window_start,
        "-",
        .data$target_window_end,
        ":",
        .data$target_strand
      ),
    
    target_sequence_length =
      .data$target_sequence_length,
    
    target_GC_percent =
      round(
        .data$target_GC_percent,
        2
      ),
    
    status =
      .data$interaction_sequence_status
  )

readr::write_csv(
  
  sequence_audit,
  
  file.path(
    audit_folder,
    "Sequence_extraction_audit.csv"
  )
)


# ------------------------------------------------------------
# 16. Build Script 18 summary
# ------------------------------------------------------------

script18_summary <- tibble::tibble(
  
  metric = c(
    
    "Validated interactions read",
    
    "Interactions with sRNA coordinates",
    
    "Interactions with target coordinates",
    
    "Interactions ready for structural prediction",
    
    "Interactions requiring sequence review",
    
    "Individual sRNA FASTA files",
    
    "Individual target FASTA files",
    
    "IntaRNA executable available",
    
    "IntaRNA predictions imported"
  ),
  
  value = c(
    
    nrow(
      interaction_sequences
    ),
    
    sum(
      interaction_sequences$
        srna_coordinates_available,
      na.rm = TRUE
    ),
    
    sum(
      interaction_sequences$
        target_coordinates_available,
      na.rm = TRUE
    ),
    
    sum(
      interaction_sequences$
        interaction_sequence_status ==
        "Ready for structural prediction",
      na.rm = TRUE
    ),
    
    sum(
      interaction_sequences$
        interaction_sequence_status !=
        "Ready for structural prediction",
      na.rm = TRUE
    ),
    
    sum(
      !is.na(
        fasta_manifest$srna_fasta
      )
    ),
    
    sum(
      !is.na(
        fasta_manifest$target_fasta
      )
    ),
    
    as.integer(
      intarna_available
    ),
    
    dplyr::n_distinct(
      combined_intarna_results$
        interaction_id
    )
  )
)

readr::write_csv(
  
  script18_summary,
  
  file.path(
    audit_folder,
    "Script18_summary.csv"
  )
)


# ------------------------------------------------------------
# 17. Verify outputs
# ------------------------------------------------------------

expected_outputs <- c(
  
  file.path(
    table_folder,
    "Interaction_sequence_metadata.csv"
  ),
  
  file.path(
    table_folder,
    "Interaction_FASTA_manifest.csv"
  ),
  
  file.path(
    table_folder,
    "IntaRNA_command_manifest.csv"
  ),
  
  file.path(
    table_folder,
    "IntaRNA_run_report.csv"
  ),
  
  file.path(
    table_folder,
    "Combined_IntaRNA_predictions.csv"
  ),
  
  file.path(
    audit_folder,
    "Sequence_extraction_audit.csv"
  ),
  
  file.path(
    audit_folder,
    "Script18_summary.csv"
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
      )$size / 1024,
      2
    )
)

print(
  output_verification,
  n = Inf
)


# ------------------------------------------------------------
# 18. Final console report
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 18 COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================\n\n"
)

cat(
  "Sequence extraction summary:\n"
)

print(
  script18_summary,
  n = Inf
)

cat(
  "\nSequence audit:\n"
)

print(
  sequence_audit,
  n = Inf
)

cat(
  "\nFASTA files:\n",
  sequence_folder,
  "\n\n"
)

cat(
  "IntaRNA commands:\n",
  file.path(
    intarna_folder,
    "IntaRNA_commands.txt"
  ),
  "\n\n"
)

cat(
  "Combined IntaRNA results:\n",
  file.path(
    table_folder,
    "Combined_IntaRNA_predictions.csv"
  ),
  "\n"
)
# ============================================================
# FIX — Import available IntaRNA results
# ============================================================

intarna_result_list <- list()

if (
  exists("intarna_command_table") &&
  nrow(intarna_command_table) > 0
) {
  
  for (
    row_index in seq_len(
      nrow(intarna_command_table)
    )
  ) {
    
    current_row <- intarna_command_table[
      row_index,
      ,
      drop = FALSE
    ]
    
    current_result_file <- as.character(
      current_row$result_file[[1]]
    )
    
    if (
      !is.na(current_result_file) &&
      file.exists(current_result_file) &&
      file.info(current_result_file)$size > 0
    ) {
      
      result_data <- tryCatch(
        
        readr::read_delim(
          current_result_file,
          delim = ";",
          show_col_types = FALSE,
          progress = FALSE,
          trim_ws = TRUE
        ),
        
        error = function(e1) {
          
          tryCatch(
            
            readr::read_csv(
              current_result_file,
              show_col_types = FALSE,
              progress = FALSE,
              trim_ws = TRUE
            ),
            
            error = function(e2) {
              
              tryCatch(
                
                readr::read_tsv(
                  current_result_file,
                  show_col_types = FALSE,
                  progress = FALSE,
                  trim_ws = TRUE
                ),
                
                error = function(e3) {
                  NULL
                }
              )
            }
          )
        }
      )
      
      if (
        !is.null(result_data) &&
        nrow(result_data) > 0
      ) {
        
        result_data$interaction_id <-
          as.character(
            current_row$interaction_id[[1]]
          )
        
        result_data$regulatory_RNA <-
          as.character(
            current_row$regulatory_RNA[[1]]
          )
        
        result_data$target_gene <-
          as.character(
            current_row$target_gene[[1]]
          )
        
        # Correct list append syntax
        intarna_result_list[[
          length(intarna_result_list) + 1L
        ]] <- result_data
      }
    }
  }
}


# ------------------------------------------------------------
# Combine imported IntaRNA results
# ------------------------------------------------------------

if (length(intarna_result_list) > 0) {
  
  combined_intarna_results <-
    dplyr::bind_rows(
      intarna_result_list
    )
  
} else {
  
  combined_intarna_results <-
    tibble::tibble(
      
      interaction_id = character(),
      regulatory_RNA = character(),
      target_gene = character(),
      
      id1 = character(),
      id2 = character(),
      
      start1 = numeric(),
      end1 = numeric(),
      
      start2 = numeric(),
      end2 = numeric(),
      
      E = numeric()
    )
}


# ------------------------------------------------------------
# Report imported results
# ------------------------------------------------------------

cat(
  "\nIntaRNA result files imported:",
  length(intarna_result_list),
  "\n"
)

cat(
  "Imported prediction rows:",
  nrow(combined_intarna_results),
  "\n"
)

print(
  combined_intarna_results,
  n = Inf
)
cat(
  "Ready interactions:",
  sum(
    interaction_sequences$interaction_sequence_status ==
      "Ready for structural prediction",
    na.rm = TRUE
  ),
  "of",
  nrow(interaction_sequences),
  "\n"
)

cat(
  "FASTA files created:",
  sum(!is.na(fasta_manifest$pair_fasta)),
  "\n"
)

cat(
  "IntaRNA available:",
  intarna_available,
  "\n"
)