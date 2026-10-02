# ============================================================
# SCRIPT 17B
# Parent-gene / genomic-context audit for the two
# putative 3'UTR-associated RNA candidates
#
# Goals:
# 1. Identify RS04120 and RS05585 in the JKD6008 Master KB.
# 2. Retrieve parent-gene coordinates, strand and product.
# 3. Compare CLASH-fragment coordinates with parent-gene 3' ends.
# 4. Re-examine Term-seq around:
#       A) CLASH candidate boundary
#       B) parent-gene 3' boundary
#
# IMPORTANT:
# No claim of independent processing is made.
# ============================================================


# ============================================================
# 01. PATHS
# ============================================================

project <- "D:/Bac-sRNA"

master_file <- file.path(
  project,
  "13_CLASH_network",
  "tables",
  "JKD6008_master_gene_knowledgebase_CLASH_updated.csv"
)

term_dir <- file.path(
  project,
  "07_processed_data",
  "downloaded_files",
  "GSE158830"
)

output_dir <- file.path(
  project,
  "17B_parent_gene_boundary_audit"
)

table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 02. PACKAGES
# ============================================================

packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "tibble",
  "purrr",
  "ggplot2"
)

missing <- setdiff(
  packages,
  rownames(installed.packages())
)

if (length(missing) > 0) {
  install.packages(missing)
}

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)
library(ggplot2)


# ============================================================
# 03. READ MASTER KB
# ============================================================

if (!file.exists(master_file)) {
  stop(
    paste(
      "Master KB not found:",
      master_file
    )
  )
}

kb <- read_csv(
  master_file,
  show_col_types = FALSE
)

cat(
  "\nMaster KB rows:",
  nrow(kb),
  "\n"
)


# ============================================================
# 04. CANDIDATE CLASH FRAGMENTS
# ============================================================

candidates <- tribble(
  
  ~candidate,
  ~candidate_start,
  ~candidate_end,
  ~candidate_strand,
  ~parent_RS_query,
  
  "3UTR-RS04120",
  811819,
  811847,
  "+",
  "RS04120",
  
  "3UTR-RS05585",
  1084557,
  1084613,
  "-",
  "RS05585"
)


candidates <- candidates %>%
  
  mutate(
    
    candidate_3prime =
      if_else(
        candidate_strand == "+",
        candidate_end,
        candidate_start
      )
  )


# ============================================================
# 05. NORMALIZATION
# ============================================================

norm <- function(x) {
  
  x %>%
    as.character() %>%
    str_to_upper() %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("-", "_")
}


# ============================================================
# 06. SEARCH MASTER KB FOR PARENT GENES
#
# Search across identifiers, not only one field.
# ============================================================

search_columns <- intersect(
  c(
    "gene_id",
    "current_locus_tag",
    "gene_symbol",
    "gene_synonym",
    "master_label",
    "display_label",
    "protein_id"
  ),
  names(kb)
)


if (length(search_columns) == 0) {
  stop("No usable identifier columns in Master KB.")
}


kb_search <- kb %>%
  
  mutate(
    master_row_id = row_number()
  ) %>%
  
  pivot_longer(
    cols = all_of(search_columns),
    names_to = "matched_field",
    values_to = "identifier"
  ) %>%
  
  mutate(
    identifier_norm = norm(identifier)
  )


parent_hits <- map_dfr(
  
  seq_len(nrow(candidates)),
  
  function(i) {
    
    q <- candidates$parent_RS_query[i]
    
    hits <- kb_search %>%
      
      filter(
        str_detect(
          identifier_norm,
          fixed(norm(q))
        )
      ) %>%
      
      mutate(
        candidate =
          candidates$candidate[i],
        
        parent_RS_query =
          q
      )
    
    hits
  }
)


write_csv(
  parent_hits,
  file.path(
    table_dir,
    "Parent_Gene_All_Identifier_Hits.csv"
  )
)


# ============================================================
# 07. COLLAPSE TO UNIQUE MASTER KB ROWS
# ============================================================

parent_rows <- parent_hits %>%
  
  distinct(
    candidate,
    parent_RS_query,
    master_row_id
  ) %>%
  
  left_join(
    kb %>%
      mutate(
        master_row_id = row_number()
      ),
    by = "master_row_id"
  )


write_csv(
  parent_rows,
  file.path(
    table_dir,
    "Parent_Gene_Candidate_MasterRows.csv"
  )
)


# ============================================================
# 08. PRINT PARENT GENE CANDIDATES BEFORE SELECTION
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "PARENT GENE MASTER-KB CANDIDATES\n"
)

cat(
  "============================================\n"
)


parent_display_columns <- intersect(
  
  c(
    "candidate",
    "parent_RS_query",
    "gene_id",
    "current_locus_tag",
    "gene_symbol",
    "gene_synonym",
    "protein_id",
    "product",
    "feature_type",
    "genomic_start",
    "genomic_end",
    "strand"
  ),
  
  names(parent_rows)
)


print(
  parent_rows %>%
    select(
      all_of(parent_display_columns)
    ),
  n = Inf,
  width = Inf
)


# ============================================================
# 09. SCORE PARENT MATCH QUALITY
#
# Prefer exact locus-tag / identifier match.
# ============================================================

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}


for (nm in c(
  "gene_id",
  "current_locus_tag",
  "gene_symbol",
  "gene_synonym",
  "master_label",
  "display_label"
)) {
  
  if (!nm %in% names(parent_rows)) {
    parent_rows[[nm]] <- ""
  }
  
  parent_rows[[nm]] <-
    safe_chr(
      parent_rows[[nm]]
    )
}


parent_rows <- parent_rows %>%
  
  rowwise() %>%
  
  mutate(
    
    exact_identifier_match =
      any(
        norm(
          c(
            gene_id,
            current_locus_tag,
            gene_symbol,
            gene_synonym,
            master_label,
            display_label
          )
        ) ==
          norm(parent_RS_query)
      ),
    
    contains_identifier_match =
      any(
        str_detect(
          norm(
            c(
              gene_id,
              current_locus_tag,
              gene_symbol,
              gene_synonym,
              master_label,
              display_label
            )
          ),
          fixed(
            norm(parent_RS_query)
          )
        )
      ),
    
    parent_match_score =
      case_when(
        exact_identifier_match ~ 2,
        contains_identifier_match ~ 1,
        TRUE ~ 0
      )
  ) %>%
  
  ungroup()


# ============================================================
# 10. SELECT BEST PARENT ROW
#
# Do not silently resolve ties.
# ============================================================

best_parent_candidates <- parent_rows %>%
  
  group_by(candidate) %>%
  
  filter(
    parent_match_score ==
      max(
        parent_match_score,
        na.rm = TRUE
      )
  ) %>%
  
  ungroup()


ambiguity_summary <- best_parent_candidates %>%
  
  count(
    candidate,
    name = "number_of_best_parent_rows"
  )


write_csv(
  ambiguity_summary,
  file.path(
    table_dir,
    "Parent_Gene_Ambiguity_Summary.csv"
  )
)


# ============================================================
# 11. KEEP ONE ROW ONLY IF UNIQUE BEST MATCH
# ============================================================

unique_parent <- best_parent_candidates %>%
  
  group_by(candidate) %>%
  
  mutate(
    number_of_best_parent_rows = n()
  ) %>%
  
  filter(
    number_of_best_parent_rows == 1
  ) %>%
  
  ungroup()


if (n_distinct(unique_parent$candidate) < 2) {
  
  warning(
    paste(
      "At least one parent gene does not have",
      "a unique best Master-KB match.",
      "Review Parent_Gene_Candidate_MasterRows.csv."
    )
  )
}


# ============================================================
# 12. BUILD GENOMIC CONTEXT
# ============================================================

context_columns <- intersect(
  
  c(
    "candidate",
    "gene_id",
    "current_locus_tag",
    "gene_symbol",
    "protein_id",
    "product",
    "feature_type",
    "genomic_start",
    "genomic_end",
    "strand"
  ),
  
  names(unique_parent)
)


parent_context <- unique_parent %>%
  
  select(
    all_of(context_columns)
  )


if (!all(
  c(
    "genomic_start",
    "genomic_end",
    "strand"
  ) %in% names(parent_context)
)) {
  
  stop(
    paste(
      "Master KB lacks required genomic coordinate",
      "or strand columns."
    )
  )
}


parent_context <- parent_context %>%
  
  mutate(
    
    genomic_start =
      as.numeric(genomic_start),
    
    genomic_end =
      as.numeric(genomic_end),
    
    strand =
      as.character(strand),
    
    parent_gene_3prime =
      if_else(
        strand == "+",
        genomic_end,
        genomic_start
      )
  )


# ============================================================
# 13. JOIN CLASH CANDIDATE COORDINATES
# ============================================================

context <- candidates %>%
  
  left_join(
    parent_context,
    by = "candidate"
  ) %>%
  
  mutate(
    
    candidate_vs_parent_strand =
      case_when(
        is.na(strand) ~
          "Parent unresolved",
        
        candidate_strand == strand ~
          "Same strand",
        
        TRUE ~
          "Opposite strand"
      ),
    
    distance_candidate_3prime_to_parent_3prime =
      candidate_3prime -
      parent_gene_3prime,
    
    absolute_distance_to_parent_3prime =
      abs(
        distance_candidate_3prime_to_parent_3prime
      )
  )


write_csv(
  context,
  file.path(
    table_dir,
    "Candidate_Parent_Gene_Genomic_Context.csv"
  )
)


# ============================================================
# 14. READ TERM-SEQ FILES
# ============================================================

term_files <- list.files(
  
  term_dir,
  
  pattern =
    "termseq.*(plus|minus).*\\.gr\\.gz$",
  
  full.names = TRUE,
  
  ignore.case = TRUE
)


term_meta <- tibble(
  
  file = term_files,
  
  filename =
    basename(term_files)
) %>%
  
  mutate(
    
    replicate =
      case_when(
        str_detect(
          filename,
          "GSM4811618"
        ) ~ "rep1",
        
        str_detect(
          filename,
          "GSM4811619"
        ) ~ "rep2",
        
        TRUE ~ "unknown"
      ),
    
    strand =
      case_when(
        str_detect(
          filename,
          "plus_strand"
        ) ~ "+",
        
        str_detect(
          filename,
          "minus_strand"
        ) ~ "-",
        
        TRUE ~ "?"
      )
  )


# ============================================================
# 15. ROBUST TERM-SEQ READER
# ============================================================

read_term <- function(
    file,
    replicate,
    strand
) {
  
  x <- read_table(
    file,
    col_names = FALSE,
    comment = "#",
    show_col_types = FALSE,
    progress = FALSE
  )
  
  
  if (ncol(x) == 2) {
    
    names(x) <-
      c(
        "position",
        "signal"
      )
    
    x <- x %>%
      mutate(
        chromosome = NA_character_
      ) %>%
      select(
        chromosome,
        position,
        signal
      )
    
  } else if (ncol(x) >= 3) {
    
    x <- x[, 1:3]
    
    names(x) <-
      c(
        "chromosome",
        "position",
        "signal"
      )
    
  } else {
    
    stop(
      paste(
        "Unexpected Term-seq format:",
        file
      )
    )
  }
  
  
  x %>%
    
    mutate(
      
      position =
        suppressWarnings(
          as.numeric(position)
        ),
      
      signal =
        suppressWarnings(
          as.numeric(signal)
        ),
      
      replicate =
        replicate,
      
      strand =
        strand,
      
      source_file =
        basename(file)
    ) %>%
    
    filter(
      !is.na(position),
      !is.na(signal)
    )
}


term_data <- pmap_dfr(
  
  term_meta,
  
  function(
    file,
    filename,
    replicate,
    strand
  ) {
    
    read_term(
      file,
      replicate,
      strand
    )
  }
)


# ============================================================
# 16. BOUNDARIES TO TEST
#
# Test both:
# A. CLASH fragment 3' coordinate
# B. parent-gene 3' coordinate
# ============================================================

boundary_table <- context %>%
  
  transmute(
    
    candidate,
    
    candidate_strand,
    
    parent_strand = strand,
    
    candidate_3prime,
    
    parent_gene_3prime,
    
    absolute_distance_to_parent_3prime
  ) %>%
  
  pivot_longer(
    
    cols =
      c(
        candidate_3prime,
        parent_gene_3prime
      ),
    
    names_to =
      "boundary_type",
    
    values_to =
      "boundary_position"
  ) %>%
  
  mutate(
    
    boundary_type =
      recode(
        boundary_type,
        
        candidate_3prime =
          "CLASH candidate boundary",
        
        parent_gene_3prime =
          "Parent-gene 3′ boundary"
      ),
    
    analysis_strand =
      case_when(
        
        boundary_type ==
          "CLASH candidate boundary" ~
          candidate_strand,
        
        boundary_type ==
          "Parent-gene 3′ boundary" ~
          parent_strand,
        
        TRUE ~
          candidate_strand
      )
  ) %>%
  
  filter(
    !is.na(boundary_position),
    !is.na(analysis_strand)
  )


write_csv(
  boundary_table,
  file.path(
    table_dir,
    "Boundaries_Tested.csv"
  )
)


# ============================================================
# 17. EXTRACT TERM-SEQ WINDOWS
# ============================================================

window_size <- 150


extract_window <- function(
    candidate_name,
    boundary_type_value,
    boundary,
    analysis_strand_value
) {
  
  term_data %>%
    
    filter(
      strand ==
        analysis_strand_value,
      
      position >=
        boundary - window_size,
      
      position <=
        boundary + window_size
    ) %>%
    
    mutate(
      
      candidate =
        candidate_name,
      
      boundary_type =
        boundary_type_value,
      
      boundary_position =
        boundary,
      
      relative_position =
        position - boundary
    )
}


windows <- pmap_dfr(
  
  boundary_table,
  
  function(
    candidate,
    candidate_strand,
    parent_strand,
    boundary_type,
    boundary_position,
    absolute_distance_to_parent_3prime,
    analysis_strand
  ) {
    
    extract_window(
      candidate,
      boundary_type,
      boundary_position,
      analysis_strand
    )
  }
)


# ============================================================
# 18. BOUNDARY METRICS
# ============================================================

boundary_radius <- 5
background_inner <- 20
background_outer <- 100


metrics <- windows %>%
  
  group_by(
    candidate,
    boundary_type,
    boundary_position,
    replicate,
    strand
  ) %>%
  
  summarise(
    
    boundary_max_signal =
      {
        z <- signal[
          abs(relative_position) <=
            boundary_radius
        ]
        
        if (length(z) == 0) 0
        else max(z, na.rm = TRUE)
      },
    
    boundary_total_signal =
      sum(
        signal[
          abs(relative_position) <=
            boundary_radius
        ],
        na.rm = TRUE
      ),
    
    background_median =
      {
        z <- signal[
          abs(relative_position) >=
            background_inner &
            abs(relative_position) <=
            background_outer
        ]
        
        if (length(z) == 0) 0
        else median(z, na.rm = TRUE)
      },
    
    strongest_signal =
      ifelse(
        length(signal) == 0,
        0,
        max(signal, na.rm = TRUE)
      ),
    
    strongest_relative_position =
      {
        if (length(signal) == 0) {
          NA_real_
        } else {
          relative_position[
            which.max(signal)
          ]
        }
      },
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    enrichment_ratio =
      (
        boundary_max_signal + 1
      ) /
      (
        background_median + 1
      ),
    
    strongest_peak_within_5nt =
      !is.na(
        strongest_relative_position
      ) &
      abs(
        strongest_relative_position
      ) <= 5,
    
    replicate_support =
      case_when(
        
        boundary_max_signal > 0 &
          enrichment_ratio >= 5 &
          strongest_peak_within_5nt ~
          "Strong local boundary peak",
        
        boundary_max_signal > 0 &
          enrichment_ratio >= 2 ~
          "Local boundary enrichment",
        
        boundary_max_signal > 0 ~
          "Boundary signal detected",
        
        TRUE ~
          "No boundary signal detected"
      )
  )


write_csv(
  metrics,
  file.path(
    table_dir,
    "Parent_vs_CLASH_TermSeq_PerReplicate.csv"
  )
)


# ============================================================
# 19. REPLICATE SUMMARY
# ============================================================

summary_table <- metrics %>%
  
  group_by(
    candidate,
    boundary_type,
    boundary_position,
    strand
  ) %>%
  
  summarise(
    
    replicates_examined =
      n_distinct(
        replicate
      ),
    
    replicates_with_signal =
      sum(
        boundary_max_signal > 0,
        na.rm = TRUE
      ),
    
    replicates_with_enrichment =
      sum(
        boundary_max_signal > 0 &
          enrichment_ratio >= 2,
        na.rm = TRUE
      ),
    
    replicates_with_strong_peak =
      sum(
        boundary_max_signal > 0 &
          enrichment_ratio >= 5 &
          strongest_peak_within_5nt,
        na.rm = TRUE
      ),
    
    median_boundary_signal =
      median(
        boundary_max_signal,
        na.rm = TRUE
      ),
    
    median_enrichment_ratio =
      median(
        enrichment_ratio,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    boundary_support_class =
      case_when(
        
        replicates_with_strong_peak >= 2 ~
          "Reproducible strong boundary peak",
        
        replicates_with_enrichment >= 2 ~
          "Reproducible local boundary enrichment",
        
        replicates_with_signal >= 2 ~
          "Reproducible signal without enrichment",
        
        replicates_with_signal == 1 ~
          "Single-replicate signal",
        
        TRUE ~
          "No boundary signal detected"
      )
  )


write_csv(
  summary_table,
  file.path(
    table_dir,
    "Parent_vs_CLASH_TermSeq_Summary.csv"
  )
)


# ============================================================
# 20. GENOMIC RELATIONSHIP CLASSIFICATION
# ============================================================

context <- context %>%
  
  mutate(
    
    genomic_relationship =
      case_when(
        
        is.na(parent_gene_3prime) ~
          "Parent gene unresolved",
        
        candidate_vs_parent_strand !=
          "Same strand" ~
          "Strand mismatch",
        
        absolute_distance_to_parent_3prime <= 10 ~
          "Candidate boundary overlaps parent 3′ end",
        
        absolute_distance_to_parent_3prime <= 50 ~
          "Candidate lies near parent 3′ end",
        
        absolute_distance_to_parent_3prime <= 200 ~
          "Candidate lies within extended 3′ region",
        
        TRUE ~
          "Candidate distant from parent 3′ end"
      )
  )


write_csv(
  context,
  file.path(
    table_dir,
    "FINAL_Candidate_Parent_Genomic_Context.csv"
  )
)


# ============================================================
# 21. FINAL INTEGRATED AUDIT
# ============================================================

integrated <- context %>%
  
  select(
    candidate,
    parent_RS_query,
    candidate_start,
    candidate_end,
    candidate_strand,
    candidate_3prime,
    gene_id,
    current_locus_tag,
    gene_symbol,
    protein_id,
    product,
    genomic_start,
    genomic_end,
    strand,
    parent_gene_3prime,
    candidate_vs_parent_strand,
    distance_candidate_3prime_to_parent_3prime,
    absolute_distance_to_parent_3prime,
    genomic_relationship
  ) %>%
  
  left_join(
    
    summary_table %>%
      
      select(
        candidate,
        boundary_type,
        boundary_position,
        boundary_support_class,
        replicates_with_signal,
        replicates_with_enrichment,
        replicates_with_strong_peak,
        median_enrichment_ratio
      ) %>%
      
      pivot_wider(
        
        names_from =
          boundary_type,
        
        values_from =
          c(
            boundary_position,
            boundary_support_class,
            replicates_with_signal,
            replicates_with_enrichment,
            replicates_with_strong_peak,
            median_enrichment_ratio
          ),
        
        names_sep =
          "__"
      ),
    
    by =
      "candidate"
  )


write_csv(
  integrated,
  file.path(
    table_dir,
    "FINAL_ParentGene_TermSeq_Integrated_Audit.csv"
  )
)


# ============================================================
# 22. FIGURE
# ============================================================

plot_data <- windows %>%
  
  filter(
    abs(relative_position) <= 100
  )


p <- ggplot(
  
  plot_data,
  
  aes(
    x = relative_position,
    y = signal
  )
) +
  
  geom_line(
    linewidth = 0.45
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  facet_grid(
    candidate + boundary_type ~ replicate,
    scales = "free_y"
  ) +
  
  labs(
    
    title =
      "Term-seq evidence at candidate and parent-gene 3′ boundaries",
    
    subtitle =
      "Position 0 denotes the boundary being tested",
    
    x =
      "Relative genomic position (nt)",
    
    y =
      "Term-seq 3′-end signal"
  ) +
  
  theme_classic(
    base_size = 9
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    strip.text =
      element_text(
        size = 8
      )
  )


print(p)


ggsave(
  file.path(
    figure_dir,
    "Figure_Parent_vs_CLASH_3prime_Boundaries.pdf"
  ),
  p,
  width = 10,
  height = 8,
  device = cairo_pdf,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_Parent_vs_CLASH_3prime_Boundaries_600dpi.png"
  ),
  p,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white"
)


# ============================================================
# 23. METHODOLOGICAL NOTES
# ============================================================

notes <- c(
  
  "SCRIPT 17B - PARENT-GENE / TERM-SEQ BOUNDARY AUDIT",
  
  "",
  
  paste(
    "This analysis distinguishes the CLASH-fragment",
    "boundary from the annotated parent-gene 3-prime end."
  ),
  
  paste(
    "Parent genes were resolved from the JKD6008 master",
    "knowledgebase using RS04120 and RS05585 identifiers."
  ),
  
  paste(
    "Term-seq was evaluated independently around both",
    "the CLASH candidate boundary and parent-gene",
    "3-prime boundary."
  ),
  
  paste(
    "A Term-seq peak supports a discrete 3-prime end",
    "but does not prove independent RNA processing."
  ),
  
  paste(
    "Absence of a peak under these experimental",
    "conditions does not prove absence of a transcript."
  ),
  
  paste(
    "No biological-priority score was calculated."
  )
)


writeLines(
  notes,
  file.path(
    output_dir,
    "Script17B_Methodological_Notes.txt"
  )
)


# ============================================================
# 24. SESSION INFO
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "sessionInfo.txt"
  )
)


# ============================================================
# 25. FINAL CONSOLE OUTPUT
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 17B COMPLETED\n"
)

cat(
  "============================================\n\n"
)


cat(
  "PARENT-GENE GENOMIC CONTEXT\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  context,
  n = Inf,
  width = Inf
)


cat(
  "\nTERM-SEQ BOUNDARY COMPARISON\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  summary_table,
  n = Inf,
  width = Inf
)


cat(
  "\nFINAL INTEGRATED AUDIT\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  integrated,
  n = Inf,
  width = Inf
)


cat(
  "\n============================================\n"
)

cat(
  "INTERPRETATION RULE\n"
)

cat(
  "============================================\n"
)

cat(
  paste(
    "\nDo not interpret a CLASH-fragment coordinate",
    "as a transcript boundary without independent evidence.",
    "\nTerm-seq supports 3′-end localization only.",
    "\nIndependent processing remains unproven without",
    "additional 5′-boundary/processing evidence.\n"
  )
)


cat(
  "\nOutputs:",
  output_dir,
  "\n"
)


cat(
  "\n============================================\n"
)

cat(
  "END SCRIPT 17B\n"
)

cat(
  "============================================\n"
)
# ============================================================
# SCRIPT 17C-A
# Resolve the true genomic neighborhood of the two
# CLASH-defined 3'UTR-associated candidates
# ============================================================

library(readr)
library(dplyr)
library(tidyr)

project <- "D:/Bac-sRNA"

kb_file <- file.path(
  project,
  "13_CLASH_network",
  "tables",
  "JKD6008_master_gene_knowledgebase_CLASH_updated.csv"
)

kb <- read_csv(kb_file, show_col_types = FALSE)

candidates <- tibble::tribble(
  ~candidate,       ~start,   ~end,     ~strand,
  "3UTR-RS04120",    811819,   811847,   "+",
  "3UTR-RS05585",   1084557,  1084613,   "-"
)

# Ensure coordinates numeric
kb <- kb %>%
  mutate(
    genomic_start = as.numeric(genomic_start),
    genomic_end   = as.numeric(genomic_end)
  )

# Search +/- 5 kb around each candidate
neighborhoods <- lapply(seq_len(nrow(candidates)), function(i) {
  
  cand <- candidates[i, ]
  
  center <- mean(c(cand$start, cand$end))
  
  kb %>%
    filter(
      !is.na(genomic_start),
      !is.na(genomic_end),
      genomic_end >= cand$start - 5000,
      genomic_start <= cand$end + 5000
    ) %>%
    mutate(
      candidate = cand$candidate,
      candidate_start = cand$start,
      candidate_end = cand$end,
      candidate_strand = cand$strand,
      
      distance_to_candidate =
        case_when(
          genomic_end < cand$start ~
            cand$start - genomic_end,
          
          genomic_start > cand$end ~
            genomic_start - cand$end,
          
          TRUE ~ 0
        ),
      
      overlaps_candidate =
        genomic_start <= cand$end &
        genomic_end >= cand$start,
      
      relationship =
        case_when(
          overlaps_candidate ~
            "OVERLAPS candidate",
          
          genomic_end < cand$start ~
            "Genomically left of candidate",
          
          genomic_start > cand$end ~
            "Genomically right of candidate",
          
          TRUE ~
            "Other"
        )
    ) %>%
    arrange(distance_to_candidate)
  
}) %>%
  bind_rows()

keep_cols <- intersect(
  c(
    "candidate",
    "candidate_start",
    "candidate_end",
    "candidate_strand",
    "distance_to_candidate",
    "relationship",
    "gene_id",
    "current_locus_tag",
    "gene_symbol",
    "protein_id",
    "product",
    "feature_type",
    "genomic_start",
    "genomic_end",
    "strand",
    "is_rna_feature",
    "is_srna",
    "srna_standard_name"
  ),
  names(neighborhoods)
)

result <- neighborhoods %>%
  select(all_of(keep_cols))

cat("\n========================================\n")
cat("TRUE GENOMIC NEIGHBORHOODS\n")
cat("========================================\n\n")

print(
  result,
  n = Inf,
  width = Inf
)

# Closest features
closest <- result %>%
  group_by(candidate) %>%
  slice_min(
    distance_to_candidate,
    n = 10,
    with_ties = TRUE
  ) %>%
  ungroup()

cat("\n========================================\n")
cat("CLOSEST FEATURES\n")
cat("========================================\n\n")

print(
  closest,
  n = Inf,
  width = Inf
)

outdir <- file.path(
  project,
  "17C_true_genomic_neighborhood"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  result,
  file.path(
    outdir,
    "Candidate_5kb_Genomic_Neighborhood.csv"
  )
)

write_csv(
  closest,
  file.path(
    outdir,
    "Candidate_Closest_Features.csv"
  )
)

cat(
  "\nSaved to:",
  outdir,
  "\n"
)