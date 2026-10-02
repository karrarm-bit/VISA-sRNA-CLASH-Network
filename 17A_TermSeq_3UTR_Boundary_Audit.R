# ============================================================
# SCRIPT 17A
# Term-seq 3' Boundary Audit for Candidate 3'UTR-associated RNAs
# Staphylococcus aureus JKD6008
#
# Candidates:
#   3UTR-RS04120 : 811819-811847 (+)
#   3UTR-RS05585 : 1084557-1084613 (-)
#
# IMPORTANT:
# Term-seq supports 3' end / boundary evidence.
# It does NOT by itself prove independent RNA processing,
# release from the parent transcript, or regulatory function.
# ============================================================


# ============================================================
# 01. PATHS
# ============================================================

project <- "D:/Bac-sRNA"

input_dir <- file.path(
  project,
  "07_processed_data",
  "downloaded_files",
  "GSE158830"
)

output_dir <- file.path(
  project,
  "17A_TermSeq_3UTR_boundary_audit"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")

dir.create(table_dir, showWarnings = FALSE)
dir.create(figure_dir, showWarnings = FALSE)


# ============================================================
# 02. PACKAGES
# ============================================================

packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "tibble",
  "purrr"
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
library(ggplot2)
library(tibble)
library(purrr)


# ============================================================
# 03. CANDIDATES
# ============================================================

candidates <- tibble::tribble(
  
  ~candidate,
  ~start,
  ~end,
  ~strand,
  
  "3UTR-RS04120",
  811819,
  811847,
  "+",
  
  "3UTR-RS05585",
  1084557,
  1084613,
  "-"
)


# ============================================================
# 04. EXPECTED 3' BOUNDARY
#
# Plus strand:
# 3' end = genomic end
#
# Minus strand:
# transcription runs opposite direction,
# therefore 3' end = genomic start
# ============================================================

candidates <- candidates |>
  
  mutate(
    
    expected_3prime =
      ifelse(
        strand == "+",
        end,
        start
      )
  )


print(candidates)


# ============================================================
# 05. LOCATE TERM-SEQ FILES
# ============================================================

term_files <- list.files(
  
  input_dir,
  
  pattern =
    "termseq.*(plus|minus).*\\.gr\\.gz$",
  
  full.names = TRUE,
  
  ignore.case = TRUE
)


cat(
  "\nTERM-SEQ FILES FOUND:",
  length(term_files),
  "\n\n"
)

print(term_files)


if (length(term_files) != 4) {
  
  warning(
    paste(
      "Expected 4 Term-seq files",
      "(2 replicates x 2 strands), found",
      length(term_files)
    )
  )
}


# ============================================================
# 06. FILE METADATA
# ============================================================

term_meta <- tibble(
  
  file = term_files,
  
  filename =
    basename(term_files)
) |>
  
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


print(term_meta)


# ============================================================
# 07. ROBUST .gr READER
#
# Expected common formats:
# chromosome position signal
# OR
# position signal
#
# We inspect number of columns automatically.
# ============================================================

read_term_file <- function(
    file,
    replicate,
    strand
) {
  
  x <- readr::read_table(
    file,
    col_names = FALSE,
    comment = "#",
    show_col_types = FALSE,
    progress = FALSE
  )
  
  
  cat(
    "\nReading:",
    basename(file),
    "\n"
  )
  
  cat(
    "Rows:",
    nrow(x),
    "Columns:",
    ncol(x),
    "\n"
  )
  
  
  if (ncol(x) == 2) {
    
    names(x) <- c(
      "position",
      "signal"
    )
    
    x <- x |>
      
      mutate(
        chromosome = NA_character_
      ) |>
      
      select(
        chromosome,
        position,
        signal
      )
    
  } else if (ncol(x) >= 3) {
    
    x <- x[, 1:3]
    
    names(x) <- c(
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
  
  
  x |>
    
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
    ) |>
    
    filter(
      !is.na(position),
      !is.na(signal)
    )
}


# ============================================================
# 08. READ ALL TERM-SEQ DATA
# ============================================================

term_data <- purrr::pmap_dfr(
  
  term_meta,
  
  function(
    file,
    filename,
    replicate,
    strand
  ) {
    
    read_term_file(
      file,
      replicate,
      strand
    )
  }
)


cat(
  "\nTotal Term-seq signal rows:",
  nrow(term_data),
  "\n"
)


# ============================================================
# 09. BASIC DATA AUDIT
# ============================================================

data_audit <- term_data |>
  
  group_by(
    replicate,
    strand
  ) |>
  
  summarise(
    
    rows = n(),
    
    min_position =
      min(position),
    
    max_position =
      max(position),
    
    total_signal =
      sum(
        signal,
        na.rm = TRUE
      ),
    
    maximum_signal =
      max(
        signal,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


write_csv(
  data_audit,
  file.path(
    table_dir,
    "TermSeq_Data_Audit.csv"
  )
)


# ============================================================
# 10. EXTRACT ±100 nt AROUND EXPECTED 3' END
# ============================================================

window_size <- 100


extract_candidate_window <- function(
    candidate_name,
    expected_boundary,
    candidate_strand
) {
  
  term_data |>
    
    filter(
      strand == candidate_strand,
      position >= expected_boundary - window_size,
      position <= expected_boundary + window_size
    ) |>
    
    mutate(
      
      candidate =
        candidate_name,
      
      expected_3prime =
        expected_boundary,
      
      relative_position =
        position -
        expected_boundary
    )
}


candidate_windows <- purrr::pmap_dfr(
  
  candidates,
  
  function(
    candidate,
    start,
    end,
    strand,
    expected_3prime
  ) {
    
    extract_candidate_window(
      candidate,
      expected_3prime,
      strand
    )
  }
)


# ============================================================
# 11. DEFINE LOCAL PEAK
#
# Boundary support window:
# +/- 5 nt from expected 3' end
#
# Local background:
# 20-100 nt away on either side.
# ============================================================

boundary_radius <- 5

background_inner <- 20
background_outer <- 100


candidate_windows <- candidate_windows |>
  
  mutate(
    
    region =
      case_when(
        
        abs(relative_position) <=
          boundary_radius ~
          "Expected boundary ±5 nt",
        
        abs(relative_position) >=
          background_inner &
          abs(relative_position) <=
          background_outer ~
          "Local background",
        
        TRUE ~
          "Intermediate region"
      )
  )


# ============================================================
# 12. PER-REPLICATE BOUNDARY METRICS
# ============================================================

boundary_metrics <- candidate_windows |>
  
  group_by(
    candidate,
    replicate,
    strand,
    expected_3prime
  ) |>
  
  summarise(
    
    boundary_max_signal =
      {
        z <- signal[
          abs(relative_position) <=
            boundary_radius
        ]
        
        if (length(z) == 0) {
          0
        } else {
          max(
            z,
            na.rm = TRUE
          )
        }
      },
    
    boundary_total_signal =
      sum(
        signal[
          abs(relative_position) <=
            boundary_radius
        ],
        na.rm = TRUE
      ),
    
    background_median_signal =
      {
        z <- signal[
          abs(relative_position) >=
            background_inner &
            abs(relative_position) <=
            background_outer
        ]
        
        if (length(z) == 0) {
          0
        } else {
          median(
            z,
            na.rm = TRUE
          )
        }
      },
    
    background_mean_signal =
      {
        z <- signal[
          abs(relative_position) >=
            background_inner &
            abs(relative_position) <=
            background_outer
        ]
        
        if (length(z) == 0) {
          0
        } else {
          mean(
            z,
            na.rm = TRUE
          )
        }
      },
    
    strongest_position =
      {
        idx <- which.max(signal)
        
        if (length(idx) == 0) {
          NA_real_
        } else {
          position[idx]
        }
      },
    
    strongest_relative_position =
      {
        idx <- which.max(signal)
        
        if (length(idx) == 0) {
          NA_real_
        } else {
          relative_position[idx]
        }
      },
    
    strongest_signal =
      ifelse(
        length(signal) == 0,
        0,
        max(
          signal,
          na.rm = TRUE
        )
      ),
    
    .groups = "drop"
  )


# ============================================================
# 13. ENRICHMENT CALCULATION
#
# Pseudocount avoids division by zero.
# ============================================================

boundary_metrics <- boundary_metrics |>
  
  mutate(
    
    boundary_to_background_ratio =
      (
        boundary_max_signal + 1
      ) /
      (
        background_median_signal + 1
      ),
    
    peak_near_expected_boundary =
      abs(
        strongest_relative_position
      ) <= boundary_radius
  )


# ============================================================
# 14. REPLICATE-LEVEL SUPPORT
#
# Conservative descriptive thresholds.
#
# IMPORTANT:
# These are audit categories, not biological proof.
# ============================================================

boundary_metrics <- boundary_metrics |>
  
  mutate(
    
    replicate_boundary_support =
      case_when(
        
        boundary_max_signal > 0 &
          boundary_to_background_ratio >= 5 &
          peak_near_expected_boundary ~
          "Strong local 3′-end signal",
        
        boundary_max_signal > 0 &
          boundary_to_background_ratio >= 2 ~
          "Local 3′-end enrichment",
        
        boundary_max_signal > 0 ~
          "3′-end signal detected",
        
        TRUE ~
          "No 3′-end signal detected"
      )
  )


# ============================================================
# 15. CROSS-REPLICATE SUMMARY
# ============================================================

candidate_summary <- boundary_metrics |>
  
  group_by(
    candidate,
    strand,
    expected_3prime
  ) |>
  
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
        boundary_to_background_ratio >= 2 &
          boundary_max_signal > 0,
        na.rm = TRUE
      ),
    
    replicates_with_strong_local_peak =
      sum(
        boundary_to_background_ratio >= 5 &
          peak_near_expected_boundary &
          boundary_max_signal > 0,
        na.rm = TRUE
      ),
    
    median_boundary_signal =
      median(
        boundary_max_signal,
        na.rm = TRUE
      ),
    
    median_enrichment_ratio =
      median(
        boundary_to_background_ratio,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ============================================================
# 16. FINAL CONSERVATIVE INTERPRETATION
# ============================================================

candidate_summary <- candidate_summary |>
  
  mutate(
    
    TermSeq_boundary_class =
      case_when(
        
        replicates_with_strong_local_peak >= 2 ~
          "Reproducible strong 3′-boundary support",
        
        replicates_with_enrichment >= 2 ~
          "Reproducible local 3′-boundary enrichment",
        
        replicates_with_signal >= 2 ~
          "Reproducible 3′-end signal without strong enrichment",
        
        replicates_with_signal == 1 ~
          "Single-replicate 3′-end signal",
        
        TRUE ~
          "No Term-seq boundary signal detected"
      ),
    
    manuscript_interpretation =
      case_when(
        
        replicates_with_strong_local_peak >= 2 ~
          paste(
            "Term-seq shows a reproducible local 3′-end peak",
            "near the CLASH-defined candidate boundary;",
            "this supports a discrete 3′ boundary but does not",
            "by itself establish independent RNA processing."
          ),
        
        replicates_with_enrichment >= 2 ~
          paste(
            "Term-seq shows reproducible enrichment near the",
            "candidate 3′ boundary; this provides boundary",
            "support but does not establish an independently",
            "processed regulatory RNA."
          ),
        
        replicates_with_signal >= 1 ~
          paste(
            "Term-seq signal is detectable near the candidate",
            "boundary, but evidence is insufficient to claim",
            "a reproducible discrete 3′ end."
          ),
        
        TRUE ~
          paste(
            "No local Term-seq signal was detected near the",
            "candidate boundary under the analyzed conditions."
          )
      )
  )


# ============================================================
# 17. SAVE TABLES
# ============================================================

write_csv(
  candidates,
  file.path(
    table_dir,
    "Candidate_Coordinates.csv"
  )
)


write_csv(
  candidate_windows,
  file.path(
    table_dir,
    "TermSeq_Local_Windows.csv"
  )
)


write_csv(
  boundary_metrics,
  file.path(
    table_dir,
    "TermSeq_Per_Replicate_Boundary_Metrics.csv"
  )
)


write_csv(
  candidate_summary,
  file.path(
    table_dir,
    "TermSeq_Candidate_Boundary_Summary.csv"
  )
)


# ============================================================
# 18. PLOT DATA
# ============================================================

plot_data <- candidate_windows |>
  
  filter(
    abs(relative_position) <= 75
  )


# ============================================================
# 19. FIGURE
# ============================================================

p <- ggplot(
  
  plot_data,
  
  aes(
    x = relative_position,
    y = signal
  )
) +
  
  geom_line(
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  facet_grid(
    candidate ~ replicate,
    scales = "free_y"
  ) +
  
  labs(
    
    title =
      "Term-seq signal around candidate 3′ boundaries",
    
    subtitle =
      paste(
        "Position 0 represents the CLASH-defined",
        "candidate 3′ boundary"
      ),
    
    x =
      "Position relative to expected 3′ boundary (nt)",
    
    y =
      "Term-seq 3′-end signal"
  ) +
  
  theme_classic(
    base_size = 10
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    strip.text =
      element_text(
        face = "bold"
      )
  )


print(p)


ggsave(
  file.path(
    figure_dir,
    "Figure_TermSeq_Candidate_3prime_Boundaries.pdf"
  ),
  p,
  width = 9,
  height = 6,
  device = cairo_pdf,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_TermSeq_Candidate_3prime_Boundaries_600dpi.png"
  ),
  p,
  width = 9,
  height = 6,
  dpi = 600,
  bg = "white"
)


# ============================================================
# 20. SANITY CHECK
# ============================================================

if (nrow(candidate_summary) != 2) {
  
  warning(
    paste(
      "Expected 2 candidate RNAs;",
      "observed",
      nrow(candidate_summary)
    )
  )
}


# ============================================================
# 21. METHODOLOGICAL NOTES
# ============================================================

notes <- c(
  
  "SCRIPT 17A - TERM-SEQ BOUNDARY AUDIT",
  
  "",
  
  paste(
    "Candidate coordinates were taken from the",
    "CLASH-derived regulatory RNA annotations."
  ),
  
  paste(
    "For plus-strand candidates, the genomic end",
    "was treated as the expected 3-prime boundary."
  ),
  
  paste(
    "For minus-strand candidates, the genomic start",
    "was treated as the expected 3-prime boundary."
  ),
  
  paste(
    "Term-seq signal was examined in the matching",
    "strand in two available control replicates."
  ),
  
  paste(
    "Signal within +/-5 nt of the expected boundary",
    "was compared descriptively with local signal",
    "20-100 nt away."
  ),
  
  paste(
    "The enrichment categories are descriptive",
    "audit thresholds and are not statistical proof",
    "of RNA processing."
  ),
  
  paste(
    "A reproducible Term-seq 3-prime end supports",
    "a discrete transcript boundary but does not",
    "establish that the candidate is independently",
    "processed or functionally regulatory."
  ),
  
  paste(
    "The available GSE158830 RAW archive used here",
    "contains Term-seq and CLASH-related files but",
    "no locally available dRNA-seq file was identified."
  )
)


writeLines(
  notes,
  file.path(
    output_dir,
    "Script17A_Methodological_Notes.txt"
  )
)


# ============================================================
# 22. SESSION INFO
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
# 23. FINAL REPORT
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "SCRIPT 17A COMPLETED\n"
)

cat(
  "============================================\n\n"
)


cat(
  "TERM-SEQ DATA AUDIT\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  data_audit,
  n = Inf,
  width = Inf
)


cat(
  "\nPER-REPLICATE BOUNDARY METRICS\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  boundary_metrics,
  n = Inf,
  width = Inf
)


cat(
  "\nCANDIDATE BOUNDARY SUMMARY\n"
)

cat(
  "--------------------------------------------\n"
)

print(
  candidate_summary,
  n = Inf,
  width = Inf
)


cat(
  "\n============================================\n"
)

cat(
  "IMPORTANT INTERPRETATION RULE\n"
)

cat(
  "============================================\n"
)

cat(
  paste(
    "\nTerm-seq evidence supports a 3′ boundary only.",
    "\nDo not call a candidate '3′UTR-derived RNA'",
    "solely from this analysis.",
    "\nIndependent processing requires additional",
    "5′-boundary/processing or transcript evidence.\n"
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
  "END SCRIPT 17A\n"
)

cat(
  "============================================\n"
)