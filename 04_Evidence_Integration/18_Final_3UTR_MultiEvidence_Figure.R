# ============================================================
# SCRIPT 18 v2
# Publication-quality Figure 4
#
# Genomic and multi-evidence characterization of putative
# 3'UTR-associated RNA candidates
#
# Panels:
# A. Genomic context of RS04120
# B. Genomic context of RS05585
# C. Term-seq profile around RS04120
# D. Term-seq profile around RS05585
# E. Independent evidence matrix
#
# IMPORTANT:
# Evidence layers remain separate.
# No integrated biological-priority score.
# No claim of independent processing.
# No claim of direct regulation.
# No causal claim.
# ============================================================


# ============================================================
# 01. PACKAGES
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "patchwork",
  "grid"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(grid)


# ============================================================
# 02. PROJECT PATHS
# ============================================================

project <- "D:/Bac-sRNA"

outdir <- file.path(
  project,
  "18_final_3UTR_multievidence_figure_v2"
)

figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tabdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 03. INPUT FILES
# ============================================================

evidence_file <- file.path(
  project,
  "17D_final_3UTR_evidence",
  "FINAL_3UTR_Candidate_Evidence_Integration.csv"
)

neighborhood_file <- file.path(
  project,
  "17C_true_genomic_neighborhood",
  "Candidate_5kb_Genomic_Neighborhood.csv"
)

term_file <- file.path(
  project,
  "17A_TermSeq_3UTR_boundary_audit",
  "tables",
  "TermSeq_Local_Windows.csv"
)


# ============================================================
# 04. INPUT CHECK
# ============================================================

input_check <- tibble(
  input = c(
    "Integrated evidence",
    "Genomic neighborhood",
    "Term-seq local windows"
  ),
  path = c(
    evidence_file,
    neighborhood_file,
    term_file
  ),
  exists = c(
    file.exists(evidence_file),
    file.exists(neighborhood_file),
    file.exists(term_file)
  )
)

cat("\n================ INPUT CHECK ================\n")

print(
  input_check,
  n = Inf,
  width = Inf
)

if (any(!input_check$exists)) {
  stop("One or more required input files are missing.")
}


# ============================================================
# 05. READ DATA
# ============================================================

evidence <- read_csv(
  evidence_file,
  show_col_types = FALSE
)

neighborhood <- read_csv(
  neighborhood_file,
  show_col_types = FALSE
)

term <- read_csv(
  term_file,
  show_col_types = FALSE
)


# ============================================================
# 06. STANDARDIZE CANDIDATE NAMES
# ============================================================

standard_candidate <- function(x) {
  
  case_when(
    
    str_detect(
      as.character(x),
      regex("RS04120", ignore_case = TRUE)
    ) ~ "3UTR-RS04120",
    
    str_detect(
      as.character(x),
      regex("RS05585", ignore_case = TRUE)
    ) ~ "3UTR-RS05585",
    
    TRUE ~ NA_character_
  )
}


evidence <- evidence %>%
  mutate(
    candidate_std =
      standard_candidate(candidate_short)
  )


neighborhood <- neighborhood %>%
  mutate(
    candidate_std =
      standard_candidate(candidate)
  )


term <- term %>%
  mutate(
    candidate_std =
      standard_candidate(candidate)
  )


# ============================================================
# 07. CANDIDATE EVIDENCE TABLE
# ============================================================

candidate_data <- evidence %>%
  
  filter(
    candidate_std %in%
      c(
        "3UTR-RS04120",
        "3UTR-RS05585"
      )
  ) %>%
  
  transmute(
    
    candidate_std,
    
    candidate_start =
      as.numeric(candidate_start),
    
    candidate_end =
      as.numeric(candidate_end),
    
    candidate_strand,
    
    putative_parent_gene,
    
    putative_parent_locus,
    
    distance_from_parent_3prime_nt =
      as.numeric(
        distance_from_parent_3prime_nt
      ),
    
    target =
      target_label,
    
    supporting_rows =
      as.numeric(
        number_of_supporting_rows
      ),
    
    hybrid_count =
      as.numeric(
        total_hybrid_count
      ),
    
    experiments =
      as.numeric(
        maximum_number_of_experiments
      ),
    
    connection_score =
      as.numeric(
        best_connection_score
      ),
    
    clash_score =
      as.numeric(
        clash_evidence_score
      ),
    
    clash_category =
      clash_support_category,
    
    target_log2FC =
      as.numeric(
        target_log2_fold_change
      ),
    
    target_padj =
      as.numeric(
        target_deseq_adjusted_p
      )
  )


# ============================================================
# 08. PREPARE GENOMIC DATA
# ============================================================

genomic <- neighborhood %>%
  
  filter(
    candidate_std %in%
      c(
        "3UTR-RS04120",
        "3UTR-RS05585"
      )
  ) %>%
  
  mutate(
    
    genomic_start =
      as.numeric(genomic_start),
    
    genomic_end =
      as.numeric(genomic_end),
    
    midpoint =
      (
        genomic_start +
          genomic_end
      ) / 2,
    
    gene_display =
      case_when(
        
        !is.na(gene_symbol) &
          gene_symbol != "" ~
          
          gene_symbol,
        
        !is.na(current_locus_tag) &
          current_locus_tag != "" ~
          
          current_locus_tag,
        
        TRUE ~
          
          gene_id
      )
  )


# ============================================================
# 09. GENOMIC PANEL FUNCTION
# ============================================================

make_genomic_panel <- function(
    candidate_name,
    panel_letter
) {
  
  cand <- candidate_data %>%
    filter(
      candidate_std == candidate_name
    )
  
  genes <- genomic %>%
    filter(
      candidate_std == candidate_name
    )
  
  if (nrow(cand) != 1) {
    stop(
      paste(
        "Unexpected candidate row count for",
        candidate_name
      )
    )
  }
  
  
  center <- mean(
    c(
      cand$candidate_start,
      cand$candidate_end
    )
  )
  
  
  # Keep a tighter genomic region for readability
  window_size <- 4000
  
  xmin <- center - window_size
  xmax <- center + window_size
  
  
  genes <- genes %>%
    
    filter(
      genomic_end >= xmin,
      genomic_start <= xmax
    ) %>%
    
    mutate(
      
      y =
        if_else(
          strand == "+",
          1.00,
          0.45
        ),
      
      label_y =
        if_else(
          strand == "+",
          1.17,
          0.27
        ),
      
      is_parent =
        gene_display ==
        cand$putative_parent_gene,
      
      is_nrdF =
        gene_display == "nrdF",
      
      emphasis =
        is_parent | is_nrdF
    )
  
  
  # Candidate display position
  candidate_y <- 1.48
  
  
  # Distance annotation
  distance_text <- paste0(
    cand$distance_from_parent_3prime_nt,
    " nt from ",
    cand$putative_parent_gene,
    " 3′ end"
  )
  
  
  # Context-specific subtitle
  subtitle_text <- if (
    candidate_name == "3UTR-RS04120"
  ) {
    
    "Putative nrdF 3′UTR-associated RNA candidate"
    
  } else {
    
    "Putative folD 3′UTR-associated RNA candidate; CLASH target = nrdF"
  }
  
  
  # --------------------------------------------
  # Plot
  # --------------------------------------------
  
  p <- ggplot() +
    
    # genomic baseline
    geom_segment(
      aes(
        x = xmin,
        xend = xmax,
        y = 0.72,
        yend = 0.72
      ),
      linewidth = 0.25
    ) +
    
    # genes
    geom_segment(
      data = genes,
      aes(
        x = genomic_start,
        xend = genomic_end,
        y = y,
        yend = y,
        linewidth = emphasis
      ),
      lineend = "butt"
    ) +
    
    scale_linewidth_manual(
      values = c(
        `FALSE` = 3.2,
        `TRUE` = 5.2
      ),
      guide = "none"
    ) +
    
    # plus-strand arrows
    geom_segment(
      data = genes %>%
        filter(strand == "+"),
      aes(
        x =
          genomic_end -
          pmin(
            100,
            pmax(
              20,
              (
                genomic_end -
                  genomic_start
              ) * 0.15
            )
          ),
        
        xend =
          genomic_end,
        
        y = y,
        yend = y
      ),
      arrow = arrow(
        length = unit(
          0.12,
          "cm"
        ),
        type = "closed"
      ),
      linewidth = 0.9
    ) +
    
    # minus-strand arrows
    geom_segment(
      data = genes %>%
        filter(strand == "-"),
      aes(
        x =
          genomic_start +
          pmin(
            100,
            pmax(
              20,
              (
                genomic_end -
                  genomic_start
              ) * 0.15
            )
          ),
        
        xend =
          genomic_start,
        
        y = y,
        yend = y
      ),
      arrow = arrow(
        length = unit(
          0.12,
          "cm"
        ),
        type = "closed"
      ),
      linewidth = 0.9
    ) +
    
    # gene labels
    geom_text(
      data = genes,
      aes(
        x = midpoint,
        y = label_y,
        label = gene_display,
        fontface =
          ifelse(
            emphasis,
            "bold",
            "plain"
          )
      ),
      size = 3.0,
      check_overlap = TRUE
    ) +
    
    # candidate fragment
    geom_segment(
      aes(
        x = cand$candidate_start,
        xend = cand$candidate_end,
        y = candidate_y,
        yend = candidate_y
      ),
      linewidth = 4.2,
      lineend = "butt"
    ) +
    
    # candidate-to-genome guide
    geom_segment(
      aes(
        x = center,
        xend = center,
        y = 1.02,
        yend = 1.41
      ),
      linewidth = 0.5
    ) +
    
    annotate(
      "text",
      x = center,
      y = 1.61,
      label = candidate_name,
      size = 3.3,
      fontface = "bold"
    ) +
    
    annotate(
      "text",
      x = center,
      y = 1.82,
      label = distance_text,
      size = 2.8
    ) +
    
    coord_cartesian(
      xlim = c(
        xmin,
        xmax
      ),
      ylim = c(
        0.05,
        2.0
      ),
      clip = "off"
    ) +
    
    labs(
      title =
        paste0(
          panel_letter,
          "  Genomic context of ",
          str_remove(
            candidate_name,
            "^3UTR-"
          )
        ),
      
      subtitle =
        subtitle_text,
      
      x =
        "JKD6008 genomic coordinate (bp)",
      
      y = NULL
    ) +
    
    theme_classic(
      base_size = 10
    ) +
    
    theme(
      
      axis.text.y =
        element_blank(),
      
      axis.ticks.y =
        element_blank(),
      
      axis.line.y =
        element_blank(),
      
      plot.title =
        element_text(
          face = "bold",
          size = 11
        ),
      
      plot.subtitle =
        element_text(
          size = 8.8
        ),
      
      plot.margin =
        margin(
          8,
          12,
          8,
          8
        )
    )
  
  
  return(p)
}


# ============================================================
# 10. PANELS A AND B
# ============================================================

panel_A <- make_genomic_panel(
  "3UTR-RS04120",
  "A"
)

panel_B <- make_genomic_panel(
  "3UTR-RS05585",
  "B"
)


# ============================================================
# 11. TERM-SEQ DATA
# ============================================================

required_term_cols <- c(
  "candidate_std",
  "relative_position",
  "signal",
  "replicate"
)

missing_term <- setdiff(
  required_term_cols,
  names(term)
)

if (length(missing_term) > 0) {
  
  print(missing_term)
  
  stop(
    "Required Term-seq columns are missing."
  )
}


term_plot <- term %>%
  
  filter(
    candidate_std %in%
      c(
        "3UTR-RS04120",
        "3UTR-RS05585"
      )
  ) %>%
  
  mutate(
    
    relative_position =
      as.numeric(relative_position),
    
    signal =
      as.numeric(signal),
    
    replicate =
      as.factor(replicate)
  ) %>%
  
  filter(
    !is.na(relative_position),
    !is.na(signal),
    relative_position >= -100,
    relative_position <= 100
  )


# ============================================================
# 12. TERM-SEQ PANEL FUNCTION
# ============================================================

make_term_panel <- function(
    candidate_name,
    panel_letter
) {
  
  df <- term_plot %>%
    
    filter(
      candidate_std == candidate_name
    )
  
  
  candidate_label <- str_remove(
    candidate_name,
    "^3UTR-"
  )
  
  
  ggplot(
    df,
    aes(
      x = relative_position,
      y = signal,
      group = replicate,
      linetype = replicate
    )
  ) +
    
    # ±5 nt boundary window
    annotate(
      "rect",
      xmin = -5,
      xmax = 5,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.08
    ) +
    
    geom_line(
      linewidth = 0.72
    ) +
    
    # exact expected boundary
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.75
    ) +
    
    annotate(
      "text",
      x = 4,
      y = Inf,
      label = "Expected boundary",
      hjust = 0,
      vjust = 1.4,
      size = 2.8
    ) +
    
    scale_x_continuous(
      breaks = c(
        -100,
        -50,
        0,
        50,
        100
      )
    ) +
    
    labs(
      title =
        paste0(
          panel_letter,
          "  Term-seq profile around ",
          candidate_label
        ),
      
      subtitle =
        "Signal within ±100 nt; shaded region indicates ±5 nt around expected boundary",
      
      x =
        "Position relative to candidate boundary (nt)",
      
      y =
        "Term-seq 3′-end signal",
      
      linetype =
        "Replicate"
    ) +
    
    theme_classic(
      base_size = 10
    ) +
    
    theme(
      
      plot.title =
        element_text(
          face = "bold",
          size = 11
        ),
      
      plot.subtitle =
        element_text(
          size = 8.5
        ),
      
      legend.position =
        "top",
      
      legend.title =
        element_text(
          size = 8
        ),
      
      legend.text =
        element_text(
          size = 8
        )
    )
}


# ============================================================
# 13. PANELS C AND D
# ============================================================

panel_C <- make_term_panel(
  "3UTR-RS04120",
  "C"
)

panel_D <- make_term_panel(
  "3UTR-RS05585",
  "D"
)


# ============================================================
# 14. PANEL E DATA
# ============================================================

evidence_matrix <- candidate_data %>%
  
  mutate(
    
    Candidate =
      case_when(
        
        candidate_std ==
          "3UTR-RS04120" ~
          "RS04120",
        
        candidate_std ==
          "3UTR-RS05585" ~
          "RS05585"
      ),
    
    Parent =
      putative_parent_gene,
    
    Distance =
      paste0(
        distance_from_parent_3prime_nt,
        " nt"
      ),
    
    Hybrids =
      as.character(
        hybrid_count
      ),
    
    Experiments =
      as.character(
        experiments
      ),
    
    Connection =
      sprintf(
        "%.3f",
        connection_score
      ),
    
    nrdF_log2FC =
      sprintf(
        "%.3f",
        target_log2FC
      ),
    
    nrdF_padj =
      format(
        target_padj,
        scientific = TRUE,
        digits = 3
      ),
    
    TermSeq =
      "No"
  ) %>%
  
  select(
    Candidate,
    Parent,
    Distance,
    Hybrids,
    Experiments,
    Connection,
    nrdF_log2FC,
    nrdF_padj,
    TermSeq
  )


# ============================================================
# 15. PANEL E — MANUAL CLEAN TABLE
# ============================================================

row_labels <- c(
  "Putative parent",
  "Distance from parent 3′ end",
  "CLASH hybrids",
  "CLASH experiments",
  "Connection score",
  "nrdF log2FC",
  "nrdF adjusted P",
  "Sharp reproducible Term-seq boundary"
)


rs04120 <- evidence_matrix %>%
  filter(Candidate == "RS04120")

rs05585 <- evidence_matrix %>%
  filter(Candidate == "RS05585")


rs04120_values <- c(
  rs04120$Parent,
  rs04120$Distance,
  rs04120$Hybrids,
  rs04120$Experiments,
  rs04120$Connection,
  rs04120$nrdF_log2FC,
  rs04120$nrdF_padj,
  rs04120$TermSeq
)


rs05585_values <- c(
  rs05585$Parent,
  rs05585$Distance,
  rs05585$Hybrids,
  rs05585$Experiments,
  rs05585$Connection,
  rs05585$nrdF_log2FC,
  rs05585$nrdF_padj,
  rs05585$TermSeq
)


table_df <- tibble(
  
  Evidence =
    factor(
      row_labels,
      levels = rev(row_labels)
    ),
  
  RS04120 =
    rs04120_values,
  
  RS05585 =
    rs05585_values
)


table_long <- table_df %>%
  
  pivot_longer(
    cols = c(
      RS04120,
      RS05585
    ),
    names_to = "Candidate",
    values_to = "Value"
  ) %>%
  
  mutate(
    
    Candidate =
      factor(
        Candidate,
        levels = c(
          "RS04120",
          "RS05585"
        )
      )
  )


# ============================================================
# 16. PANEL E PLOT
# ============================================================

panel_E <- ggplot(
  table_long,
  aes(
    x = Candidate,
    y = Evidence
  )
) +
  
  geom_tile(
    fill = "white",
    linewidth = 0.45
  ) +
  
  geom_text(
    aes(
      label = Value
    ),
    size = 3.4
  ) +
  
  labs(
    title =
      "E  Independent evidence layers",
    
    subtitle =
      "Evidence is reported separately; no integrated biological-priority score was calculated",
    
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 10
  ) +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_text(
        face = "bold",
        size = 9.5
      ),
    
    axis.text.y =
      element_text(
        size = 8.8
      ),
    
    axis.ticks =
      element_blank(),
    
    plot.title =
      element_text(
        face = "bold",
        size = 11
      ),
    
    plot.subtitle =
      element_text(
        size = 8.7
      ),
    
    plot.margin =
      margin(
        8,
        8,
        4,
        8
      )
  )


# ============================================================
# 17. INTERPRETATION STRIP
# ============================================================

interpretation_panel <- ggplot() +
  
  annotate(
    "text",
    x = 0,
    y = 1,
    hjust = 0,
    vjust = 0.5,
    size = 3.6,
    fontface = "bold",
    label =
      paste(
        "Interpretation:",
        "genomic context and CLASH support candidate 3′UTR association;",
        "Term-seq does not establish independent processing."
      )
  ) +
  
  xlim(
    0,
    1
  ) +
  
  ylim(
    0,
    2
  ) +
  
  theme_void()


# ============================================================
# 18. COMBINE FIGURE
# ============================================================

top_row <- panel_A + panel_B +
  plot_layout(
    ncol = 2
  )


middle_row <- panel_C + panel_D +
  plot_layout(
    ncol = 2
  )


figure4 <- (
  
  top_row /
    
    middle_row /
    
    panel_E /
    
    interpretation_panel
  
) +
  
  plot_layout(
    
    heights = c(
      0.95,
      1.20,
      1.10,
      0.18
    )
  ) +
  
  plot_annotation(
    
    title =
      paste(
        "Genomic and multi-evidence characterization",
        "of putative 3′UTR-associated RNA candidates"
      ),
    
    subtitle =
      paste(
        "Genomic context, RNase III-CLASH interaction evidence,",
        "target transcriptomic response and Term-seq boundary evidence",
        "are evaluated as independent evidence layers."
      ),
    
    theme =
      theme(
        
        plot.title =
          element_text(
            face = "bold",
            size = 15
          ),
        
        plot.subtitle =
          element_text(
            size = 10
          )
      )
  )


# ============================================================
# 19. DISPLAY
# ============================================================

print(figure4)


# ============================================================
# 20. EXPORT PDF
# ============================================================

ggsave(
  
  file.path(
    figdir,
    "Figure4_3UTR_MultiEvidence_FINAL_v2.pdf"
  ),
  
  plot =
    figure4,
  
  width =
    12,
  
  height =
    12.7,
  
  units =
    "in",
  
  device =
    cairo_pdf,
  
  limitsize =
    FALSE
)


# ============================================================
# 21. EXPORT PNG — 600 DPI
# ============================================================

ggsave(
  
  file.path(
    figdir,
    "Figure4_3UTR_MultiEvidence_FINAL_v2_600dpi.png"
  ),
  
  plot =
    figure4,
  
  width =
    12,
  
  height =
    12.7,
  
  units =
    "in",
  
  dpi =
    600,
  
  limitsize =
    FALSE
)


# ============================================================
# 22. EXPORT TIFF — 600 DPI
# ============================================================

ggsave(
  
  file.path(
    figdir,
    "Figure4_3UTR_MultiEvidence_FINAL_v2_600dpi.tiff"
  ),
  
  plot =
    figure4,
  
  width =
    12,
  
  height =
    12.7,
  
  units =
    "in",
  
  dpi =
    600,
  
  compression =
    "lzw",
  
  limitsize =
    FALSE
)


# ============================================================
# 23. SAVE SOURCE TABLES
# ============================================================

write_csv(
  
  candidate_data,
  
  file.path(
    tabdir,
    "Figure4_Candidate_Evidence.csv"
  )
)


write_csv(
  
  term_plot,
  
  file.path(
    tabdir,
    "Figure4_TermSeq_Profile_Data.csv"
  )
)


write_csv(
  
  table_df,
  
  file.path(
    tabdir,
    "Figure4_Evidence_Matrix.csv"
  )
)


# ============================================================
# 24. FINAL FIGURE LEGEND
# ============================================================

figure_legend <- paste(
  
  "Figure 4. Genomic and multi-evidence characterization of two",
  "putative 3′UTR-associated RNA candidates interacting with nrdF.",
  
  "(A) The RS04120-labelled CLASH RNA fragment is located 53 nt",
  "from the annotated 3′ end of nrdF on the same strand, supporting",
  "its classification as a putative nrdF 3′UTR-associated RNA candidate.",
  
  "(B) The RS05585-labelled fragment is located 177 nt from the",
  "annotated 3′ end of folD on the same strand, supporting a putative",
  "folD 3′UTR-associated origin; its CLASH-supported interaction target",
  "is nrdF.",
  
  "(C,D) Term-seq 3′-end signal within ±100 nt of each candidate",
  "boundary. Dashed lines denote the expected candidate boundary and",
  "the shaded regions indicate ±5 nt around that position.",
  
  "(E) Independent evidence summary. RS04120 showed greater CLASH",
  "hybrid abundance and experimental recurrence than RS05585, whereas",
  "both interactions converged on nrdF, which was significantly",
  "decreased following vancomycin exposure.",
  
  "No reproducible sharp Term-seq boundary was established for either",
  "candidate. The combined evidence therefore supports conservative",
  "classification as putative 3′UTR-associated RNA candidates but does",
  "not establish independent RNA processing, direct regulation of nrdF,",
  "or a causal role in the vancomycin response."
)


writeLines(
  
  figure_legend,
  
  file.path(
    outdir,
    "Figure4_FINAL_Legend.txt"
  )
)


# ============================================================
# 25. MANUSCRIPT RESULTS PARAGRAPH
# ============================================================

results_text <- paste(
  
  "Genomic-context analysis refined the interpretation of the two",
  "CLASH-defined RNA fragments associated with nrdF. The RS04120-labelled",
  "fragment was located 53 nt from the annotated 3′ end of nrdF on the",
  "same strand, supporting its classification as a putative nrdF",
  "3′UTR-associated RNA candidate. In contrast, the RS05585-labelled",
  "fragment was located 177 nt from the annotated 3′ end of folD on the",
  "same strand, consistent with a putative folD 3′UTR-associated RNA",
  "candidate. RNase III-CLASH recovered interactions of both fragments",
  "with nrdF, with greater support for RS04120 (15 hybrids across two",
  "experiments; connection score 0.289) than for RS05585 (7 hybrids in",
  "one experiment; connection score 0.106). Independently, nrdF was",
  "significantly decreased following vancomycin exposure",
  "(log2FC = -1.118; adjusted P = 3.04 × 10^-11). Term-seq analysis",
  "did not establish a reproducible sharp 3′ boundary for either",
  "candidate. Accordingly, the combined evidence supports their",
  "conservative designation as putative 3′UTR-associated RNA candidates,",
  "while independent processing and direct regulation of nrdF remain",
  "unproven."
)


writeLines(
  
  results_text,
  
  file.path(
    outdir,
    "Figure4_FINAL_Results_Paragraph.txt"
  )
)


# ============================================================
# 26. SESSION INFO
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
# 27. FINAL CONSOLE AUDIT
# ============================================================

cat("\n\n")
cat("============================================\n")
cat("SCRIPT 18 v2 COMPLETED\n")
cat("============================================\n")


cat("\nFINAL EVIDENCE MATRIX\n")
cat("--------------------------------------------\n")

print(
  table_df,
  n = Inf,
  width = Inf
)


cat("\n\nFINAL INTERPRETATION\n")
cat("--------------------------------------------\n")

cat(
  paste0(
    
    "RS04120:\n",
    "  Genomic context: nrdF-associated\n",
    "  Distance: 53 nt from nrdF 3′ end\n",
    "  CLASH: 15 hybrids / 2 experiments\n",
    "  Connection score: 0.289\n\n",
    
    "RS05585:\n",
    "  Genomic context: folD-associated\n",
    "  Distance: 177 nt from folD 3′ end\n",
    "  CLASH target: nrdF\n",
    "  CLASH: 7 hybrids / 1 experiment\n",
    "  Connection score: 0.106\n\n",
    
    "Shared target response:\n",
    "  nrdF log2FC = -1.118\n",
    "  adjusted P = 3.04e-11\n\n",
    
    "Term-seq:\n",
    "  No reproducible sharp candidate boundary established.\n\n",
    
    "Manuscript terminology:\n",
    "  putative 3′UTR-associated RNA candidate\n\n",
    
    "Do NOT claim:\n",
    "  independent processing\n",
    "  direct regulation\n",
    "  causal vancomycin-response mechanism\n"
  )
)


cat("\n\nOUTPUT DIRECTORY:\n")
cat(outdir, "\n")


cat("\nMain publication files:\n")

cat(
  "Figure4_3UTR_MultiEvidence_FINAL_v2.pdf\n"
)

cat(
  "Figure4_3UTR_MultiEvidence_FINAL_v2_600dpi.png\n"
)

cat(
  "Figure4_3UTR_MultiEvidence_FINAL_v2_600dpi.tiff\n"
)

cat(
  "Figure4_FINAL_Legend.txt\n"
)

cat("\n============================================\n")
cat("END SCRIPT 18 v2\n")
cat("============================================\n")