# ============================================================
# SCRIPT 22 — FINAL CORRECTED FIGURE 3
# Publication-ready version
#
# Genomic context and independent evidence layers for two
# putative 3′UTR-associated RNA candidates
#
# IMPORTANT:
# - CLASH evidence and transcriptomic evidence remain separate.
# - Term-seq is descriptive boundary evidence only.
# - No integrated biological-priority score is calculated.
# - No claim of independent RNA processing or direct regulation.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "patchwork",
  "grid"
)

missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(grid)


# ============================================================
# 2. PROJECT PATHS
# ============================================================

project <- "D:/Bac-sRNA"

termseq_file <- file.path(
  project,
  "17A_TermSeq_3UTR_boundary_audit",
  "tables",
  "TermSeq_Local_Windows.csv"
)

coord_file <- file.path(
  project,
  "17A_TermSeq_3UTR_boundary_audit",
  "tables",
  "Candidate_Coordinates.csv"
)

neighborhood_file <- file.path(
  project,
  "17C_true_genomic_neighborhood",
  "Candidate_5kb_Genomic_Neighborhood.csv"
)

nrdF_file <- file.path(
  project,
  "14_revised_evidence_layered_network",
  "tables",
  "nrdF_Candidate_Interactions.csv"
)


# NEW OUTPUT FOLDER
outdir <- file.path(
  project,
  "22_final_Figure3_corrected"
)

figdir <- file.path(
  outdir,
  "figures"
)

tabdir <- file.path(
  outdir,
  "tables"
)

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
# 3. CHECK INPUT FILES
# ============================================================

required_files <- c(
  termseq_file,
  coord_file,
  neighborhood_file,
  nrdF_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    paste(
      "Missing required file(s):",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


# ============================================================
# 4. READ DATA
# ============================================================

termseq <- read_csv(
  termseq_file,
  show_col_types = FALSE
)

coords <- read_csv(
  coord_file,
  show_col_types = FALSE
)

neighborhood <- read_csv(
  neighborhood_file,
  show_col_types = FALSE
)

nrdF_net <- read_csv(
  nrdF_file,
  show_col_types = FALSE
)


# ============================================================
# 5. HELPER FUNCTION
# ============================================================

first_existing <- function(df, candidates) {
  
  hit <- candidates[
    candidates %in% names(df)
  ]
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  hit[1]
}


# ============================================================
# 6. STANDARDIZE nrdF INTERACTION DATA
# ============================================================

logfc_col <- first_existing(
  nrdF_net,
  c(
    "target_log2_fold_change",
    "target_log2FoldChange",
    "target_log2FC",
    "target_log2fc"
  )
)

padj_col <- first_existing(
  nrdF_net,
  c(
    "target_deseq_adjusted_p",
    "target_padj",
    "target_padj_final",
    "target_adjusted_p_value"
  )
)

if (is.na(logfc_col)) {
  stop("Could not detect the nrdF target log2FC column.")
}

if (is.na(padj_col)) {
  stop("Could not detect the nrdF target adjusted-P column.")
}


nrdF_std <- nrdF_net %>%
  
  mutate(
    
    candidate = case_when(
      
      str_detect(
        source_label,
        regex(
          "RS04120",
          ignore_case = TRUE
        )
      ) ~ "RS04120",
      
      str_detect(
        source_label,
        regex(
          "RS05585",
          ignore_case = TRUE
        )
      ) ~ "RS05585",
      
      TRUE ~ NA_character_
    ),
    
    target_log2FC =
      as.numeric(.data[[logfc_col]]),
    
    target_padj =
      as.numeric(.data[[padj_col]])
  ) %>%
  
  filter(
    candidate %in% c(
      "RS04120",
      "RS05585"
    )
  )


# ============================================================
# 7. AUDITED CANDIDATE METADATA
# ============================================================

candidate_meta <- tibble(
  
  candidate = c(
    "RS04120",
    "RS05585"
  ),
  
  putative_parent = c(
    "nrdF",
    "folD"
  ),
  
  distance_nt = c(
    53,
    177
  ),
  
  description = c(
    "Putative nrdF 3′UTR-associated RNA candidate",
    "Putative folD 3′UTR-associated RNA candidate"
  )
)


evidence <- candidate_meta %>%
  
  left_join(
    
    nrdF_std %>%
      
      select(
        candidate,
        total_hybrid_count,
        maximum_number_of_experiments,
        best_connection_score,
        target_log2FC,
        target_padj
      ),
    
    by = "candidate"
  )


# ============================================================
# 8. STANDARDIZE CANDIDATE COORDINATES
# ============================================================

coord_candidate <- first_existing(
  coords,
  c(
    "candidate",
    "candidate_name",
    "Candidate"
  )
)

coord_start <- first_existing(
  coords,
  c(
    "candidate_start",
    "start",
    "genomic_start"
  )
)

coord_end <- first_existing(
  coords,
  c(
    "candidate_end",
    "end",
    "genomic_end"
  )
)

if (
  any(
    is.na(
      c(
        coord_candidate,
        coord_start,
        coord_end
      )
    )
  )
) {
  
  stop(
    "Could not identify candidate-coordinate columns."
  )
}


coords_std <- coords %>%
  
  transmute(
    
    raw_candidate =
      as.character(
        .data[[coord_candidate]]
      ),
    
    candidate_start =
      as.numeric(
        .data[[coord_start]]
      ),
    
    candidate_end =
      as.numeric(
        .data[[coord_end]]
      )
  ) %>%
  
  mutate(
    
    candidate = case_when(
      
      str_detect(
        raw_candidate,
        regex(
          "RS04120",
          ignore_case = TRUE
        )
      ) ~ "RS04120",
      
      str_detect(
        raw_candidate,
        regex(
          "RS05585",
          ignore_case = TRUE
        )
      ) ~ "RS05585",
      
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(
    !is.na(candidate)
  )


# ============================================================
# 9. STANDARDIZE GENOMIC NEIGHBORHOOD
# ============================================================

nb_candidate <- first_existing(
  neighborhood,
  c(
    "candidate",
    "candidate_name",
    "Candidate"
  )
)

nb_start <- first_existing(
  neighborhood,
  c(
    "genomic_start",
    "feature_start",
    "start"
  )
)

nb_end <- first_existing(
  neighborhood,
  c(
    "genomic_end",
    "feature_end",
    "end"
  )
)

nb_strand <- first_existing(
  neighborhood,
  c(
    "strand",
    "feature_strand"
  )
)

if (
  any(
    is.na(
      c(
        nb_candidate,
        nb_start,
        nb_end,
        nb_strand
      )
    )
  )
) {
  
  stop(
    "Could not identify genomic-neighborhood columns."
  )
}


# ============================================================
# 10. CREATE BEST GENE LABEL
# ============================================================

label_priority <- c(
  "gene_symbol",
  "display_label",
  "current_locus_tag",
  "gene_id"
)

label_priority <- label_priority[
  label_priority %in% names(neighborhood)
]


best_label <- function(df) {
  
  out <- rep(
    NA_character_,
    nrow(df)
  )
  
  for (nm in label_priority) {
    
    z <- as.character(
      df[[nm]]
    )
    
    usable <-
      (is.na(out) | out == "") &
      !is.na(z) &
      z != ""
    
    out[usable] <-
      z[usable]
  }
  
  out[is.na(out)] <- ""
  
  out
}


neighborhood_std <- neighborhood %>%
  
  mutate(
    
    candidate_raw =
      as.character(
        .data[[nb_candidate]]
      ),
    
    candidate = case_when(
      
      str_detect(
        candidate_raw,
        regex(
          "RS04120",
          ignore_case = TRUE
        )
      ) ~ "RS04120",
      
      str_detect(
        candidate_raw,
        regex(
          "RS05585",
          ignore_case = TRUE
        )
      ) ~ "RS05585",
      
      TRUE ~ NA_character_
    ),
    
    feature_start =
      as.numeric(
        .data[[nb_start]]
      ),
    
    feature_end =
      as.numeric(
        .data[[nb_end]]
      ),
    
    strand =
      as.character(
        .data[[nb_strand]]
      ),
    
    gene_label =
      best_label(neighborhood),
    
    gene_mid =
      (
        feature_start +
          feature_end
      ) / 2
  ) %>%
  
  filter(
    !is.na(candidate)
  )


# ============================================================
# 11. CLEAN GENOMIC LABELS
# ============================================================

important_symbols <- c(
  "nrdI",
  "nrdE",
  "nrdF",
  "sstA",
  "sstB",
  "sstC",
  "qoxC",
  "qoxB",
  "qoxA",
  "folD",
  "purE",
  "purK",
  "purC",
  "purQ",
  "purL"
)


neighborhood_std <- neighborhood_std %>%
  
  mutate(
    
    gene_label_clean = case_when(
      
      gene_label %in%
        important_symbols ~
        gene_label,
      
      str_detect(
        gene_label,
        "^SAA6008_|^Saa-6008"
      ) ~ "",
      
      TRUE ~
        gene_label
    )
  )


# ============================================================
# 12. GENOMIC PANEL FUNCTION
# ============================================================

genomic_panel <- function(
    candidate_id,
    panel_title,
    panel_subtitle,
    distance_label) {
  
  
  dat <- neighborhood_std %>%
    filter(
      candidate ==
        candidate_id
    )
  
  
  cc <- coords_std %>%
    filter(
      candidate ==
        candidate_id
    ) %>%
    slice(1)
  
  
  if (nrow(dat) == 0) {
    
    stop(
      paste(
        "No genomic neighborhood data for",
        candidate_id
      )
    )
  }
  
  
  if (nrow(cc) == 0) {
    
    stop(
      paste(
        "No candidate coordinate data for",
        candidate_id
      )
    )
  }
  
  
  candidate_mid <- mean(
    c(
      cc$candidate_start,
      cc$candidate_end
    )
  )
  
  
  dat <- dat %>%
    
    mutate(
      
      gene_y =
        ifelse(
          strand == "+",
          0.28,
          -0.28
        ),
      
      label_y =
        ifelse(
          strand == "+",
          0.49,
          -0.49
        )
    )
  
  
  xmin <- min(
    dat$feature_start,
    na.rm = TRUE
  )
  
  xmax <- max(
    dat$feature_end,
    na.rm = TRUE
  )
  
  
  ggplot() +
    
    # genomic baseline
    annotate(
      "segment",
      x = xmin,
      xend = xmax,
      y = 0,
      yend = 0,
      linewidth = 0.30
    ) +
    
    # genomic features
    geom_segment(
      data = dat,
      aes(
        x = feature_start,
        xend = feature_end,
        y = gene_y,
        yend = gene_y
      ),
      linewidth = 5.0,
      lineend = "butt"
    ) +
    
    # feature labels
    geom_text(
      data = dat %>%
        filter(
          gene_label_clean != ""
        ),
      aes(
        x = gene_mid,
        y = label_y,
        label = gene_label_clean
      ),
      size = 2.65,
      check_overlap = TRUE
    ) +
    
    # candidate position
    annotate(
      "segment",
      x = candidate_mid,
      xend = candidate_mid,
      y = 0.14,
      yend = 0.91,
      linewidth = 0.65
    ) +
    
    # candidate name
    annotate(
      "text",
      x = candidate_mid,
      y = 1.07,
      label = paste0(
        "3UTR-",
        candidate_id
      ),
      fontface = "bold",
      size = 3.15
    ) +
    
    # distance annotation
    annotate(
      "text",
      x = candidate_mid,
      y = 1.34,
      label = distance_label,
      size = 2.65
    ) +
    
    labs(
      title =
        panel_title,
      
      subtitle =
        panel_subtitle,
      
      x =
        "JKD6008 genomic coordinate (bp)",
      
      y = NULL
    ) +
    
    coord_cartesian(
      ylim = c(
        -0.72,
        1.52
      ),
      clip = "off"
    ) +
    
    theme_classic(
      base_size = 9
    ) +
    
    theme(
      
      axis.line.y =
        element_blank(),
      
      axis.ticks.y =
        element_blank(),
      
      axis.text.y =
        element_blank(),
      
      axis.title.y =
        element_blank(),
      
      axis.text.x =
        element_text(
          size = 7.4
        ),
      
      axis.title.x =
        element_text(
          size = 8
        ),
      
      plot.title =
        element_text(
          face = "bold",
          size = 10.3
        ),
      
      plot.subtitle =
        element_text(
          size = 7.8,
          margin = margin(
            b = 5
          )
        ),
      
      plot.margin =
        margin(
          6,
          14,
          6,
          14
        )
    )
}


# ============================================================
# 13. PANELS A AND B
# ============================================================

panel_A <- genomic_panel(
  
  candidate_id =
    "RS04120",
  
  panel_title =
    "A  Genomic context of RS04120",
  
  panel_subtitle =
    "Putative nrdF 3′UTR-associated RNA candidate",
  
  distance_label =
    "53 nt from the annotated nrdF 3′ end"
)


panel_B <- genomic_panel(
  
  candidate_id =
    "RS05585",
  
  panel_title =
    "B  Genomic context of RS05585",
  
  panel_subtitle =
    "Putative folD 3′UTR-associated RNA candidate; CLASH target = nrdF",
  
  distance_label =
    "177 nt from the annotated folD 3′ end"
)


# ============================================================
# 14. STANDARDIZE TERM-SEQ DATA
# ============================================================

term_candidate <- first_existing(
  termseq,
  c(
    "candidate",
    "candidate_name",
    "Candidate",
    "Candidate_ID"
  )
)

term_position <- first_existing(
  termseq,
  c(
    "relative_position",
    "relative_position_nt",
    "position_relative_to_boundary",
    "relative_nt",
    "position_relative",
    "relative_pos"
  )
)

term_signal <- first_existing(
  termseq,
  c(
    "signal",
    "termseq_signal",
    "term_seq_signal",
    "three_prime_signal",
    "value",
    "count"
  )
)

term_rep <- first_existing(
  termseq,
  c(
    "replicate",
    "Replicate",
    "sample",
    "sample_id"
  )
)


if (is.na(term_candidate)) {
  
  term_candidate <-
    names(termseq)[
      str_detect(
        names(termseq),
        regex(
          "candidate",
          ignore_case = TRUE
        )
      )
    ][1]
}


if (is.na(term_position)) {
  
  term_position <-
    names(termseq)[
      str_detect(
        names(termseq),
        regex(
          "relative|position",
          ignore_case = TRUE
        )
      )
    ][1]
}


if (is.na(term_signal)) {
  
  term_signal <-
    names(termseq)[
      str_detect(
        names(termseq),
        regex(
          "signal|count|value",
          ignore_case = TRUE
        )
      )
    ][1]
}


if (is.na(term_rep)) {
  
  term_rep <-
    names(termseq)[
      str_detect(
        names(termseq),
        regex(
          "rep|sample",
          ignore_case = TRUE
        )
      )
    ][1]
}


if (
  any(
    is.na(
      c(
        term_candidate,
        term_position,
        term_signal,
        term_rep
      )
    )
  )
) {
  
  stop(
    "Could not identify required Term-seq columns."
  )
}


termseq_std <- termseq %>%
  
  transmute(
    
    candidate_raw =
      as.character(
        .data[[term_candidate]]
      ),
    
    position =
      as.numeric(
        .data[[term_position]]
      ),
    
    signal =
      as.numeric(
        .data[[term_signal]]
      ),
    
    replicate_raw =
      as.character(
        .data[[term_rep]]
      )
  ) %>%
  
  mutate(
    
    candidate = case_when(
      
      str_detect(
        candidate_raw,
        regex(
          "RS04120",
          ignore_case = TRUE
        )
      ) ~ "RS04120",
      
      str_detect(
        candidate_raw,
        regex(
          "RS05585",
          ignore_case = TRUE
        )
      ) ~ "RS05585",
      
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(
    !is.na(candidate),
    position >= -100,
    position <= 100
  )


# ============================================================
# 15. STANDARDIZE REPLICATE NAMES
# ============================================================

termseq_std <- termseq_std %>%
  
  group_by(candidate) %>%
  
  mutate(
    
    replicate =
      paste0(
        "rep",
        match(
          replicate_raw,
          unique(
            replicate_raw
          )
        )
      )
  ) %>%
  
  ungroup()


# ============================================================
# 16. FINAL TERM-SEQ PANEL FUNCTION
#
# IMPORTANT FIX:
# NO TEXT IS PLACED INSIDE THE DATA FIELD.
#
# x = 0 is represented only by:
#   - dashed vertical line
#   - ±5 nt shaded region
#
# Meaning is explained in subtitle.
# ============================================================

termseq_panel <- function(
    candidate_id,
    panel_title) {
  
  
  dat <- termseq_std %>%
    
    filter(
      candidate ==
        candidate_id
    )
  
  
  if (nrow(dat) == 0) {
    
    stop(
      paste(
        "No Term-seq data found for",
        candidate_id
      )
    )
  }
  
  
  ymax <- max(
    dat$signal,
    na.rm = TRUE
  )
  
  
  # modest headroom only
  upper_limit <-
    ymax * 1.08
  
  
  ggplot(
    dat,
    aes(
      x = position,
      y = signal,
      group = replicate,
      linetype = replicate
    )
  ) +
    
    # ±5 nt region
    annotate(
      "rect",
      xmin = -5,
      xmax = 5,
      ymin = 0,
      ymax = Inf,
      fill = "grey92",
      alpha = 0.60
    ) +
    
    # Term-seq signals
    geom_line(
      linewidth = 0.60
    ) +
    
    # CLASH-defined candidate boundary
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.60
    ) +
    
    scale_linetype_manual(
      values = c(
        rep1 = "solid",
        rep2 = "dashed"
      ),
      name = "Replicate"
    ) +
    
    scale_x_continuous(
      
      breaks = c(
        -100,
        -50,
        0,
        50,
        100
      ),
      
      limits = c(
        -105,
        105
      )
    ) +
    
    scale_y_continuous(
      
      limits = c(
        0,
        upper_limit
      ),
      
      expand = expansion(
        mult = c(
          0,
          0.02
        )
      )
    ) +
    
    labs(
      
      title =
        panel_title,
      
      subtitle =
        paste0(
          "Term-seq 3′-end signal within ±100 nt; ",
          "dashed vertical line = candidate boundary; ",
          "shading = ±5 nt"
        ),
      
      x =
        "Position relative to candidate boundary (nt)",
      
      y =
        "Term-seq 3′-end signal"
    ) +
    
    theme_classic(
      base_size = 9
    ) +
    
    theme(
      
      legend.position =
        "top",
      
      legend.justification =
        "center",
      
      legend.direction =
        "horizontal",
      
      legend.title =
        element_text(
          size = 7.4
        ),
      
      legend.text =
        element_text(
          size = 7.4
        ),
      
      legend.key.width =
        unit(
          0.65,
          "cm"
        ),
      
      plot.title =
        element_text(
          face = "bold",
          size = 10.3
        ),
      
      plot.subtitle =
        element_text(
          size = 7.5,
          margin = margin(
            b = 5
          )
        ),
      
      axis.title =
        element_text(
          size = 8
        ),
      
      axis.text =
        element_text(
          size = 7.4
        ),
      
      plot.margin =
        margin(
          6,
          14,
          6,
          14
        )
    )
}


# ============================================================
# 17. PANELS C AND D
# ============================================================

panel_C <- termseq_panel(
  
  candidate_id =
    "RS04120",
  
  panel_title =
    "C  Term-seq profile around RS04120"
)


panel_D <- termseq_panel(
  
  candidate_id =
    "RS05585",
  
  panel_title =
    "D  Term-seq profile around RS05585"
)


# ============================================================
# 18. PREPARE EVIDENCE SUMMARY
# ============================================================

evidence_one <- evidence %>%
  
  group_by(candidate) %>%
  
  slice(1) %>%
  
  ungroup()


get_e <- function(id, col) {
  
  x <- evidence_one %>%
    
    filter(
      candidate == id
    )
  
  x[[col]][1]
}


# ============================================================
# 19. MANUSCRIPT-READY EVIDENCE TABLE
# ============================================================

summary_table <- tibble(
  
  Evidence = c(
    
    "Putative parent",
    
    "Distance from parent 3′ end",
    
    "RNase III-CLASH hybrids",
    
    "Independent CLASH experiments",
    
    "CLASH connection score",
    
    "nrdF target log2FC",
    
    "nrdF adjusted P",
    
    "Sharp reproducible Term-seq boundary"
  ),
  
  
  RS04120 = c(
    
    "nrdF",
    
    "53 nt",
    
    as.character(
      get_e(
        "RS04120",
        "total_hybrid_count"
      )
    ),
    
    as.character(
      get_e(
        "RS04120",
        "maximum_number_of_experiments"
      )
    ),
    
    sprintf(
      "%.3f",
      get_e(
        "RS04120",
        "best_connection_score"
      )
    ),
    
    sprintf(
      "%.3f",
      get_e(
        "RS04120",
        "target_log2FC"
      )
    ),
    
    "3.04 × 10⁻¹¹",
    
    "No"
  ),
  
  
  RS05585 = c(
    
    "folD",
    
    "177 nt",
    
    as.character(
      get_e(
        "RS05585",
        "total_hybrid_count"
      )
    ),
    
    as.character(
      get_e(
        "RS05585",
        "maximum_number_of_experiments"
      )
    ),
    
    sprintf(
      "%.3f",
      get_e(
        "RS05585",
        "best_connection_score"
      )
    ),
    
    sprintf(
      "%.3f",
      get_e(
        "RS05585",
        "target_log2FC"
      )
    ),
    
    "3.04 × 10⁻¹¹",
    
    "No"
  )
)


# ============================================================
# 20. PANEL E DATA
# ============================================================

table_long <- summary_table %>%
  
  mutate(
    row_id =
      row_number()
  ) %>%
  
  pivot_longer(
    
    cols = c(
      RS04120,
      RS05585
    ),
    
    names_to =
      "candidate",
    
    values_to =
      "value"
  ) %>%
  
  mutate(
    
    Evidence =
      factor(
        Evidence,
        levels =
          rev(
            summary_table$Evidence
          )
      ),
    
    candidate =
      factor(
        candidate,
        levels = c(
          "RS04120",
          "RS05585"
        )
      ),
    
    alternate =
      row_id %% 2 == 0
  )


# ============================================================
# 21. PANEL E
# ============================================================

panel_E <- ggplot(
  table_long,
  aes(
    x = candidate,
    y = Evidence
  )
) +
  
  geom_tile(
    
    aes(
      alpha = alternate
    ),
    
    fill = "grey88",
    
    linewidth = 0.20
  ) +
  
  scale_alpha_manual(
    
    values = c(
      "FALSE" = 0.10,
      "TRUE" = 0.30
    ),
    
    guide = "none"
  ) +
  
  geom_text(
    
    aes(
      label = value
    ),
    
    size = 3.0
  ) +
  
  scale_x_discrete(
    position = "top"
  ) +
  
  labs(
    
    title =
      "E  Independent evidence layers",
    
    subtitle =
      paste0(
        "Evidence layers are reported separately; ",
        "no integrated biological-priority score was calculated"
      ),
    
    x = NULL,
    
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 9
  ) +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_text(
        face = "bold",
        size = 8.5,
        margin = margin(
          b = 5
        )
      ),
    
    axis.text.y =
      element_text(
        size = 7.9
      ),
    
    axis.ticks =
      element_blank(),
    
    plot.title =
      element_text(
        face = "bold",
        size = 10.3
      ),
    
    plot.subtitle =
      element_text(
        size = 7.7,
        margin = margin(
          b = 6
        )
      ),
    
    plot.margin =
      margin(
        5,
        30,
        5,
        30
      )
  )


# ============================================================
# 22. ASSEMBLE FIGURE
# ============================================================

top_row <-
  
  panel_A +
  
  panel_B +
  
  plot_layout(
    widths = c(
      1,
      1
    )
  )


middle_row <-
  
  panel_C +
  
  panel_D +
  
  plot_layout(
    widths = c(
      1,
      1
    )
  )


figure3_final <-
  
  top_row /
  
  middle_row /
  
  panel_E +
  
  plot_layout(
    
    heights = c(
      0.88,
      1.20,
      0.94
    )
  ) +
  
  plot_annotation(
    
    title =
      paste0(
        "Genomic context and independent evidence layers for two ",
        "putative 3′UTR-associated RNA candidates"
      ),
    
    subtitle =
      paste0(
        "Genomic proximity, RNase III-CLASH interaction evidence, ",
        "target transcriptomic response and Term-seq boundary evidence ",
        "are evaluated independently."
      ),
    
    theme =
      theme(
        
        plot.title =
          element_text(
            face = "bold",
            size = 15,
            margin = margin(
              b = 4
            )
          ),
        
        plot.subtitle =
          element_text(
            size = 9.4,
            margin = margin(
              b = 8
            )
          ),
        
        plot.margin =
          margin(
            12,
            18,
            12,
            18
          )
      )
  )


# ============================================================
# 23. DISPLAY
# ============================================================

print(
  figure3_final
)


# ============================================================
# 24. OUTPUT FILE NAMES
# ============================================================

pdf_file <- file.path(
  figdir,
  "Figure3_FINAL_CORRECTED.pdf"
)

png_file <- file.path(
  figdir,
  "Figure3_FINAL_CORRECTED_600dpi.png"
)

tiff_file <- file.path(
  figdir,
  "Figure3_FINAL_CORRECTED_600dpi.tiff"
)


# ============================================================
# 25. EXPORT PDF
# ============================================================

ggsave(
  
  filename =
    pdf_file,
  
  plot =
    figure3_final,
  
  width =
    13.5,
  
  height =
    12.5,
  
  units =
    "in",
  
  device =
    cairo_pdf,
  
  limitsize =
    FALSE
)


# ============================================================
# 26. EXPORT PNG
# ============================================================

ggsave(
  
  filename =
    png_file,
  
  plot =
    figure3_final,
  
  width =
    13.5,
  
  height =
    12.5,
  
  units =
    "in",
  
  dpi =
    600,
  
  limitsize =
    FALSE
)


# ============================================================
# 27. EXPORT TIFF
# ============================================================

ggsave(
  
  filename =
    tiff_file,
  
  plot =
    figure3_final,
  
  width =
    13.5,
  
  height =
    12.5,
  
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
# 28. SAVE DATA USED IN FIGURE
# ============================================================

write_csv(
  
  summary_table,
  
  file.path(
    tabdir,
    "Figure3_Final_Evidence_Summary.csv"
  )
)


write_csv(
  
  termseq_std,
  
  file.path(
    tabdir,
    "Figure3_Final_TermSeq_Data.csv"
  )
)


write_csv(
  
  neighborhood_std,
  
  file.path(
    tabdir,
    "Figure3_Final_Genomic_Context.csv"
  )
)


# ============================================================
# 29. FINAL FIGURE LEGEND
# ============================================================

legend_text <- paste0(
  
  "Figure 3. Genomic context and independent evidence layers for two ",
  "putative 3′UTR-associated RNA candidates. ",
  
  "(A) The RS04120-labelled CLASH RNA fragment is located 53 nt from ",
  "the annotated nrdF 3′ end on the same strand, supporting its ",
  "description as a putative nrdF 3′UTR-associated RNA candidate. ",
  
  "(B) The RS05585-labelled CLASH RNA fragment is located 177 nt from ",
  "the annotated folD 3′ end on the same strand, consistent with a ",
  "putative extended folD 3′UTR-associated RNA candidate; the ",
  "RNase III-CLASH-associated mRNA target considered here is nrdF. ",
  
  "(C-D) Term-seq 3′-end signal within ±100 nt of each CLASH-defined ",
  "candidate boundary. The dashed vertical line denotes the candidate ",
  "boundary and the shaded region denotes ±5 nt around that position. ",
  "Neither candidate showed a sharp, reproducible Term-seq boundary ",
  "across replicates at the CLASH-defined fragment boundary. ",
  
  "(E) Independent evidence summary. The RS04120-nrdF interaction was ",
  "represented by 15 CLASH hybrids across two experiments with a ",
  "connection score of 0.289, whereas RS05585-nrdF was represented by ",
  "seven hybrids in one experiment with a connection score of 0.106. ",
  "Independently, nrdF transcript abundance decreased under vancomycin ",
  "exposure (log2FC = -1.118; adjusted P = 3.04 × 10^-11). ",
  
  "Genomic context, CLASH interaction evidence, target transcriptomic ",
  "response and Term-seq boundary evidence were evaluated as separate ",
  "evidence layers. These observations support the description of the ",
  "fragments as putative 3′UTR-associated RNA candidates but do not ",
  "establish independent RNA processing, direct regulation of nrdF, ",
  "or a causal role in the vancomycin response."
)


writeLines(
  
  legend_text,
  
  file.path(
    outdir,
    "Figure3_FINAL_Legend.txt"
  )
)


# ============================================================
# 30. SAVE METHODOLOGICAL NOTE
# ============================================================

method_note <- c(
  
  "FIGURE 3 — INTERPRETATION RULES",
  "================================",
  "",
  
  "RS04120:",
  "- Candidate fragment lies 53 nt from the annotated nrdF 3′ end.",
  "- Same-strand genomic context supports 3′UTR association.",
  "- RNase III-CLASH interaction with nrdF: 15 hybrids, 2 experiments.",
  "- Connection score: 0.289.",
  "- nrdF transcript response: log2FC = -1.118; adjusted P = 3.04e-11.",
  "- Term-seq does not demonstrate a sharp reproducible boundary.",
  "",
  
  "RS05585:",
  "- Candidate fragment lies 177 nt from the annotated folD 3′ end.",
  "- Same-strand genomic context supports a putative extended folD 3′UTR association.",
  "- RNase III-CLASH interaction with nrdF: 7 hybrids, 1 experiment.",
  "- Connection score: 0.106.",
  "- nrdF transcript response: log2FC = -1.118; adjusted P = 3.04e-11.",
  "- Term-seq does not demonstrate a sharp reproducible boundary.",
  "",
  
  "Preferred terminology:",
  "putative 3′UTR-associated RNA candidate",
  "",
  
  "Do NOT claim:",
  "- independently processed RNA",
  "- direct regulation of nrdF",
  "- causal role in vancomycin response",
  "- experimentally validated regulatory axis",
  "",
  
  "Evidence layers remain analytically separate:",
  "- genomic context",
  "- RNase III-CLASH interaction evidence",
  "- target transcriptomic response",
  "- descriptive Term-seq boundary evidence",
  "",
  
  "No integrated biological-priority score was used."
)


writeLines(
  
  method_note,
  
  file.path(
    outdir,
    "Figure3_Methodological_Notes.txt"
  )
)


# ============================================================
# 31. SESSION INFO
# ============================================================

capture.output(
  
  sessionInfo(),
  
  file =
    file.path(
      outdir,
      "sessionInfo.txt"
    )
)


# ============================================================
# 32. FINAL CONSOLE REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FIGURE 3 FINAL CORRECTED VERSION COMPLETED\n")
cat("============================================================\n\n")

cat("Main corrections:\n")
cat("1. No annotation text is placed inside Term-seq peak fields.\n")
cat("2. Dashed line at x = 0 identifies the candidate boundary.\n")
cat("3. Shaded region identifies ±5 nt around the boundary.\n")
cat("4. Boundary meaning is explained in the panel subtitle.\n")
cat("5. Adjusted P is displayed as 3.04 × 10⁻¹¹ in Panel E.\n")
cat("6. CLASH and transcriptomic evidence remain separate.\n")
cat("7. No processing/direct-regulation claim is made.\n\n")

cat("PNG:\n")
cat(png_file, "\n\n")

cat("PDF:\n")
cat(pdf_file, "\n\n")

cat("TIFF:\n")
cat(tiff_file, "\n\n")

cat("Legend:\n")
cat(
  file.path(
    outdir,
    "Figure3_FINAL_Legend.txt"
  ),
  "\n\n"
)

cat("============================================================\n")
project <- "D:/Bac-sRNA"

files <- list.files(
  project,
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

intarna_files <- files[
  grepl(
    "intarna|interaction.*energy|binding.*site|duplex|seed",
    basename(files),
    ignore.case = TRUE
  )
]

cat("\n===== POSSIBLE IntaRNA FILES =====\n")
cat(intarna_files, sep = "\n")
# ============================================================
# AUDIT THE TWO FINAL IntaRNA TABLES BEFORE REVISED FIGURE 4
# ============================================================

f1 <- "D:/Bac-sRNA/19_IntaRNA_integration_validated_coloured/tables/Manuscript_IntaRNA_structural_summary.csv"

f2 <- "D:/Bac-sRNA/29_IntaRNA_structural_evidence/Table5_IntaRNA_structural_characteristics.csv"

cat("\n====================================================\n")
cat("FILE 1 — VALIDATED COLOURED PIPELINE\n")
cat("====================================================\n")

x1 <- read.csv(f1, check.names = FALSE)

cat("\nDimensions:\n")
print(dim(x1))

cat("\nColumn names:\n")
print(names(x1))

cat("\nData:\n")
print(x1)


cat("\n\n====================================================\n")
cat("FILE 2 — LATER SCRIPT 29 PIPELINE\n")
cat("====================================================\n")

x2 <- read.csv(f2, check.names = FALSE)

cat("\nDimensions:\n")
print(dim(x2))

cat("\nColumn names:\n")
print(names(x2))

cat("\nData:\n")
print(x2)


cat("\n\n====================================================\n")
cat("RS04120 / RS05585 CHECK\n")
cat("====================================================\n")

cat("\nFILE 1:\n")
print(
  x1[
    apply(
      x1,
      1,
      function(z)
        any(
          grepl(
            "RS04120|RS05585|nrdF",
            z,
            ignore.case = TRUE
          )
        )
    ),
    ,
    drop = FALSE
  ]
)

cat("\nFILE 2:\n")
print(
  x2[
    apply(
      x2,
      1,
      function(z)
        any(
          grepl(
            "RS04120|RS05585|nrdF",
            z,
            ignore.case = TRUE
          )
        )
    ),
    ,
    drop = FALSE
  ]
)