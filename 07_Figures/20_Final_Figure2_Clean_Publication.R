# ============================================================
# SCRIPT 20 — FINAL REPLACEMENT
# Publication-ready Figure 2
#
# Evidence-layered reconstruction of the
# vancomycin-associated RNase III-CLASH network
#
# A. Independent evidence-layer architecture
# B. Vancomycin-responsive target subset
# C. Benchmark-associated interactions and nrdF candidates
# D. Target transcriptomic evidence across complete network
#
# ANALYTICAL PRINCIPLES
# ------------------------------------------------------------
# 1. RNase III-CLASH defines the interaction network.
# 2. Target differential expression is NOT an inclusion criterion.
# 3. The 100 Higher/Intermediate CLASH interactions and the
#    10 interactions with significant target DEGs are separate
#    summaries derived independently from the 211-edge network.
# 4. Target transcriptomic response is an independent layer.
# 5. No integrated biological-priority score is used.
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
  "scales",
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
library(scales)
library(grid)


# ============================================================
# 02. PROJECT PATHS
# ============================================================

project <- "D:/Bac-sRNA"

input_dir <- file.path(
  project,
  "14_revised_evidence_layered_network",
  "tables"
)

outdir <- file.path(
  project,
  "20_final_Figure2_clean_publication"
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
# 03. INPUT FILES
# ============================================================

all_file <- file.path(
  input_dir,
  "Evidence_Layered_CLASH_Network_ALL.csv"
)

responsive_file <- file.path(
  input_dir,
  "Vancomycin_Responsive_Target_Subset.csv"
)

benchmark_file <- file.path(
  input_dir,
  "Benchmark_vs_nrdF_Descriptive_Comparison.csv"
)


# ============================================================
# 04. INPUT CHECK
# ============================================================

input_check <- tibble(
  input = c(
    "Complete CLASH network",
    "Responsive target subset",
    "Benchmark/nrdF comparison"
  ),
  
  path = c(
    all_file,
    responsive_file,
    benchmark_file
  ),
  
  exists = c(
    file.exists(all_file),
    file.exists(responsive_file),
    file.exists(benchmark_file)
  )
)

cat("\n")
cat("============================================\n")
cat("INPUT CHECK\n")
cat("============================================\n")

print(
  input_check,
  n = Inf,
  width = Inf
)

if (any(!input_check$exists)) {
  stop("One or more required Script 14 files are missing.")
}


# ============================================================
# 05. READ DATA
# ============================================================

net <- read_csv(
  all_file,
  show_col_types = FALSE
)

responsive <- read_csv(
  responsive_file,
  show_col_types = FALSE
)

benchmark <- read_csv(
  benchmark_file,
  show_col_types = FALSE
)


# ============================================================
# 06. CORE COUNTS
# ============================================================

n_total <- nrow(net)


n_hi_int <- net %>%
  
  filter(
    str_detect(
      as.character(
        clash_support_category
      ),
      regex(
        "Higher|Intermediate",
        ignore_case = TRUE
      )
    )
  ) %>%
  
  nrow()


n_deg_edges <- nrow(
  responsive
)


n_deg_targets <- responsive %>%
  
  summarise(
    n = n_distinct(
      target_label
    )
  ) %>%
  
  pull(n)


cat("\n")
cat("============================================\n")
cat("CORE COUNTS\n")
cat("============================================\n")

cat(
  "Complete CLASH network: ",
  n_total,
  "\n",
  sep = ""
)

cat(
  "Higher/Intermediate CLASH support: ",
  n_hi_int,
  "\n",
  sep = ""
)

cat(
  "Interactions with significant target DEGs: ",
  n_deg_edges,
  "\n",
  sep = ""
)

cat(
  "Unique significant target genes: ",
  n_deg_targets,
  "\n",
  sep = ""
)


if (n_total != 211) {
  
  warning(
    paste(
      "Expected 211 complete CLASH interactions;",
      "observed",
      n_total
    )
  )
}


if (n_hi_int != 100) {
  
  warning(
    paste(
      "Expected 100 Higher/Intermediate interactions;",
      "observed",
      n_hi_int
    )
  )
}


if (n_deg_edges != 10) {
  
  warning(
    paste(
      "Expected 10 DEG-associated interactions;",
      "observed",
      n_deg_edges
    )
  )
}


# ============================================================
# 07. CLEAN DISPLAY LABELS
# ============================================================

clean_source_label <- function(x) {
  
  x <- as.character(x)
  
  case_when(
    
    str_detect(
      x,
      regex(
        "RS04120",
        ignore_case = TRUE
      )
    ) ~ "RS04120",
    
    str_detect(
      x,
      regex(
        "RS05585",
        ignore_case = TRUE
      )
    ) ~ "RS05585",
    
    str_detect(
      x,
      regex(
        "RsaOI|Sau-6477",
        ignore_case = TRUE
      )
    ) ~ "RsaOI",
    
    str_detect(
      x,
      regex(
        "SprA2/SprAs2",
        ignore_case = TRUE
      )
    ) ~ "SprA2/SprAs2",
    
    str_detect(
      x,
      regex(
        "SprA2/RsaJ",
        ignore_case = TRUE
      )
    ) ~ "SprA2/RsaJ",
    
    TRUE ~ x
  )
}


clean_target_label <- function(x) {
  
  x <- as.character(x)
  
  case_when(
    
    str_detect(
      x,
      regex(
        "^nrdF$",
        ignore_case = TRUE
      )
    ) ~ "nrdF",
    
    str_detect(
      x,
      regex(
        "^guaB$",
        ignore_case = TRUE
      )
    ) ~ "guaB",
    
    str_detect(
      x,
      regex(
        "^spa$",
        ignore_case = TRUE
      )
    ) ~ "spa",
    
    str_detect(
      x,
      regex(
        "^qoxB$",
        ignore_case = TRUE
      )
    ) ~ "qoxB",
    
    str_detect(
      x,
      regex(
        "^gltB$",
        ignore_case = TRUE
      )
    ) ~ "gltB",
    
    str_detect(
      x,
      regex(
        "^dhaK$",
        ignore_case = TRUE
      )
    ) ~ "dhaK",
    
    str_detect(
      x,
      regex(
        "^modA$",
        ignore_case = TRUE
      )
    ) ~ "modA",
    
    TRUE ~ x
  )
}


# ============================================================
# 08. PANEL A
# CLEAN INDEPENDENT-BRANCH ARCHITECTURE
# ============================================================

panel_A <- ggplot() +
  
  # ----------------------------------------------------------
# LEFT BOX — COMPLETE NETWORK
# ----------------------------------------------------------

annotate(
  "rect",
  xmin = 0.55,
  xmax = 2.45,
  ymin = 1.35,
  ymax = 2.30,
  fill = "white",
  linewidth = 0.75
) +
  
  annotate(
    "text",
    x = 1.50,
    y = 2.03,
    label = as.character(
      n_total
    ),
    fontface = "bold",
    size = 5.4
  ) +
  
  annotate(
    "text",
    x = 1.50,
    y = 1.67,
    label = "Complete CLASH network",
    fontface = "bold",
    size = 3.9
  ) +
  
  
  # ----------------------------------------------------------
# UPPER ARROW
# ----------------------------------------------------------

annotate(
  "segment",
  x = 2.45,
  xend = 4.15,
  y = 2.05,
  yend = 2.68,
  linewidth = 0.75,
  arrow = arrow(
    length = unit(
      0.15,
      "cm"
    ),
    type = "closed"
  )
) +
  
  annotate(
    "text",
    x = 3.25,
    y = 2.80,
    label = "CLASH-support\nstratification",
    size = 2.9,
    lineheight = 1.05
  ) +
  
  
  # ----------------------------------------------------------
# UPPER BOX — 100
# ----------------------------------------------------------

annotate(
  "rect",
  xmin = 4.25,
  xmax = 6.65,
  ymin = 2.18,
  ymax = 3.30,
  fill = "white",
  linewidth = 0.75
) +
  
  annotate(
    "text",
    x = 5.45,
    y = 3.00,
    label = as.character(
      n_hi_int
    ),
    fontface = "bold",
    size = 5.5
  ) +
  
  annotate(
    "text",
    x = 5.45,
    y = 2.55,
    label =
      "Higher / Intermediate\nCLASH support",
    fontface = "bold",
    size = 3.75,
    lineheight = 1.08
  ) +
  
  
  # ----------------------------------------------------------
# LOWER ARROW
# ----------------------------------------------------------

annotate(
  "segment",
  x = 2.45,
  xend = 4.15,
  y = 1.62,
  yend = 0.97,
  linewidth = 0.75,
  arrow = arrow(
    length = unit(
      0.15,
      "cm"
    ),
    type = "closed"
  )
) +
  
  annotate(
    "text",
    x = 3.25,
    y = 0.83,
    label =
      "Independent target\ntranscriptomic layer",
    size = 2.9,
    lineheight = 1.05
  ) +
  
  
  # ----------------------------------------------------------
# LOWER BOX — 10
# ----------------------------------------------------------

annotate(
  "rect",
  xmin = 4.25,
  xmax = 6.65,
  ymin = 0.30,
  ymax = 1.45,
  fill = "white",
  linewidth = 0.75
) +
  
  annotate(
    "text",
    x = 5.45,
    y = 1.14,
    label = as.character(
      n_deg_edges
    ),
    fontface = "bold",
    size = 5.5
  ) +
  
  annotate(
    "text",
    x = 5.45,
    y = 0.68,
    label =
      "CLASH interactions with\nsignificant target DEGs",
    fontface = "bold",
    size = 3.55,
    lineheight = 1.08
  ) +
  
  
  # ----------------------------------------------------------
# KEY PRINCIPLE
# ----------------------------------------------------------

annotate(
  "text",
  x = 3.55,
  y = -0.18,
  label =
    paste(
      "Target differential expression was not required",
      "for inclusion in the CLASH network."
    ),
  fontface = "bold",
  size = 3.35
) +
  
  coord_cartesian(
    xlim = c(
      0.20,
      7.00
    ),
    ylim = c(
      -0.35,
      3.48
    ),
    clip = "off"
  ) +
  
  labs(
    title =
      "A  Independent evidence layers derived from the complete CLASH network",
    
    subtitle =
      paste(
        "CLASH-support stratification and target transcriptomic",
        "response are evaluated separately"
      )
  ) +
  
  theme_void(
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
        size = 8.8
      ),
    
    plot.margin =
      margin(
        8,
        15,
        10,
        15
      )
  )


# ============================================================
# 09. PANEL B DATA
# ============================================================

responsive_plot <- responsive %>%
  
  mutate(
    
    source_display =
      clean_source_label(
        source_label
      ),
    
    target_display =
      clean_target_label(
        target_label
      ),
    
    hybrid_count =
      as.numeric(
        total_hybrid_count
      )
  )


# ============================================================
# 10. PANEL B SOURCE ORDER
# ============================================================

preferred_source_order_B <- c(
  "INT-0-1153",
  "INT-0-1448",
  "INT-0-332",
  "INT-0-419",
  "INT-0-421",
  "RS04120",
  "RS05585",
  "SprA2/RsaJ",
  "sRNA-00053",
  "sRNA-00246"
)


existing_sources_B <- unique(
  responsive_plot$source_display
)


source_levels_B <- c(
  
  preferred_source_order_B[
    preferred_source_order_B %in%
      existing_sources_B
  ],
  
  setdiff(
    existing_sources_B,
    preferred_source_order_B
  )
)


source_positions_B <- tibble(
  
  source_display =
    source_levels_B,
  
  source_y =
    rev(
      seq_along(
        source_levels_B
      )
    )
)


# ============================================================
# 11. PANEL B TARGET ORDER
# ============================================================

preferred_target_order_B <- c(
  "dhaK",
  "gltB",
  "guaB",
  "modA",
  "nrdF",
  "qoxB",
  "spa"
)


existing_targets_B <- unique(
  responsive_plot$target_display
)


target_levels_B <- c(
  
  preferred_target_order_B[
    preferred_target_order_B %in%
      existing_targets_B
  ],
  
  setdiff(
    existing_targets_B,
    preferred_target_order_B
  )
)


target_positions_B <- tibble(
  
  target_display =
    target_levels_B,
  
  target_y =
    seq(
      from =
        length(
          source_levels_B
        ),
      
      to = 1,
      
      length.out =
        length(
          target_levels_B
        )
    )
)


responsive_plot <- responsive_plot %>%
  
  left_join(
    source_positions_B,
    by = "source_display"
  ) %>%
  
  left_join(
    target_positions_B,
    by = "target_display"
  )


# ============================================================
# 12. PANEL B EDGE APPEARANCE
# ============================================================

responsive_plot <- responsive_plot %>%
  
  mutate(
    
    edge_linetype =
      case_when(
        
        str_detect(
          as.character(
            clash_support_category
          ),
          regex(
            "Higher",
            ignore_case = TRUE
          )
        ) ~ "solid",
        
        str_detect(
          as.character(
            clash_support_category
          ),
          regex(
            "Intermediate",
            ignore_case = TRUE
          )
        ) ~ "dotted",
        
        TRUE ~ "dashed"
      ),
    
    edge_width =
      scales::rescale(
        hybrid_count,
        to = c(
          0.45,
          2.30
        )
      )
  )


# ============================================================
# 13. PANEL B
# ============================================================

panel_B <- ggplot() +
  
  geom_curve(
    data = responsive_plot,
    aes(
      x = 1.40,
      y = source_y,
      xend = 3.60,
      yend = target_y,
      linewidth = edge_width,
      linetype = edge_linetype
    ),
    curvature = 0.07,
    alpha = 0.78
  ) +
  
  scale_linewidth_identity() +
  
  scale_linetype_identity() +
  
  geom_point(
    data = source_positions_B,
    aes(
      x = 1.40,
      y = source_y
    ),
    shape = 22,
    size = 3.2,
    fill = "white",
    stroke = 0.85
  ) +
  
  geom_point(
    data = target_positions_B,
    aes(
      x = 3.60,
      y = target_y
    ),
    shape = 21,
    size = 3.2,
    fill = "white",
    stroke = 0.85
  ) +
  
  geom_text(
    data = source_positions_B,
    aes(
      x = 1.25,
      y = source_y,
      label = source_display
    ),
    hjust = 1,
    size = 3.0
  ) +
  
  geom_text(
    data = target_positions_B,
    aes(
      x = 3.75,
      y = target_y,
      label = target_display
    ),
    hjust = 0,
    size = 3.0,
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = 1.40,
    y =
      max(
        source_positions_B$source_y
      ) + 0.85,
    label =
      "Regulatory RNAs",
    fontface = "bold",
    size = 3.45
  ) +
  
  annotate(
    "text",
    x = 3.60,
    y =
      max(
        source_positions_B$source_y
      ) + 0.85,
    label =
      "Responsive mRNA targets",
    fontface = "bold",
    size = 3.45
  ) +
  
  coord_cartesian(
    xlim = c(
      0.40,
      4.55
    ),
    ylim = c(
      0.40,
      max(
        source_positions_B$source_y
      ) + 1.25
    ),
    clip = "off"
  ) +
  
  labs(
    title =
      "B  Vancomycin-responsive target subset",
    
    subtitle =
      paste0(
        n_deg_edges,
        " CLASH-supported interactions involving ",
        n_deg_targets,
        " significantly responsive mRNA targets"
      )
  ) +
  
  theme_void(
    base_size = 10
  ) +
  
  theme(
    
    legend.position =
      "none",
    
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
        42,
        8,
        42
      )
  )


# ============================================================
# 14. PANEL C DATA
# ============================================================

benchmark_plot <- benchmark %>%
  
  mutate(
    
    source_display =
      clean_source_label(
        source_label
      ),
    
    target_display =
      clean_target_label(
        target_label
      ),
    
    hybrid_count =
      as.numeric(
        total_hybrid_count
      )
  )


# ============================================================
# 15. PANEL C ORDER
# ============================================================

benchmark_plot <- benchmark_plot %>%
  
  mutate(
    
    order_key =
      case_when(
        
        source_display == "RsaOI" ~ 1,
        
        source_display == "SprA2/RsaJ" &
          target_display == "guaB" ~ 2,
        
        source_display == "SprA2/RsaJ" &
          target_display ==
          "WP_000482650.1" ~ 3,
        
        source_display ==
          "SprA2/SprAs2" ~ 4,
        
        source_display ==
          "RS04120" ~ 5,
        
        source_display ==
          "RS05585" ~ 6,
        
        TRUE ~ 99
      )
  ) %>%
  
  arrange(
    order_key
  ) %>%
  
  mutate(
    
    interaction_y =
      rev(
        seq_len(
          n()
        )
      ),
    
    edge_linetype =
      case_when(
        
        str_detect(
          as.character(
            clash_support_category
          ),
          regex(
            "Higher",
            ignore_case = TRUE
          )
        ) ~ "solid",
        
        str_detect(
          as.character(
            clash_support_category
          ),
          regex(
            "Intermediate",
            ignore_case = TRUE
          )
        ) ~ "dotted",
        
        TRUE ~ "dashed"
      ),
    
    edge_width =
      scales::rescale(
        hybrid_count,
        to = c(
          0.50,
          2.50
        )
      )
  )


# ============================================================
# 16. PANEL C
# ============================================================

panel_C <- ggplot() +
  
  geom_curve(
    data = benchmark_plot,
    aes(
      x = 1.40,
      y = interaction_y,
      xend = 3.60,
      yend = interaction_y,
      linewidth = edge_width,
      linetype = edge_linetype
    ),
    curvature = 0.06,
    alpha = 0.80
  ) +
  
  scale_linewidth_identity() +
  
  scale_linetype_identity() +
  
  geom_point(
    data = benchmark_plot,
    aes(
      x = 1.40,
      y = interaction_y
    ),
    shape = 22,
    size = 3.2,
    fill = "white",
    stroke = 0.85
  ) +
  
  geom_point(
    data = benchmark_plot,
    aes(
      x = 3.60,
      y = interaction_y
    ),
    shape = 21,
    size = 3.2,
    fill = "white",
    stroke = 0.85
  ) +
  
  geom_text(
    data = benchmark_plot,
    aes(
      x = 1.25,
      y = interaction_y,
      label = source_display
    ),
    hjust = 1,
    size = 3.0,
    fontface =
      ifelse(
        benchmark_plot$source_display %in%
          c(
            "RsaOI",
            "SprA2/RsaJ",
            "SprA2/SprAs2"
          ),
        "bold",
        "plain"
      )
  ) +
  
  geom_text(
    data = benchmark_plot,
    aes(
      x = 3.75,
      y = interaction_y,
      label = target_display
    ),
    hjust = 0,
    size = 3.0,
    fontface =
      ifelse(
        benchmark_plot$target_display ==
          "nrdF",
        "bold",
        "plain"
      )
  ) +
  
  annotate(
    "text",
    x = 1.40,
    y =
      max(
        benchmark_plot$interaction_y
      ) + 0.80,
    label =
      "Regulatory RNAs",
    fontface = "bold",
    size = 3.45
  ) +
  
  annotate(
    "text",
    x = 3.60,
    y =
      max(
        benchmark_plot$interaction_y
      ) + 0.80,
    label =
      "CLASH-associated mRNA targets",
    fontface = "bold",
    size = 3.45
  ) +
  
  coord_cartesian(
    xlim = c(
      0.35,
      4.65
    ),
    ylim = c(
      0.40,
      max(
        benchmark_plot$interaction_y
      ) + 1.20
    ),
    clip = "off"
  ) +
  
  labs(
    title =
      "C  Benchmark-associated interactions and nrdF candidates",
    
    subtitle =
      "Selected CLASH interactions are shown independently of target DEG status"
  ) +
  
  theme_void(
    base_size = 10
  ) +
  
  theme(
    
    legend.position =
      "none",
    
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
        45,
        8,
        45
      )
  )


# ============================================================
# 17. PANEL D DATA
# ============================================================

transcriptomic_distribution <- net %>%
  
  mutate(
    
    transcript_category =
      case_when(
        
        target_is_significant_deg %in%
          TRUE ~
          
          "Significant target DEG",
        
        str_detect(
          as.character(
            target_transcriptomic_layer
          ),
          regex(
            "measured|below",
            ignore_case = TRUE
          )
        ) ~
          
          "Measured, below DEG threshold",
        
        TRUE ~
          
          "No mapped target transcriptomic evidence"
      )
  ) %>%
  
  count(
    transcript_category,
    name = "n"
  )


desired_categories <- tibble(
  
  transcript_category = c(
    "Significant target DEG",
    "Measured, below DEG threshold",
    "No mapped target transcriptomic evidence"
  )
)


transcriptomic_distribution <- desired_categories %>%
  
  left_join(
    transcriptomic_distribution,
    by = "transcript_category"
  ) %>%
  
  mutate(
    
    n =
      replace_na(
        n,
        0
      ),
    
    percentage =
      100 *
      n /
      sum(n),
    
    transcript_category =
      factor(
        transcript_category,
        levels = c(
          "Significant target DEG",
          "Measured, below DEG threshold",
          "No mapped target transcriptomic evidence"
        )
      )
  )


# ============================================================
# 18. PANEL D
# ============================================================

max_bar <- max(
  transcriptomic_distribution$n
)


panel_D <- ggplot(
  transcriptomic_distribution,
  aes(
    x = transcript_category,
    y = n
  )
) +
  
  geom_col(
    width = 0.60,
    fill = "grey82",
    linewidth = 0.6
  ) +
  
  geom_text(
    aes(
      label =
        paste0(
          n,
          "\n(",
          sprintf(
            "%.1f",
            percentage
          ),
          "%)"
        )
    ),
    vjust = -0.55,
    fontface = "bold",
    size = 3.7
  ) +
  
  scale_y_continuous(
    
    limits = c(
      0,
      max_bar * 1.30
    ),
    
    breaks =
      pretty(
        c(
          0,
          max_bar * 1.20
        ),
        n = 5
      ),
    
    expand = c(
      0,
      0
    )
  ) +
  
  labs(
    title =
      "D  Target transcriptomic evidence across the complete CLASH network",
    
    subtitle =
      paste0(
        "All ",
        n_total,
        " CLASH-supported interactions are retained regardless of target DEG status"
      ),
    
    x = NULL,
    
    y =
      "Number of CLASH-supported interactions"
  ) +
  
  theme_classic(
    base_size = 10
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        size = 8.7
      ),
    
    axis.title.y =
      element_text(
        size = 9
      ),
    
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
        25,
        12,
        25
      )
  )


# ============================================================
# 19. COMBINE FINAL FIGURE
# ============================================================

figure2 <- (
  
  panel_A /
    
    (
      panel_B |
        panel_C
    ) /
    
    panel_D
  
) +
  
  plot_layout(
    
    heights = c(
      1.05,
      1.55,
      1.10
    )
  ) +
  
  plot_annotation(
    
    title =
      paste(
        "Evidence-layered reconstruction of the",
        "vancomycin-associated RNase III-CLASH network"
      ),
    
    subtitle =
      paste(
        "Physical RNA-RNA interaction evidence and target transcriptomic",
        "response are evaluated as analytically distinct evidence layers."
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
          ),
        
        plot.margin =
          margin(
            12,
            15,
            12,
            15
          )
      )
  )


# ============================================================
# 20. DISPLAY
# ============================================================

print(
  figure2
)


# ============================================================
# 21. EXPORT PDF
# ============================================================

ggsave(
  
  filename = file.path(
    figdir,
    "Figure2_FINAL_CLEAN_v2.pdf"
  ),
  
  plot = figure2,
  
  width = 13.5,
  
  height = 12.2,
  
  units = "in",
  
  device = cairo_pdf,
  
  limitsize = FALSE
)


# ============================================================
# 22. EXPORT PNG — 600 DPI
# ============================================================

ggsave(
  
  filename = file.path(
    figdir,
    "Figure2_FINAL_CLEAN_v2_600dpi.png"
  ),
  
  plot = figure2,
  
  width = 13.5,
  
  height = 12.2,
  
  units = "in",
  
  dpi = 600,
  
  limitsize = FALSE
)


# ============================================================
# 23. EXPORT TIFF — 600 DPI
# ============================================================

ggsave(
  
  filename = file.path(
    figdir,
    "Figure2_FINAL_CLEAN_v2_600dpi.tiff"
  ),
  
  plot = figure2,
  
  width = 13.5,
  
  height = 12.2,
  
  units = "in",
  
  dpi = 600,
  
  compression = "lzw",
  
  limitsize = FALSE
)


# ============================================================
# 24. SAVE FIGURE SOURCE TABLES
# ============================================================

write_csv(
  
  responsive_plot,
  
  file.path(
    tabdir,
    "Figure2B_Responsive_Target_Subset_v2.csv"
  )
)


write_csv(
  
  benchmark_plot,
  
  file.path(
    tabdir,
    "Figure2C_Benchmark_nrdF_Interactions_v2.csv"
  )
)


write_csv(
  
  transcriptomic_distribution,
  
  file.path(
    tabdir,
    "Figure2D_Transcriptomic_Evidence_Distribution_v2.csv"
  )
)


# ============================================================
# 25. FINAL FIGURE LEGEND
# ============================================================

figure_legend <- paste(
  
  "Figure 2. Evidence-layered reconstruction of the",
  "vancomycin-associated RNase III-CLASH network.",
  
  "(A) The complete collapsed RNase III-CLASH network comprised",
  "211 sRNA-mRNA interactions. Two analytically independent summaries",
  "were derived from this interaction universe: 100 interactions",
  "classified as having higher or intermediate descriptive CLASH",
  "support and 10 interactions involving mRNA targets that met the",
  "prespecified differential-expression criterion. Target differential",
  "expression was not required for inclusion in the CLASH network.",
  
  "(B) The vancomycin-responsive target subset comprised 10",
  "CLASH-supported interactions involving seven significantly",
  "responsive mRNA targets.",
  
  "(C) Selected benchmark-associated CLASH interactions and the two",
  "putative 3′UTR-associated nrdF candidate interactions. RsaOI was",
  "retained in the complete CLASH network independently of target",
  "differential-expression status. This panel does not imply recovery",
  "of all previously reported RsaOI targets.",
  
  "(D) Across the complete network, 10 interactions involved",
  "significant target DEGs, 71 involved targets that were measured",
  "but remained below the prespecified DEG threshold, and 130 lacked",
  "mapped target transcriptomic evidence in the integrated analysis.",
  
  "In panels B and C, edge width represents total CLASH hybrid count,",
  "whereas line type denotes the descriptive CLASH-support category.",
  "Physical interaction evidence and target transcriptional response",
  "were therefore retained as analytically distinct evidence layers."
)


writeLines(
  
  figure_legend,
  
  file.path(
    outdir,
    "Figure2_FINAL_CLEAN_v2_Legend.txt"
  )
)


# ============================================================
# 26. MANUSCRIPT RESULTS PARAGRAPH
# ============================================================

results_text <- paste(
  
  "RNase III-CLASH and target transcriptomic response were analyzed",
  "as independent evidence layers to avoid restricting the interaction",
  "network to targets showing changes in steady-state mRNA abundance.",
  "The complete collapsed CLASH network comprised 211 sRNA-mRNA",
  "interactions. Of these, 100 were classified as having higher or",
  "intermediate descriptive CLASH support. Independently, 10",
  "CLASH-supported interactions involved seven mRNA targets meeting",
  "the prespecified vancomycin-responsive differential-expression",
  "criterion. Across the complete network, 71 additional interactions",
  "involved targets that were measured but remained below the",
  "differential-expression threshold, whereas 130 lacked mapped target",
  "transcriptomic evidence in the integrated dataset. RsaOI-associated",
  "CLASH evidence was retained independently of target differential-",
  "expression status. Thus, absence from the responsive-target subset",
  "was not interpreted as absence of an RNA-RNA interaction."
)


writeLines(
  
  results_text,
  
  file.path(
    outdir,
    "Figure2_FINAL_CLEAN_v2_Results.txt"
  )
)


# ============================================================
# 27. REVIEWER RESPONSE TEXT
# ============================================================

reviewer_text <- paste(
  
  "In the revised analysis, target mRNA differential expression is",
  "no longer used as a prerequisite for retaining RNase III-CLASH-",
  "supported interactions. The complete CLASH network is analyzed",
  "first, while target transcriptional response is overlaid as an",
  "independent evidence layer. This distinction accommodates regulatory",
  "interactions that may affect translation without producing a",
  "qualifying change in steady-state target mRNA abundance. RsaOI is",
  "detected and retained in the revised CLASH network. However, the",
  "previously reported RsaOI-atl and RsaOI-HPr interactions were not",
  "present in the analyzed raw CLASH interaction table; their",
  "non-recovery is therefore treated as a dataset/recovery limitation",
  "rather than evidence of biological absence."
)


writeLines(
  
  reviewer_text,
  
  file.path(
    outdir,
    "Figure2_FINAL_CLEAN_v2_Reviewer_Response.txt"
  )
)


# ============================================================
# 28. SESSION INFO
# ============================================================

writeLines(
  
  capture.output(
    sessionInfo()
  ),
  
  file.path(
    outdir,
    "sessionInfo_v2.txt"
  )
)


# ============================================================
# 29. FINAL CONSOLE AUDIT
# ============================================================

cat("\n\n")
cat("============================================\n")
cat("SCRIPT 20 FINAL v2 COMPLETED\n")
cat("============================================\n")


cat("\nFINAL COUNTS\n")
cat("--------------------------------------------\n")

cat(
  "Complete CLASH network: ",
  n_total,
  "\n",
  sep = ""
)

cat(
  "Higher/Intermediate CLASH support: ",
  n_hi_int,
  "\n",
  sep = ""
)

cat(
  "Interactions with significant target DEGs: ",
  n_deg_edges,
  "\n",
  sep = ""
)

cat(
  "Unique significant target genes: ",
  n_deg_targets,
  "\n",
  sep = ""
)


cat("\nTRANSCRIPTOMIC DISTRIBUTION\n")
cat("--------------------------------------------\n")

print(
  transcriptomic_distribution,
  n = Inf,
  width = Inf
)


cat("\nANALYTICAL INTERPRETATION\n")
cat("--------------------------------------------\n")

cat(
  paste0(
    
    "211 = complete CLASH-supported network.\n",
    
    "100 = Higher/Intermediate descriptive CLASH support.\n",
    
    "10 = CLASH interactions involving significant target DEGs.\n\n",
    
    "IMPORTANT:\n",
    
    "100 and 10 are independent summaries derived from 211.\n",
    
    "The 10 interactions are NOT a downstream filter of the 100.\n",
    
    "Target DEG status is NOT required for network inclusion.\n",
    
    "CLASH and transcriptomic response remain separate evidence layers.\n"
  )
)


cat("\nOUTPUT DIRECTORY\n")
cat("--------------------------------------------\n")

cat(
  outdir,
  "\n"
)


cat("\nMAIN FIGURE TO REVIEW\n")
cat("--------------------------------------------\n")

cat(
  file.path(
    figdir,
    "Figure2_FINAL_CLEAN_v2_600dpi.png"
  ),
  "\n"
)


cat("\n============================================\n")
cat("END SCRIPT 20 FINAL v2\n")
cat("============================================\n")