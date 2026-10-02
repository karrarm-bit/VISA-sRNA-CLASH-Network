# ============================================================
# SCRIPT 19
# FINAL FIGURE 2 — Evidence-layered RNase III-CLASH network
#
# Panels:
# A. Network architecture / evidence-layer logic
# B. Vancomycin-responsive target subset (10 interactions)
# C. Benchmark sRNA-associated CLASH interactions + nrdF candidates
# D. Target transcriptomic evidence across all 211 CLASH edges
#
# CENTRAL PRINCIPLE:
# RNase III-CLASH defines the interaction network.
# Target transcriptional response is an independent evidence layer.
# Target DEG status is NOT an inclusion criterion.
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
  "19_final_Figure2_evidence_layered_network"
)

figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")

dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)


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

summary_file <- file.path(
  input_dir,
  "Evidence_Layer_Summary.csv"
)


# ============================================================
# 04. INPUT CHECK
# ============================================================

input_check <- tibble(
  file = c(
    "Complete CLASH network",
    "Responsive target subset",
    "Benchmark/nrdF subset",
    "Evidence summary"
  ),
  path = c(
    all_file,
    responsive_file,
    benchmark_file,
    summary_file
  ),
  exists = c(
    file.exists(all_file),
    file.exists(responsive_file),
    file.exists(benchmark_file),
    file.exists(summary_file)
  )
)

cat("\n============================================\n")
cat("INPUT CHECK\n")
cat("============================================\n")

print(input_check, n = Inf, width = Inf)

if (any(!input_check$exists)) {
  stop("One or more required Script 14 outputs are missing.")
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

evidence_summary <- read_csv(
  summary_file,
  show_col_types = FALSE
)


# ============================================================
# 06. BASIC AUDIT
# ============================================================

cat("\n============================================\n")
cat("NETWORK AUDIT\n")
cat("============================================\n")

cat("Complete CLASH network:", nrow(net), "\n")
cat("Responsive-target subset:", nrow(responsive), "\n")
cat("Benchmark/nrdF subset:", nrow(benchmark), "\n")


if (nrow(net) != 211) {
  warning(
    paste(
      "Expected 211 complete CLASH interactions;",
      "observed", nrow(net)
    )
  )
}

if (nrow(responsive) != 10) {
  warning(
    paste(
      "Expected 10 responsive-target interactions;",
      "observed", nrow(responsive)
    )
  )
}


# ============================================================
# 07. STANDARDIZE LABELS
# ============================================================

clean_source_label <- function(x) {
  
  x <- as.character(x)
  
  case_when(
    
    str_detect(
      x,
      regex("RS04120", ignore_case = TRUE)
    ) ~ "RS04120",
    
    str_detect(
      x,
      regex("RS05585", ignore_case = TRUE)
    ) ~ "RS05585",
    
    str_detect(
      x,
      regex("RsaOI|Sau-6477", ignore_case = TRUE)
    ) ~ "RsaOI",
    
    str_detect(
      x,
      regex("SprA2/SprAs2", ignore_case = TRUE)
    ) ~ "SprA2/SprAs2",
    
    str_detect(
      x,
      regex("SprA2/RsaJ", ignore_case = TRUE)
    ) ~ "SprA2/RsaJ",
    
    TRUE ~ x
  )
}


clean_target_label <- function(x) {
  
  x <- as.character(x)
  
  case_when(
    
    str_detect(
      x,
      regex("^nrdF$", ignore_case = TRUE)
    ) ~ "nrdF",
    
    str_detect(
      x,
      regex("^guaB$", ignore_case = TRUE)
    ) ~ "guaB",
    
    str_detect(
      x,
      regex("^spa$", ignore_case = TRUE)
    ) ~ "spa",
    
    str_detect(
      x,
      regex("^qoxB$", ignore_case = TRUE)
    ) ~ "qoxB",
    
    str_detect(
      x,
      regex("^gltB$", ignore_case = TRUE)
    ) ~ "gltB",
    
    str_detect(
      x,
      regex("^dhaK$", ignore_case = TRUE)
    ) ~ "dhaK",
    
    str_detect(
      x,
      regex("^modA$", ignore_case = TRUE)
    ) ~ "modA",
    
    TRUE ~ x
  )
}


# ============================================================
# 08. PANEL A DATA
# ============================================================

higher_intermediate_n <- net %>%
  filter(
    clash_support_category %in%
      c(
        "Higher CLASH support",
        "Intermediate CLASH support",
        "Higher",
        "Intermediate"
      )
  ) %>%
  nrow()


# fallback to established result if wording differs
if (higher_intermediate_n == 0) {
  
  higher_intermediate_n <- net %>%
    filter(
      str_detect(
        clash_support_category,
        regex(
          "Higher|Intermediate",
          ignore_case = TRUE
        )
      )
    ) %>%
    nrow()
}


cat(
  "Higher/intermediate CLASH support:",
  higher_intermediate_n,
  "\n"
)


architecture <- tibble(
  
  stage = factor(
    c(
      "Complete CLASH network",
      "Higher / Intermediate\nCLASH support",
      "Significant target-DEG\nassociated subset"
    ),
    levels = c(
      "Complete CLASH network",
      "Higher / Intermediate\nCLASH support",
      "Significant target-DEG\nassociated subset"
    )
  ),
  
  n = c(
    nrow(net),
    higher_intermediate_n,
    nrow(responsive)
  ),
  
  x = c(
    1,
    2,
    3
  )
)


# ============================================================
# 09. PANEL A
# ============================================================

panel_A <- ggplot(
  architecture,
  aes(
    x = x,
    y = 1
  )
) +
  
  geom_segment(
    aes(
      x = 1.20,
      xend = 1.80,
      y = 1,
      yend = 1
    ),
    arrow = arrow(
      length = unit(
        0.16,
        "cm"
      ),
      type = "closed"
    ),
    linewidth = 0.8
  ) +
  
  geom_segment(
    aes(
      x = 2.20,
      xend = 2.80,
      y = 1,
      yend = 1
    ),
    arrow = arrow(
      length = unit(
        0.16,
        "cm"
      ),
      type = "closed"
    ),
    linewidth = 0.8
  ) +
  
  geom_label(
    aes(
      label = paste0(
        n,
        "\n",
        stage
      )
    ),
    size = 4,
    fontface = "bold",
    label.size = 0.5,
    label.padding = unit(
      0.35,
      "lines"
    )
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = 0.72,
    label =
      "descriptive CLASH-support stratification",
    size = 3.0
  ) +
  
  annotate(
    "text",
    x = 2.5,
    y = 0.72,
    label =
      "independent transcriptomic layer",
    size = 3.0
  ) +
  
  annotate(
    "text",
    x = 2,
    y = 0.37,
    label =
      paste(
        "Target differential expression was NOT required",
        "for inclusion in the CLASH interaction network."
      ),
    fontface = "bold",
    size = 3.5
  ) +
  
  coord_cartesian(
    xlim = c(
      0.55,
      3.45
    ),
    ylim = c(
      0.20,
      1.40
    ),
    clip = "off"
  ) +
  
  labs(
    title =
      "A  Evidence-layer architecture of the RNase III-CLASH network",
    
    subtitle =
      paste(
        "CLASH defines the interaction universe;",
        "target transcriptional response is evaluated separately"
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
        size = 9
      ),
    
    plot.margin =
      margin(
        10,
        20,
        10,
        20
      )
  )


# ============================================================
# 10. PREPARE PANEL B
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
      ),
    
    source_class =
      as.character(
        source_RNA_class
      ),
    
    target_response =
      as.character(
        target_vancomycin_response
      )
  )


# source positions
source_order_B <- responsive_plot %>%
  
  distinct(
    source_display
  ) %>%
  
  arrange(
    source_display
  ) %>%
  
  mutate(
    source_y =
      rev(
        seq_len(n())
      )
  )


# target positions
target_order_B <- responsive_plot %>%
  
  distinct(
    target_display
  ) %>%
  
  arrange(
    target_display
  ) %>%
  
  mutate(
    target_y =
      seq(
        from =
          max(
            1,
            nrow(source_order_B)
          ),
        
        to =
          1,
        
        length.out =
          n()
      )
  )


responsive_plot <- responsive_plot %>%
  
  left_join(
    source_order_B,
    by = "source_display"
  ) %>%
  
  left_join(
    target_order_B,
    by = "target_display"
  )


# ============================================================
# 11. PANEL B
# ============================================================

panel_B <- ggplot() +
  
  geom_curve(
    data = responsive_plot,
    aes(
      x = 1,
      y = source_y,
      xend = 3,
      yend = target_y,
      linewidth = hybrid_count,
      linetype = clash_support_category
    ),
    curvature = 0.08,
    alpha = 0.75
  ) +
  
  geom_point(
    data = source_order_B,
    aes(
      x = 1,
      y = source_y
    ),
    shape = 22,
    size = 3.5,
    fill = "white",
    stroke = 0.8
  ) +
  
  geom_point(
    data = target_order_B,
    aes(
      x = 3,
      y = target_y
    ),
    shape = 21,
    size = 3.5,
    fill = "white",
    stroke = 0.8
  ) +
  
  geom_text(
    data = source_order_B,
    aes(
      x = 0.92,
      y = source_y,
      label = source_display
    ),
    hjust = 1,
    size = 3.1
  ) +
  
  geom_text(
    data = target_order_B,
    aes(
      x = 3.08,
      y = target_y,
      label = target_display
    ),
    hjust = 0,
    size = 3.1,
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = 1,
    y =
      max(source_order_B$source_y) +
      0.8,
    label =
      "Regulatory RNAs",
    fontface = "bold",
    size = 3.5
  ) +
  
  annotate(
    "text",
    x = 3,
    y =
      max(source_order_B$source_y) +
      0.8,
    label =
      "Responsive mRNA targets",
    fontface = "bold",
    size = 3.5
  ) +
  
  scale_linewidth_continuous(
    range = c(
      0.4,
      2.0
    ),
    name =
      "Hybrid count"
  ) +
  
  labs(
    title =
      "B  Vancomycin-responsive target subset",
    
    subtitle =
      paste0(
        nrow(responsive_plot),
        " CLASH-supported interactions involving ",
        n_distinct(
          responsive_plot$target_display
        ),
        " significantly responsive mRNA targets"
      ),
    
    linetype =
      "CLASH support"
  ) +
  
  coord_cartesian(
    xlim = c(
      0.25,
      3.75
    ),
    clip = "off"
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
        size = 9
      ),
    
    legend.position =
      "bottom",
    
    plot.margin =
      margin(
        10,
        45,
        10,
        45
      )
  )


# ============================================================
# 12. PREPARE PANEL C
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


# make duplicate target labels visually distinct only when needed
benchmark_plot <- benchmark_plot %>%
  
  group_by(
    target_display
  ) %>%
  
  mutate(
    target_display_plot =
      if_else(
        n() > 1 &
          target_display ==
          "WP_000482650.1",
        paste0(
          target_display,
          " [",
          row_number(),
          "]"
        ),
        target_display
      )
  ) %>%
  
  ungroup()


source_order_C <- benchmark_plot %>%
  
  mutate(
    source_display_plot =
      if_else(
        duplicated(source_display) |
          duplicated(
            source_display,
            fromLast = TRUE
          ),
        paste0(
          source_display,
          " [",
          row_number(),
          "]"
        ),
        source_display
      )
  )


# Use one row per interaction vertically
benchmark_plot <- benchmark_plot %>%
  
  mutate(
    interaction_y =
      rev(
        seq_len(n())
      )
  )


# ============================================================
# 13. PANEL C
# ============================================================

panel_C <- ggplot() +
  
  geom_curve(
    data = benchmark_plot,
    aes(
      x = 1,
      y = interaction_y,
      xend = 3,
      yend = interaction_y,
      linewidth = hybrid_count,
      linetype = clash_support_category
    ),
    curvature = 0.08,
    alpha = 0.80
  ) +
  
  geom_point(
    data = benchmark_plot,
    aes(
      x = 1,
      y = interaction_y
    ),
    shape = 22,
    size = 3.7,
    fill = "white",
    stroke = 0.9
  ) +
  
  geom_point(
    data = benchmark_plot,
    aes(
      x = 3,
      y = interaction_y
    ),
    shape = 21,
    size = 3.7,
    fill = "white",
    stroke = 0.9
  ) +
  
  geom_text(
    data = benchmark_plot,
    aes(
      x = 0.92,
      y = interaction_y,
      label = source_display
    ),
    hjust = 1,
    size = 3.2,
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
      x = 3.08,
      y = interaction_y,
      label = target_display
    ),
    hjust = 0,
    size = 3.2,
    fontface =
      ifelse(
        benchmark_plot$target_display == "nrdF",
        "bold",
        "plain"
      )
  ) +
  
  annotate(
    "text",
    x = 1,
    y =
      max(
        benchmark_plot$interaction_y
      ) +
      0.8,
    label =
      "Regulatory RNAs",
    fontface = "bold",
    size = 3.5
  ) +
  
  annotate(
    "text",
    x = 3,
    y =
      max(
        benchmark_plot$interaction_y
      ) +
      0.8,
    label =
      "CLASH-associated mRNA targets",
    fontface = "bold",
    size = 3.5
  ) +
  
  scale_linewidth_continuous(
    range = c(
      0.45,
      2.1
    ),
    name =
      "Hybrid count"
  ) +
  
  labs(
    title =
      "C  Benchmark sRNA-associated CLASH interactions and nrdF candidates",
    
    subtitle =
      paste(
        "RsaOI is retained in the complete CLASH network despite its",
        "recovered target not meeting the target-DEG criterion"
      ),
    
    linetype =
      "CLASH support"
  ) +
  
  coord_cartesian(
    xlim = c(
      0.25,
      3.75
    ),
    clip = "off"
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
        size = 8.7
      ),
    
    legend.position =
      "bottom",
    
    plot.margin =
      margin(
        10,
        55,
        10,
        55
      )
  )


# ============================================================
# 14. PANEL D DATA
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


# force desired order
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
      n /
      sum(n) *
      100,
    
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


cat("\n============================================\n")
cat("TARGET TRANSCRIPTOMIC DISTRIBUTION\n")
cat("============================================\n")

print(
  transcriptomic_distribution,
  n = Inf,
  width = Inf
)


# ============================================================
# 15. PANEL D
# ============================================================

panel_D <- ggplot(
  transcriptomic_distribution,
  aes(
    x = transcript_category,
    y = n
  )
) +
  
  geom_col(
    width = 0.65,
    fill = "grey80",
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
    vjust = -0.35,
    fontface = "bold",
    size = 3.6
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.14
      )
    )
  ) +
  
  labs(
    title =
      "D  Target transcriptomic evidence across the complete CLASH network",
    
    subtitle =
      paste0(
        "All ",
        nrow(net),
        " CLASH-supported interactions are retained regardless of target-DEG status"
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
        size = 8.5
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 11
      ),
    
    plot.subtitle =
      element_text(
        size = 9
      )
  )


# ============================================================
# 16. CENTRAL INTERPRETATION STRIP
# ============================================================

interpretation <- ggplot() +
  
  annotate(
    "text",
    x = 0,
    y = 0.5,
    hjust = 0,
    vjust = 0.5,
    fontface = "bold",
    size = 3.6,
    label =
      paste(
        "Interpretation:",
        "RNase III-CLASH defines the interaction network;",
        "target transcriptional response is evaluated as an independent evidence layer."
      )
  ) +
  
  xlim(
    0,
    1
  ) +
  
  ylim(
    0,
    1
  ) +
  
  theme_void()


# ============================================================
# 17. COMBINE FIGURE
# ============================================================

figure2 <- (
  
  panel_A /
    
    (panel_B | panel_C) /
    
    panel_D /
    
    interpretation
  
) +
  
  plot_layout(
    heights = c(
      0.75,
      1.55,
      1.05,
      0.18
    )
  ) +
  
  plot_annotation(
    
    title =
      "Evidence-layered reconstruction of the vancomycin-associated RNase III-CLASH network",
    
    subtitle =
      paste(
        "Physical RNA-RNA interaction evidence and target transcriptional response",
        "are retained as analytically distinct evidence layers."
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
# 18. DISPLAY
# ============================================================

print(figure2)


# ============================================================
# 19. EXPORT PDF
# ============================================================

ggsave(
  filename = file.path(
    figdir,
    "Figure2_Evidence_Layered_CLASH_Network_FINAL.pdf"
  ),
  plot = figure2,
  width = 13,
  height = 13,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)


# ============================================================
# 20. EXPORT PNG 600 DPI
# ============================================================

ggsave(
  filename = file.path(
    figdir,
    "Figure2_Evidence_Layered_CLASH_Network_FINAL_600dpi.png"
  ),
  plot = figure2,
  width = 13,
  height = 13,
  units = "in",
  dpi = 600,
  limitsize = FALSE
)


# ============================================================
# 21. EXPORT TIFF 600 DPI
# ============================================================

ggsave(
  filename = file.path(
    figdir,
    "Figure2_Evidence_Layered_CLASH_Network_FINAL_600dpi.tiff"
  ),
  plot = figure2,
  width = 13,
  height = 13,
  units = "in",
  dpi = 600,
  compression = "lzw",
  limitsize = FALSE
)


# ============================================================
# 22. SAVE FIGURE SOURCE TABLES
# ============================================================

write_csv(
  architecture,
  file.path(
    tabdir,
    "Figure2A_Network_Architecture.csv"
  )
)

write_csv(
  responsive_plot,
  file.path(
    tabdir,
    "Figure2B_Responsive_Target_Subset.csv"
  )
)

write_csv(
  benchmark_plot,
  file.path(
    tabdir,
    "Figure2C_Benchmark_nrdF_Context.csv"
  )
)

write_csv(
  transcriptomic_distribution,
  file.path(
    tabdir,
    "Figure2D_Target_Transcriptomic_Distribution.csv"
  )
)


# ============================================================
# 23. FIGURE LEGEND
# ============================================================

figure_legend <- paste(
  
  "Figure 2. Evidence-layered reconstruction of the vancomycin-associated",
  "RNase III-CLASH regulatory network.",
  
  "(A) RNase III-CLASH defined a complete network of 211 collapsed",
  "sRNA-mRNA interactions. One hundred interactions showed higher or",
  "intermediate descriptive CLASH support, whereas only 10 interactions",
  "involved mRNA targets meeting the prespecified differential-expression",
  "criterion. Target differential expression was not required for",
  "inclusion in the CLASH network.",
  
  "(B) The vancomycin-responsive target subset comprised 10 CLASH-supported",
  "interactions involving seven significantly responsive mRNA targets.",
  
  "(C) Selected benchmark sRNA-associated CLASH interactions and the two",
  "putative 3′UTR-associated nrdF candidate interactions. RsaOI was",
  "recovered in the CLASH dataset and retained in the complete network",
  "despite its recovered target not meeting the target differential-expression",
  "criterion. This panel does not imply recovery of all previously reported",
  "RsaOI targets.",
  
  "(D) Distribution of target transcriptomic evidence across all 211",
  "CLASH-supported interactions: 10 interactions involved significant",
  "target DEGs, 71 involved targets that were measured but remained below",
  "the prespecified DEG threshold, and 130 lacked mapped target",
  "transcriptomic evidence in the integrated analysis.",
  
  "Thus, physical interaction evidence and target transcriptional response",
  "were treated as independent evidence layers rather than requiring",
  "mRNA differential expression for network inclusion."
)


writeLines(
  figure_legend,
  file.path(
    outdir,
    "Figure2_FINAL_Legend.txt"
  )
)


# ============================================================
# 24. MANUSCRIPT-SAFE RESULTS TEXT
# ============================================================

results_text <- paste(
  
  "To avoid restricting the regulatory network to interactions whose",
  "targets changed at the steady-state mRNA level, RNase III-CLASH and",
  "target transcriptomic response were analyzed as independent evidence",
  "layers. The complete collapsed CLASH network contained 211 sRNA-mRNA",
  "interactions, of which 100 showed higher or intermediate descriptive",
  "CLASH support. Only 10 interactions involved targets meeting the",
  "prespecified vancomycin-responsive differential-expression criterion,",
  "corresponding to seven responsive mRNAs. Across the complete network,",
  "71 additional interactions involved targets that were measured but",
  "remained below the differential-expression threshold, whereas 130",
  "lacked mapped target transcriptomic evidence in the integrated dataset.",
  "Importantly, known sRNA-associated CLASH signals such as RsaOI were",
  "retained independently of target differential-expression status.",
  "Accordingly, absence from the responsive-target subset was not",
  "interpreted as absence of an RNA-RNA interaction."
)


writeLines(
  results_text,
  file.path(
    outdir,
    "Figure2_FINAL_Results_Paragraph.txt"
  )
)


# ============================================================
# 25. REVIEWER-SAFE NOTE
# ============================================================

reviewer_note <- paste(
  
  "The revised analysis no longer uses target mRNA differential expression",
  "as a prerequisite for retention of CLASH-supported interactions.",
  "RNase III-CLASH defines the interaction network, while target",
  "transcriptional response is evaluated independently. This distinction",
  "is important because sRNA-mediated regulation may occur at the",
  "translational level without a corresponding change in steady-state",
  "target mRNA abundance. RsaOI was detected in the raw CLASH dataset",
  "and is retained in the complete revised network. However, the specific",
  "previously reported RsaOI-atl and RsaOI-HPr pairs were not present in",
  "the analyzed raw CLASH interaction table; their non-recovery is",
  "therefore treated as a dataset/recovery limitation rather than evidence",
  "of biological absence."
)


writeLines(
  reviewer_note,
  file.path(
    outdir,
    "Figure2_Reviewer_Response_Note.txt"
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
# 27. FINAL AUDIT
# ============================================================

cat("\n\n============================================\n")
cat("SCRIPT 19 COMPLETED\n")
cat("============================================\n")

cat("\nNETWORK COUNTS\n")
cat("--------------------------------------------\n")

cat(
  "Complete CLASH interactions:",
  nrow(net),
  "\n"
)

cat(
  "Higher/intermediate CLASH interactions:",
  higher_intermediate_n,
  "\n"
)

cat(
  "Significant target-DEG interactions:",
  nrow(responsive),
  "\n"
)

cat(
  "Unique significant target genes:",
  n_distinct(
    responsive_plot$target_display
  ),
  "\n"
)


cat("\nTARGET TRANSCRIPTOMIC DISTRIBUTION\n")
cat("--------------------------------------------\n")

print(
  transcriptomic_distribution,
  n = Inf,
  width = Inf
)


cat("\nBENCHMARK/nrdF CONTEXT\n")
cat("--------------------------------------------\n")

print(
  benchmark_plot %>%
    select(
      source_display,
      target_display,
      total_hybrid_count,
      maximum_number_of_experiments,
      best_connection_score,
      clash_evidence_score,
      clash_support_category,
      target_transcriptomic_layer
    ),
  n = Inf,
  width = Inf
)


cat("\nINTERPRETATION\n")
cat("--------------------------------------------\n")

cat(
  paste0(
    "CLASH defines the network.\n",
    "Target DEG status is an independent evidence layer.\n",
    "RsaOI remains in the complete network.\n",
    "Absence from the responsive-target subset does not imply absence of interaction.\n",
    "No mechanistic or causal inference is made from target transcriptomic response alone.\n"
  )
)


cat("\nOUTPUT DIRECTORY\n")
cat("--------------------------------------------\n")

cat(
  outdir,
  "\n"
)

cat("\nMain files:\n")

cat(
  "1. Figure2_Evidence_Layered_CLASH_Network_FINAL.pdf\n"
)

cat(
  "2. Figure2_Evidence_Layered_CLASH_Network_FINAL_600dpi.png\n"
)

cat(
  "3. Figure2_Evidence_Layered_CLASH_Network_FINAL_600dpi.tiff\n"
)

cat(
  "4. Figure2_FINAL_Legend.txt\n"
)

cat(
  "5. Figure2_FINAL_Results_Paragraph.txt\n"
)

cat(
  "6. Figure2_Reviewer_Response_Note.txt\n"
)

cat("\n============================================\n")
cat("END SCRIPT 19\n")
cat("============================================\n")