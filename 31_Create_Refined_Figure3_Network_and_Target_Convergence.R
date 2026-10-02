# ============================================================
# 31_Create_Refined_Figure3_Fixed.R
#
# Purpose:
# Create a refined Figure 3 from the complete prioritised
# regulatory RNA-mRNA edge table.
#
# Panel A:
# Complete prioritised network of 10 interactions.
#
# Panel B:
# Target-convergence summary showing the number of distinct
# regulatory RNAs linked to each target gene.
#
# Input:
# D:/Bac-sRNA/15_network_figures/tables/
# Figure2B_complete_network_edges.csv
#
# This version does not require a separate node table.
# ============================================================


# ------------------------------------------------------------
# 1. Install and load required packages
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "readr",
  "stringr",
  "ggplot2",
  "patchwork",
  "svglite",
  "scales"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(patchwork)
library(svglite)
library(scales)


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

input_folder <- file.path(
  project_folder,
  "15_network_figures",
  "tables"
)

input_file <- file.path(
  input_folder,
  "Figure2B_complete_network_edges.csv"
)

output_folder <- file.path(
  project_folder,
  "15_network_figures",
  "Figure3_refined_fixed"
)

if (!dir.exists(output_folder)) {
  dir.create(
    output_folder,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 3. Detect input file automatically
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  
  possible_files <- list.files(
    path = input_folder,
    pattern = "Figure2B.*complete.*network.*edges.*\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(possible_files) == 0) {
    
    possible_files <- list.files(
      path = input_folder,
      pattern = "complete.*network.*edges.*\\.csv$",
      full.names = TRUE,
      ignore.case = TRUE
    )
  }
  
  if (length(possible_files) == 0) {
    stop(
      paste0(
        "Complete network edge file was not found in:\n",
        input_folder
      )
    )
  }
  
  input_file <- possible_files[1]
  
  message(
    "Detected edge file:\n",
    input_file
  )
}


# ------------------------------------------------------------
# 4. Output files
# ------------------------------------------------------------

output_png <- file.path(
  output_folder,
  "Figure3_refined_network_and_target_convergence.png"
)

output_tiff <- file.path(
  output_folder,
  "Figure3_refined_network_and_target_convergence.tiff"
)

output_pdf <- file.path(
  output_folder,
  "Figure3_refined_network_and_target_convergence.pdf"
)

output_svg <- file.path(
  output_folder,
  "Figure3_refined_network_and_target_convergence.svg"
)

output_panel_a <- file.path(
  output_folder,
  "Figure3A_complete_prioritised_network.png"
)

output_panel_b <- file.path(
  output_folder,
  "Figure3B_target_convergence_summary.png"
)

output_edges_csv <- file.path(
  output_folder,
  "Figure3_refined_edge_data.csv"
)

output_nodes_csv <- file.path(
  output_folder,
  "Figure3_refined_node_data.csv"
)

output_convergence_csv <- file.path(
  output_folder,
  "Figure3_target_convergence_summary.csv"
)


# ------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------

clean_column_names <- function(data) {
  
  names(data) <- names(data) |>
    str_trim() |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("^_|_$", "")
  
  data
}


clean_text <- function(x) {
  
  x <- as.character(x)
  x <- str_trim(x)
  
  x[
    is.na(x) |
      x == "" |
      str_to_upper(x) %in% c(
        "NA",
        "N/A",
        "NULL",
        "NONE",
        "-"
      )
  ] <- NA_character_
  
  x
}


safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}


normalise_id <- function(x) {
  
  x |>
    clean_text() |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("^_|_$", "")
}


# ------------------------------------------------------------
# 6. Read edge file
# ------------------------------------------------------------

edges_raw <- read_csv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
) |>
  clean_column_names()


# ------------------------------------------------------------
# 7. Confirm expected columns
# ------------------------------------------------------------

required_columns <- c(
  "source_label",
  "target_label",
  "source_class",
  "target_response",
  "target_log2_fold_change",
  "hybrid_count",
  "confidence_score",
  "confidence_class"
)

missing_columns <- setdiff(
  required_columns,
  names(edges_raw)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "Required columns are missing:\n",
      paste(
        missing_columns,
        collapse = ", "
      ),
      "\n\nAvailable columns:\n",
      paste(
        names(edges_raw),
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 8. Add optional columns when absent
# ------------------------------------------------------------

if (!"number_of_experiments" %in% names(edges_raw)) {
  edges_raw$number_of_experiments <- 1
}

if (!"target_adjusted_p_value" %in% names(edges_raw)) {
  edges_raw$target_adjusted_p_value <- NA_real_
}

if (!"included_in_high_confidence" %in% names(edges_raw)) {
  edges_raw$included_in_high_confidence <- FALSE
}


# ------------------------------------------------------------
# 9. Clean edge data
# ------------------------------------------------------------

edges <- edges_raw |>
  transmute(
    
    edge_id = if (
      "edge_id" %in% names(edges_raw)
    ) {
      clean_text(edge_id)
    } else {
      paste0(
        "EDGE_",
        sprintf(
          "%03d",
          row_number()
        )
      )
    },
    
    source =
      clean_text(source_label),
    
    target =
      clean_text(target_label),
    
    source_class =
      clean_text(source_class),
    
    target_response =
      clean_text(target_response),
    
    target_log2_fold_change =
      safe_numeric(target_log2_fold_change),
    
    target_adjusted_p_value =
      safe_numeric(target_adjusted_p_value),
    
    hybrid_count =
      safe_numeric(hybrid_count),
    
    number_of_experiments =
      safe_numeric(number_of_experiments),
    
    confidence_score =
      safe_numeric(confidence_score),
    
    confidence_class =
      clean_text(confidence_class),
    
    included_in_high_confidence =
      as.logical(included_in_high_confidence)
  ) |>
  
  filter(
    !is.na(source),
    !is.na(target)
  ) |>
  
  mutate(
    
    source_key =
      normalise_id(source),
    
    target_key =
      normalise_id(target),
    
    confidence_class_plot = case_when(
      
      str_detect(
        str_to_lower(confidence_class),
        "high"
      ) ~
        "High confidence",
      
      str_detect(
        str_to_lower(confidence_class),
        "moderate"
      ) ~
        "Moderate confidence",
      
      str_detect(
        str_to_lower(confidence_class),
        "supported"
      ) ~
        "Supported candidate",
      
      TRUE ~
        "Supported candidate"
    ),
    
    target_direction = case_when(
      
      target_log2_fold_change < 0 ~
        "Downregulated target",
      
      target_log2_fold_change > 0 ~
        "Upregulated target",
      
      str_detect(
        str_to_lower(target_response),
        "repress|down"
      ) ~
        "Downregulated target",
      
      str_detect(
        str_to_lower(target_response),
        "induc|up"
      ) ~
        "Upregulated target",
      
      TRUE ~
        "Target response unavailable"
    ),
    
    rna_class_plot = case_when(
      
      str_detect(
        str_to_lower(source_class),
        "intergenic"
      ) ~
        "Intergenic candidate RNA",
      
      str_detect(
        str_to_lower(source_class),
        "named"
      ) ~
        "Named sRNA",
      
      str_detect(
        str_to_lower(source_class),
        "numbered"
      ) ~
        "Numbered sRNA candidate",
      
      str_detect(
        str_to_lower(source_class),
        "utr"
      ) ~
        "3′UTR-derived RNA",
      
      TRUE ~
        source_class
    )
  )


# ------------------------------------------------------------
# 10. Calculate scaled edge widths
# ------------------------------------------------------------

edges <- edges |>
  mutate(
    
    hybrid_count_for_scaling = case_when(
      is.na(hybrid_count) ~ 1,
      hybrid_count <= 0 ~ 1,
      TRUE ~ hybrid_count
    ),
    
    edge_width = scales::rescale(
      log1p(hybrid_count_for_scaling),
      to = c(
        0.55,
        2.40
      )
    )
  )


# ------------------------------------------------------------
# 11. Build regulatory RNA node table directly from edges
# ------------------------------------------------------------

regulatory_nodes <- edges |>
  group_by(
    source,
    source_key,
    source_class,
    rna_class_plot
  ) |>
  summarise(
    
    source_degree =
      n_distinct(target_key),
    
    total_hybrid_support =
      sum(
        hybrid_count,
        na.rm = TRUE
      ),
    
    maximum_confidence_score = if (
      all(is.na(confidence_score))
    ) {
      NA_real_
    } else {
      max(
        confidence_score,
        na.rm = TRUE
      )
    },
    
    .groups = "drop"
  ) |>
  
  arrange(
    desc(maximum_confidence_score),
    desc(total_hybrid_support),
    source
  ) |>
  
  mutate(
    node_id =
      source,
    
    node_key =
      source_key,
    
    node_label =
      source,
    
    node_role =
      "Regulatory RNA",
    
    x =
      0,
    
    y =
      rev(
        seq_len(n())
      )
  )


# ------------------------------------------------------------
# 12. Calculate target convergence
# ------------------------------------------------------------

target_convergence <- edges |>
  group_by(
    target,
    target_key
  ) |>
  
  summarise(
    
    number_of_regulatory_rnas =
      n_distinct(source_key),
    
    linked_regulatory_rnas =
      paste(
        unique(source),
        collapse = "; "
      ),
    
    target_log2_fold_change =
      first(
        target_log2_fold_change
      ),
    
    target_adjusted_p_value =
      first(
        target_adjusted_p_value
      ),
    
    target_direction =
      first(
        target_direction
      ),
    
    total_hybrid_support =
      sum(
        hybrid_count,
        na.rm = TRUE
      ),
    
    maximum_confidence_score = if (
      all(is.na(confidence_score))
    ) {
      NA_real_
    } else {
      max(
        confidence_score,
        na.rm = TRUE
      )
    },
    
    strongest_confidence_class = case_when(
      
      any(
        confidence_class_plot ==
          "High confidence"
      ) ~
        "High confidence",
      
      any(
        confidence_class_plot ==
          "Moderate confidence"
      ) ~
        "Moderate confidence",
      
      TRUE ~
        "Supported candidate"
    ),
    
    .groups = "drop"
  ) |>
  
  arrange(
    desc(number_of_regulatory_rnas),
    desc(total_hybrid_support),
    target
  )


# ------------------------------------------------------------
# 13. Build target node table
# ------------------------------------------------------------

target_nodes <- target_convergence |>
  mutate(
    
    node_id =
      target,
    
    node_key =
      target_key,
    
    node_label =
      target,
    
    node_role =
      "Target gene",
    
    expression_direction =
      target_direction,
    
    x =
      1,
    
    y = seq(
      from = nrow(regulatory_nodes),
      to = 1,
      length.out = n()
    )
  )


# ------------------------------------------------------------
# 14. Combine node tables
# ------------------------------------------------------------

node_positions <- bind_rows(
  
  regulatory_nodes |>
    transmute(
      node_id,
      node_key,
      node_label,
      node_role,
      rna_class_plot,
      expression_direction =
        "Regulatory RNA",
      number_of_regulatory_rnas =
        source_degree,
      total_hybrid_support,
      x,
      y
    ),
  
  target_nodes |>
    transmute(
      node_id,
      node_key,
      node_label,
      node_role,
      rna_class_plot =
        "Target gene",
      expression_direction,
      number_of_regulatory_rnas,
      total_hybrid_support,
      x,
      y
    )
)


# ------------------------------------------------------------
# 15. Attach coordinates to edges
#
# This is where the previous target_key error occurred.
# Both source_key and target_key now exist explicitly.
# ------------------------------------------------------------

edge_plot_data <- edges |>
  
  left_join(
    regulatory_nodes |>
      select(
        source_key,
        source_x = x,
        source_y = y
      ),
    by = "source_key"
  ) |>
  
  left_join(
    target_nodes |>
      select(
        target_key,
        target_x = x,
        target_y = y
      ),
    by = "target_key"
  )


# ------------------------------------------------------------
# 16. Verify coordinates
# ------------------------------------------------------------

if (
  any(
    is.na(edge_plot_data$source_x)
  ) ||
  any(
    is.na(edge_plot_data$target_x)
  )
) {
  
  unmatched_edges <- edge_plot_data |>
    filter(
      is.na(source_x) |
        is.na(target_x)
    )
  
  write_csv(
    unmatched_edges,
    file.path(
      output_folder,
      "Unmatched_network_edges_audit.csv"
    ),
    na = ""
  )
  
  stop(
    paste0(
      "Some edge coordinates could not be matched.\n",
      "See Unmatched_network_edges_audit.csv in:\n",
      output_folder
    )
  )
}


# ------------------------------------------------------------
# 17. Save processed data
# ------------------------------------------------------------

write_csv(
  edge_plot_data,
  output_edges_csv,
  na = ""
)

write_csv(
  node_positions,
  output_nodes_csv,
  na = ""
)

write_csv(
  target_convergence,
  output_convergence_csv,
  na = ""
)


# ============================================================
# PANEL A
# Complete prioritised interaction network
# ============================================================


# ------------------------------------------------------------
# 18. Create Panel A
# ------------------------------------------------------------

panel_a <- ggplot() +
  
  geom_curve(
    data = edge_plot_data,
    aes(
      x = source_x + 0.035,
      y = source_y,
      xend = target_x - 0.035,
      yend = target_y,
      linewidth = edge_width,
      linetype = confidence_class_plot
    ),
    curvature = 0.08,
    colour = "#666666",
    alpha = 0.85,
    arrow = arrow(
      length = grid::unit(
        0.085,
        "inches"
      ),
      type = "closed"
    )
  ) +
  
  geom_point(
    data = regulatory_nodes,
    aes(
      x = x,
      y = y,
      shape = rna_class_plot
    ),
    size = 5.2,
    fill = "#E9C46A",
    colour = "#303030",
    stroke = 0.75
  ) +
  
  geom_point(
    data = target_nodes,
    aes(
      x = x,
      y = y,
      fill = expression_direction,
      size = number_of_regulatory_rnas
    ),
    shape = 21,
    colour = "#202020",
    stroke = 0.75
  ) +
  
  geom_text(
    data = regulatory_nodes,
    aes(
      x = x - 0.045,
      y = y,
      label = source
    ),
    hjust = 1,
    size = 3.25,
    fontface = "bold",
    family = "Arial"
  ) +
  
  geom_text(
    data = target_nodes,
    aes(
      x = x + 0.045,
      y = y,
      label = target
    ),
    hjust = 0,
    size = 3.4,
    fontface = "italic",
    family = "Arial"
  ) +
  
  annotate(
    "text",
    x = 0,
    y = max(node_positions$y) + 1.1,
    label = "Regulatory RNAs",
    fontface = "bold",
    size = 4.4
  ) +
  
  annotate(
    "text",
    x = 1,
    y = max(node_positions$y) + 1.1,
    label = "Vancomycin-responsive targets",
    fontface = "bold",
    size = 4.4
  ) +
  
  scale_shape_manual(
    name = "Regulatory RNA class",
    values = c(
      "Intergenic candidate RNA" = 24,
      "Named sRNA" = 22,
      "Numbered sRNA candidate" = 23,
      "3′UTR-derived RNA" = 21
    )
  ) +
  
  scale_fill_manual(
    name = "Target response",
    values = c(
      "Downregulated target" = "#3B83BD",
      "Upregulated target" = "#E67E22",
      "Target response unavailable" = "#BDBDBD"
    )
  ) +
  
  scale_linetype_manual(
    name = "Interaction confidence",
    values = c(
      "High confidence" = "solid",
      "Moderate confidence" = "longdash",
      "Supported candidate" = "dotted"
    )
  ) +
  
  scale_linewidth_identity() +
  
  scale_size_continuous(
    name = "No. of targeting RNAs",
    range = c(
      4.8,
      7.4
    ),
    breaks = sort(
      unique(
        target_nodes$
          number_of_regulatory_rnas
      )
    )
  ) +
  
  coord_cartesian(
    xlim = c(
      -0.62,
      1.60
    ),
    ylim = c(
      0.3,
      max(node_positions$y) + 1.5
    ),
    clip = "off"
  ) +
  
  theme_void(
    base_size = 11
  ) +
  
  theme(
    text = element_text(
      family = "Arial"
    ),
    
    legend.position = "bottom",
    
    legend.box = "vertical",
    
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    plot.margin = margin(
      24,
      38,
      18,
      38
    )
  )


# ============================================================
# PANEL B
# Target-convergence summary
# ============================================================


# ------------------------------------------------------------
# 19. Prepare Panel B data
# ------------------------------------------------------------

convergence_plot_data <- target_convergence |>
  mutate(
    
    target_factor = factor(
      target,
      levels = rev(target)
    ),
    
    expression_direction = case_when(
      
      target_log2_fold_change < 0 ~
        "Downregulated",
      
      target_log2_fold_change > 0 ~
        "Upregulated",
      
      TRUE ~
        "Unavailable"
    ),
    
    expression_label = case_when(
      
      expression_direction ==
        "Downregulated" ~
        paste0(
          "\u2193  log2FC = ",
          sprintf(
            "%.2f",
            target_log2_fold_change
          )
        ),
      
      expression_direction ==
        "Upregulated" ~
        paste0(
          "\u2191  log2FC = ",
          sprintf(
            "%.2f",
            target_log2_fold_change
          )
        ),
      
      TRUE ~
        "Expression unavailable"
    )
  )


# ------------------------------------------------------------
# 20. Create Panel B
# ------------------------------------------------------------

maximum_target_count <- max(
  convergence_plot_data$
    number_of_regulatory_rnas,
  na.rm = TRUE
)


panel_b <- ggplot(
  convergence_plot_data,
  aes(
    x = number_of_regulatory_rnas,
    y = target_factor
  )
) +
  
  geom_segment(
    aes(
      x = 0,
      xend = number_of_regulatory_rnas,
      yend = target_factor
    ),
    linewidth = 1.1,
    colour = "#C9C9C9"
  ) +
  
  geom_point(
    aes(
      fill = expression_direction,
      size = total_hybrid_support
    ),
    shape = 21,
    colour = "#202020",
    stroke = 0.75
  ) +
  
  geom_text(
    aes(
      x = number_of_regulatory_rnas + 0.10,
      label = expression_label
    ),
    hjust = 0,
    size = 3.45,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    name = "Target response",
    values = c(
      "Downregulated" = "#3B83BD",
      "Upregulated" = "#E67E22",
      "Unavailable" = "#BDBDBD"
    )
  ) +
  
  scale_size_continuous(
    name = "Cumulative CLASH hybrid support",
    range = c(
      4.2,
      9
    )
  ) +
  
  scale_x_continuous(
    breaks = seq(
      0,
      maximum_target_count + 1,
      by = 1
    ),
    
    limits = c(
      0,
      maximum_target_count + 1.25
    ),
    
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  labs(
    x = "Number of distinct regulatory RNAs",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    text = element_text(
      family = "Arial"
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 12.5,
      margin = margin(
        t = 10
      )
    ),
    
    axis.text.x = element_text(
      size = 10.5,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      colour = "black",
      face = "italic"
    ),
    
    axis.line.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    panel.grid = element_blank(),
    
    legend.position = "bottom",
    
    legend.box = "vertical",
    
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    plot.margin = margin(
      24,
      40,
      18,
      15
    )
  )


# ============================================================
# COMBINE FIGURE
# ============================================================


# ------------------------------------------------------------
# 21. Combine panels
# ------------------------------------------------------------

figure3_refined <- (
  panel_a |
    panel_b
) +
  
  plot_layout(
    widths = c(
      1.20,
      0.80
    )
  ) +
  
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        family = "Arial",
        face = "bold",
        size = 18
      )
    )
  )


# ------------------------------------------------------------
# 22. Display figures
# ------------------------------------------------------------

print(panel_a)
print(panel_b)
print(figure3_refined)


# ------------------------------------------------------------
# 23. Save outputs
# ------------------------------------------------------------

ggsave(
  filename = output_panel_a,
  plot = panel_a,
  width = 9,
  height = 8.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = output_panel_b,
  plot = panel_b,
  width = 7,
  height = 8.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = output_png,
  plot = figure3_refined,
  width = 16,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = output_tiff,
  plot = figure3_refined,
  width = 16,
  height = 9,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = output_pdf,
  plot = figure3_refined,
  width = 16,
  height = 9,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = output_svg,
  plot = figure3_refined,
  width = 16,
  height = 9,
  units = "in",
  device = svglite,
  bg = "white"
)


# ------------------------------------------------------------
# 24. Final report
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "Refined Figure 3 created successfully\n",
  "============================================\n\n",
  
  "Input file:\n",
  normalizePath(
    input_file,
    winslash = "/"
  ),
  "\n\n",
  
  "Interactions displayed: ",
  nrow(edge_plot_data),
  "\n",
  
  "Regulatory RNAs displayed: ",
  nrow(regulatory_nodes),
  "\n",
  
  "Unique targets displayed: ",
  nrow(target_nodes),
  "\n",
  
  "Targets with multiple regulatory RNAs: ",
  sum(
    target_convergence$
      number_of_regulatory_rnas > 1
  ),
  "\n\n",
  
  "Combined PNG:\n",
  normalizePath(
    output_png,
    winslash = "/"
  ),
  "\n\n",
  
  "Combined TIFF:\n",
  normalizePath(
    output_tiff,
    winslash = "/"
  ),
  "\n\n",
  
  "Editable SVG:\n",
  normalizePath(
    output_svg,
    winslash = "/"
  ),
  "\n\n"
)

print(
  target_convergence |>
    select(
      target,
      number_of_regulatory_rnas,
      linked_regulatory_rnas,
      target_direction,
      target_log2_fold_change,
      total_hybrid_support,
      strongest_confidence_class
    )
)