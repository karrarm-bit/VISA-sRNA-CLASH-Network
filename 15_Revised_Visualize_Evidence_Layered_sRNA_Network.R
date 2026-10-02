# ============================================================
# BacRegRNA Project
# Script 15 REVISED
#
# Evidence-layered visualization of the CLASH-supported
# regulatory RNA-mRNA network in VISA Staphylococcus aureus
#
# IMPORTANT PRINCIPLES
# ------------------------------------------------------------
# 1. CLASH defines the interaction network.
# 2. Target differential expression is an independent layer.
# 3. No interaction is excluded merely because target mRNA
#    is not significantly differentially expressed.
# 4. "High confidence" terminology is avoided.
# 5. UTR candidates are described conservatively as
#    "putative 3UTR-associated RNA".
# 6. Known sRNAs are benchmarks, not scoring bonuses.
# ============================================================


# ------------------------------------------------------------
# 1. PROJECT PATHS
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

script14_folder <- file.path(
  project_folder,
  "14_revised_evidence_layered_network"
)

script14_table_folder <- file.path(
  script14_folder,
  "tables"
)

network_file <- file.path(
  script14_table_folder,
  "Evidence_Layered_CLASH_Network_ALL.csv"
)

responsive_file <- file.path(
  script14_table_folder,
  "Vancomycin_Responsive_Target_Subset.csv"
)

benchmark_file <- file.path(
  script14_table_folder,
  "Benchmark_vs_nrdF_Descriptive_Comparison.csv"
)

nrdf_file <- file.path(
  script14_table_folder,
  "nrdF_Candidate_Interactions.csv"
)

output_folder <- file.path(
  project_folder,
  "15_revised_network_figures"
)

figure_folder <- file.path(
  output_folder,
  "figures"
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
  figure_folder,
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


# ------------------------------------------------------------
# 2. PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "tidyr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "patchwork",
  "scales"
)

installed_names <- rownames(
  installed.packages()
)

missing_packages <- setdiff(
  required_packages,
  installed_names
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
library(ggrepel)
library(patchwork)
library(scales)


# ------------------------------------------------------------
# 3. VERIFY INPUTS
# ------------------------------------------------------------

required_files <- c(
  network_file,
  responsive_file,
  benchmark_file,
  nrdf_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    paste0(
      "Required Script 14 Revised files missing:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}

cat(
  "\n========================================\n",
  "SCRIPT 15 REVISED\n",
  "EVIDENCE-LAYERED NETWORK VISUALIZATION\n",
  "========================================\n"
)


# ------------------------------------------------------------
# 4. READ DATA
# ------------------------------------------------------------

network <- read_csv(
  network_file,
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

nrdf <- read_csv(
  nrdf_file,
  show_col_types = FALSE
)

cat(
  "\nInput dataset sizes:\n"
)

cat(
  "All CLASH-supported edges:",
  nrow(network),
  "\n"
)

cat(
  "Responsive-target subset:",
  nrow(responsive),
  "\n"
)

cat(
  "Benchmark/nrdF rows:",
  nrow(benchmark),
  "\n"
)

cat(
  "nrdF candidate rows:",
  nrow(nrdf),
  "\n"
)


# ------------------------------------------------------------
# 5. HELPERS
# ------------------------------------------------------------

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
  
  value <- str_to_lower(
    str_squish(
      as.character(x)
    )
  )
  
  value %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}


short_label <- function(
    x,
    max_length = 25
) {
  
  x <- safe_chr(x)
  
  ifelse(
    str_length(x) > max_length,
    paste0(
      str_sub(
        x,
        1,
        max_length - 3
      ),
      "..."
    ),
    x
  )
}


# ------------------------------------------------------------
# 6. VERIFY CRITICAL COLUMNS
# ------------------------------------------------------------

required_network_columns <- c(
  "source_node",
  "source_label",
  "source_raw_names",
  "source_RNA_class",
  "target_node",
  "target_label",
  "target_raw_names",
  "number_of_supporting_rows",
  "total_hybrid_count",
  "maximum_number_of_experiments",
  "best_adjusted_p_value",
  "best_connection_score",
  "clash_evidence_score",
  "clash_support_category",
  "target_log2_fold_change",
  "target_deseq_adjusted_p",
  "target_is_significant_deg",
  "target_strong_response",
  "target_transcriptomic_layer",
  "benchmark_annotation",
  "evidence_profile"
)

missing_columns <- setdiff(
  required_network_columns,
  names(network)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "Missing required columns in revised Script 14 output:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 7. STANDARDIZE MAIN NETWORK
# ------------------------------------------------------------

network <- network |>
  mutate(
    
    source_node =
      safe_chr(source_node),
    
    target_node =
      safe_chr(target_node),
    
    source_label =
      safe_chr(source_label),
    
    target_label =
      safe_chr(target_label),
    
    source_raw_names =
      safe_chr(source_raw_names),
    
    target_raw_names =
      safe_chr(target_raw_names),
    
    source_RNA_class =
      safe_chr(source_RNA_class),
    
    clash_support_category =
      safe_chr(clash_support_category),
    
    target_transcriptomic_layer =
      safe_chr(target_transcriptomic_layer),
    
    benchmark_annotation =
      safe_chr(benchmark_annotation),
    
    evidence_profile =
      safe_chr(evidence_profile),
    
    number_of_supporting_rows =
      safe_num(number_of_supporting_rows),
    
    total_hybrid_count =
      safe_num(total_hybrid_count),
    
    maximum_number_of_experiments =
      safe_num(maximum_number_of_experiments),
    
    best_adjusted_p_value =
      safe_num(best_adjusted_p_value),
    
    best_connection_score =
      safe_num(best_connection_score),
    
    clash_evidence_score =
      safe_num(clash_evidence_score),
    
    target_log2_fold_change =
      safe_num(target_log2_fold_change),
    
    target_deseq_adjusted_p =
      safe_num(target_deseq_adjusted_p),
    
    target_is_significant_deg =
      safe_logical(target_is_significant_deg),
    
    target_strong_response =
      safe_logical(target_strong_response)
  )


# ------------------------------------------------------------
# 8. BUILD ROBUST DISPLAY LABELS
# ------------------------------------------------------------

network <- network |>
  mutate(
    
    source_display =
      case_when(
        
        source_label != "" ~
          source_label,
        
        source_raw_names != "" ~
          source_raw_names,
        
        TRUE ~
          source_node
      ),
    
    target_display =
      case_when(
        
        target_label != "" ~
          target_label,
        
        target_raw_names != "" ~
          target_raw_names,
        
        TRUE ~
          target_node
      ),
    
    source_display =
      short_label(
        source_display,
        26
      ),
    
    target_display =
      short_label(
        target_display,
        24
      )
  )


# ------------------------------------------------------------
# 9. CONSERVATIVE RNA CLASS LABELS
# ------------------------------------------------------------

network <- network |>
  mutate(
    
    source_class_display =
      case_when(
        
        str_detect(
          source_RNA_class,
          regex(
            "Named sRNA",
            ignore_case = TRUE
          )
        ) ~
          "Named sRNA",
        
        str_detect(
          source_RNA_class,
          regex(
            "3UTR",
            ignore_case = TRUE
          )
        ) ~
          "Putative 3′UTR-associated RNA",
        
        str_detect(
          source_RNA_class,
          regex(
            "5UTR",
            ignore_case = TRUE
          )
        ) ~
          "Putative 5′UTR-associated RNA",
        
        str_detect(
          source_RNA_class,
          regex(
            "Intergenic",
            ignore_case = TRUE
          )
        ) ~
          "Intergenic RNA candidate",
        
        str_detect(
          source_RNA_class,
          regex(
            "Numbered",
            ignore_case = TRUE
          )
        ) ~
          "Numbered sRNA candidate",
        
        TRUE ~
          "Other regulatory RNA candidate"
      )
  )


# ------------------------------------------------------------
# 10. TARGET TRANSCRIPTOMIC STATUS
# ------------------------------------------------------------

network <- network |>
  mutate(
    
    target_status =
      case_when(
        
        target_strong_response ~
          "Strong target DEG",
        
        target_is_significant_deg ~
          "Significant target DEG",
        
        !is.na(target_log2_fold_change) ~
          "Measured, below DEG threshold",
        
        TRUE ~
          "No mapped transcriptomic evidence"
      )
  )


# ------------------------------------------------------------
# 11. TARGET RESPONSE DIRECTION
# ------------------------------------------------------------

network <- network |>
  mutate(
    
    target_direction =
      case_when(
        
        target_is_significant_deg &
          target_log2_fold_change > 0 ~
          "Induced",
        
        target_is_significant_deg &
          target_log2_fold_change < 0 ~
          "Repressed",
        
        !is.na(target_log2_fold_change) ~
          "Measured non-DEG",
        
        TRUE ~
          "Not mapped"
      )
  )


# ------------------------------------------------------------
# 12. NORMALIZE CLASH SUPPORT ORDER
# ------------------------------------------------------------

network <- network |>
  mutate(
    
    clash_support_category =
      factor(
        clash_support_category,
        levels = c(
          "Higher CLASH support",
          "Intermediate CLASH support",
          "Limited CLASH support",
          "Minimal CLASH support"
        )
      )
  )


# ------------------------------------------------------------
# 13. DEFINE FIGURE SUBSETS
#
# Panel A:
# Higher + Intermediate CLASH support
#
# Panel B:
# All significant target-DE interactions
#
# Panel C:
# RsaOI + SprA2/RsaJ + nrdF candidates
# ------------------------------------------------------------

panel_A_edges <- network |>
  filter(
    clash_support_category %in% c(
      "Higher CLASH support",
      "Intermediate CLASH support"
    )
  )


panel_B_edges <- network |>
  filter(
    target_is_significant_deg
  )


panel_C_edges <- network |>
  filter(
    
    str_detect(
      str_to_lower(
        paste(
          source_display,
          source_node,
          source_raw_names
        )
      ),
      "rsaoi|sau[-_ ]?6477|spra2|rsaj"
    ) |
      
      str_detect(
        str_to_lower(
          paste(
            target_display,
            target_node,
            target_raw_names
          )
        ),
        "\\bnrdf\\b"
      )
  )


# ------------------------------------------------------------
# 14. CHECK SUBSETS
# ------------------------------------------------------------

if (nrow(panel_A_edges) == 0) {
  
  stop(
    "Panel A contains zero interactions."
  )
}

if (nrow(panel_B_edges) == 0) {
  
  stop(
    "Panel B contains zero target-DE interactions."
  )
}

if (nrow(panel_C_edges) == 0) {
  
  stop(
    "Panel C contains zero benchmark/candidate interactions."
  )
}


cat(
  "\nFigure subsets:\n"
)

cat(
  "Panel A — Higher/intermediate CLASH support:",
  nrow(panel_A_edges),
  "\n"
)

cat(
  "Panel B — Significant target DEG:",
  nrow(panel_B_edges),
  "\n"
)

cat(
  "Panel C — Benchmark + nrdF:",
  nrow(panel_C_edges),
  "\n"
)


# ------------------------------------------------------------
# 15. CREATE UNIQUE EDGE IDs
# ------------------------------------------------------------

network <- network |>
  mutate(
    
    edge_id =
      paste0(
        source_node,
        "__TO__",
        target_node
      )
  )

panel_A_edges <- panel_A_edges |>
  mutate(
    edge_id =
      paste0(
        source_node,
        "__TO__",
        target_node
      )
  )

panel_B_edges <- panel_B_edges |>
  mutate(
    edge_id =
      paste0(
        source_node,
        "__TO__",
        target_node
      )
  )

panel_C_edges <- panel_C_edges |>
  mutate(
    edge_id =
      paste0(
        source_node,
        "__TO__",
        target_node
      )
  )


# ------------------------------------------------------------
# 16. BUILD MANUAL BIPARTITE LAYOUT
# ------------------------------------------------------------

build_layout <- function(edge_data) {
  
  source_summary <- edge_data |>
    group_by(
      source_node,
      source_display,
      source_class_display
    ) |>
    summarise(
      
      maximum_CLASH_score =
        max(
          clash_evidence_score,
          na.rm = TRUE
        ),
      
      total_hybrids =
        sum(
          total_hybrid_count,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    ) |>
    arrange(
      desc(maximum_CLASH_score),
      desc(total_hybrids),
      source_display
    )
  
  
  target_summary <- edge_data |>
    group_by(
      target_node,
      target_display
    ) |>
    summarise(
      
      target_status =
        first(target_status),
      
      target_direction =
        first(target_direction),
      
      target_log2FC =
        first(target_log2_fold_change),
      
      maximum_CLASH_score =
        max(
          clash_evidence_score,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    ) |>
    arrange(
      
      desc(
        target_status %in%
          c(
            "Strong target DEG",
            "Significant target DEG"
          )
      ),
      
      desc(
        abs(
          ifelse(
            is.na(target_log2FC),
            0,
            target_log2FC
          )
        )
      ),
      
      desc(maximum_CLASH_score),
      
      target_display
    )
  
  
  source_positions <- source_summary |>
    mutate(
      
      x = 0,
      
      y =
        seq(
          from = nrow(source_summary),
          to = 1,
          length.out =
            nrow(source_summary)
        ),
      
      node_type =
        "Regulatory RNA",
      
      node_class =
        source_class_display,
      
      transcriptomic_status =
        NA_character_,
      
      response_direction =
        NA_character_
    ) |>
    transmute(
      
      node_id =
        source_node,
      
      display_label =
        source_display,
      
      x,
      y,
      node_type,
      node_class,
      transcriptomic_status,
      response_direction
    )
  
  
  target_positions <- target_summary |>
    mutate(
      
      x = 1,
      
      y =
        seq(
          from = nrow(target_summary),
          to = 1,
          length.out =
            nrow(target_summary)
        ),
      
      node_type =
        "mRNA target",
      
      node_class =
        "mRNA target",
      
      transcriptomic_status =
        target_status,
      
      response_direction =
        target_direction
    ) |>
    transmute(
      
      node_id =
        target_node,
      
      display_label =
        target_display,
      
      x,
      y,
      node_type,
      node_class,
      transcriptomic_status,
      response_direction
    )
  
  
  nodes <- bind_rows(
    source_positions,
    target_positions
  )
  
  
  plot_edges <- edge_data |>
    left_join(
      
      nodes |>
        select(
          source_node = node_id,
          x_from = x,
          y_from = y
        ),
      
      by = "source_node"
    ) |>
    
    left_join(
      
      nodes |>
        select(
          target_node = node_id,
          x_to = x,
          y_to = y
        ),
      
      by = "target_node"
    )
  
  
  list(
    nodes = nodes,
    edges = plot_edges
  )
}


# ------------------------------------------------------------
# 17. BUILD THREE LAYOUTS
# ------------------------------------------------------------

layout_A <- build_layout(
  panel_A_edges
)

layout_B <- build_layout(
  panel_B_edges
)

layout_C <- build_layout(
  panel_C_edges
)


# ------------------------------------------------------------
# 18. PLOT FUNCTION
# ------------------------------------------------------------

create_network_plot <- function(
    layout_data,
    title_text,
    subtitle_text
) {
  
  edges_plot <- layout_data$edges |>
    mutate(
      
      edge_width =
        sqrt(
          pmax(
            total_hybrid_count,
            1
          )
        ),
      
      support =
        factor(
          clash_support_category,
          levels = c(
            "Higher CLASH support",
            "Intermediate CLASH support",
            "Limited CLASH support",
            "Minimal CLASH support"
          )
        )
    )
  
  
  nodes_plot <- layout_data$nodes |>
    mutate(
      
      plot_group =
        case_when(
          
          node_type ==
            "Regulatory RNA" ~
            node_class,
          
          response_direction ==
            "Induced" ~
            "Induced target",
          
          response_direction ==
            "Repressed" ~
            "Repressed target",
          
          response_direction ==
            "Measured non-DEG" ~
            "Measured non-DEG target",
          
          TRUE ~
            "Target without mapped transcriptomic evidence"
        ),
      
      plot_shape =
        case_when(
          
          node_type ==
            "mRNA target" ~
            "mRNA target",
          
          node_class ==
            "Named sRNA" ~
            "Named sRNA",
          
          str_detect(
            node_class,
            "3′UTR"
          ) ~
            "Putative 3′UTR-associated RNA",
          
          str_detect(
            node_class,
            "5′UTR"
          ) ~
            "Putative 5′UTR-associated RNA",
          
          str_detect(
            node_class,
            "Intergenic"
          ) ~
            "Intergenic RNA candidate",
          
          str_detect(
            node_class,
            "Numbered"
          ) ~
            "Numbered sRNA candidate",
          
          TRUE ~
            "Other regulatory RNA"
        )
    )
  
  
  maximum_y <- max(
    nodes_plot$y,
    na.rm = TRUE
  )
  
  
  ggplot() +
    
    geom_curve(
      
      data = edges_plot,
      
      aes(
        x = x_from,
        y = y_from,
        xend = x_to,
        yend = y_to,
        linewidth = edge_width,
        linetype = support
      ),
      
      curvature = 0.10,
      alpha = 0.72,
      lineend = "round",
      
      arrow = grid::arrow(
        type = "closed",
        length =
          grid::unit(
            2,
            "mm"
          )
      )
    ) +
    
    geom_point(
      
      data = nodes_plot,
      
      aes(
        x = x,
        y = y,
        shape = plot_shape,
        fill = plot_group
      ),
      
      size = 4.8,
      stroke = 0.55
    ) +
    
    geom_text(
      
      data =
        nodes_plot |>
        filter(
          x == 0
        ),
      
      aes(
        x = x,
        y = y,
        label = display_label
      ),
      
      hjust = 1.15,
      size = 3.1,
      fontface = "bold"
    ) +
    
    geom_text(
      
      data =
        nodes_plot |>
        filter(
          x == 1
        ),
      
      aes(
        x = x,
        y = y,
        label = display_label
      ),
      
      hjust = -0.15,
      size = 3.1
    ) +
    
    annotate(
      
      "text",
      
      x = 0,
      y = maximum_y + 0.75,
      
      label =
        "Regulatory RNAs",
      
      fontface = "bold",
      size = 3.8
    ) +
    
    annotate(
      
      "text",
      
      x = 1,
      y = maximum_y + 0.75,
      
      label =
        "mRNA targets",
      
      fontface = "bold",
      size = 3.8
    ) +
    
    scale_shape_manual(
      
      values = c(
        
        "Named sRNA" = 22,
        
        "Putative 3′UTR-associated RNA" = 21,
        
        "Putative 5′UTR-associated RNA" = 24,
        
        "Intergenic RNA candidate" = 23,
        
        "Numbered sRNA candidate" = 25,
        
        "Other regulatory RNA" = 24,
        
        "mRNA target" = 21
      ),
      
      name =
        "Node type"
    ) +
    
    scale_linetype_manual(
      
      values = c(
        
        "Higher CLASH support" =
          "solid",
        
        "Intermediate CLASH support" =
          "solid",
        
        "Limited CLASH support" =
          "longdash",
        
        "Minimal CLASH support" =
          "dotted"
      ),
      
      name =
        "CLASH evidence"
    ) +
    
    scale_linewidth_continuous(
      
      range = c(
        0.35,
        2.0
      ),
      
      guide =
        "none"
    ) +
    
    scale_fill_discrete(
      name =
        "RNA / target status"
    ) +
    
    coord_cartesian(
      
      xlim = c(
        -0.48,
        1.48
      ),
      
      clip = "off"
    ) +
    
    labs(
      
      title =
        title_text,
      
      subtitle =
        subtitle_text
    ) +
    
    theme_void(
      base_size = 10
    ) +
    
    theme(
      
      plot.background =
        element_rect(
          fill = "white",
          colour = NA
        ),
      
      panel.background =
        element_rect(
          fill = "white",
          colour = NA
        ),
      
      legend.background =
        element_rect(
          fill = "white",
          colour = NA
        ),
      
      legend.key =
        element_rect(
          fill = "white",
          colour = NA
        ),
      
      legend.position =
        "bottom",
      
      legend.box =
        "vertical",
      
      plot.title =
        element_text(
          face = "bold",
          size = 12,
          hjust = 0.5
        ),
      
      plot.subtitle =
        element_text(
          size = 9.5,
          hjust = 0.5
        ),
      
      plot.margin =
        margin(
          25,
          120,
          35,
          120,
          unit = "pt"
        )
    )
}


# ------------------------------------------------------------
# 19. PANEL A
# ------------------------------------------------------------

panel_A <- create_network_plot(
  
  layout_A,
  
  title_text =
    "Evidence-supported CLASH regulatory network",
  
  subtitle_text =
    paste0(
      "Higher and intermediate CLASH support; ",
      "target transcriptomic response is not an inclusion criterion ",
      "(n = ",
      nrow(panel_A_edges),
      " interactions)"
    )
)


# ------------------------------------------------------------
# 20. PANEL B
# ------------------------------------------------------------

panel_B <- create_network_plot(
  
  layout_B,
  
  title_text =
    "Vancomycin-responsive target subset",
  
  subtitle_text =
    paste0(
      "CLASH-supported interactions whose mRNA targets ",
      "meet the prespecified differential-expression criterion ",
      "(n = ",
      nrow(panel_B_edges),
      " interactions)"
    )
)


# ------------------------------------------------------------
# 21. PANEL C
# ------------------------------------------------------------

panel_C <- create_network_plot(
  
  layout_C,
  
  title_text =
    "Benchmark RNAs and nrdF candidate interactions",
  
  subtitle_text =
    paste0(
      "Descriptive comparison of recovered known sRNAs ",
      "and putative 3′UTR-associated nrdF interactions ",
      "(n = ",
      nrow(panel_C_edges),
      " interactions)"
    )
)


# ------------------------------------------------------------
# 22. DISPLAY
# ------------------------------------------------------------

print(panel_A)

print(panel_B)

print(panel_C)


# ------------------------------------------------------------
# 23. EXPORT INDIVIDUAL FIGURES
# ------------------------------------------------------------

save_panel <- function(
    plot_object,
    base_name,
    width,
    height
) {
  
  ggsave(
    
    filename =
      file.path(
        figure_folder,
        paste0(
          base_name,
          ".pdf"
        )
      ),
    
    plot =
      plot_object,
    
    width =
      width,
    
    height =
      height,
    
    units =
      "in",
    
    device =
      cairo_pdf,
    
    bg =
      "white",
    
    limitsize =
      FALSE
  )
  
  
  ggsave(
    
    filename =
      file.path(
        figure_folder,
        paste0(
          base_name,
          "_600dpi.png"
        )
      ),
    
    plot =
      plot_object,
    
    width =
      width,
    
    height =
      height,
    
    units =
      "in",
    
    dpi =
      600,
    
    bg =
      "white",
    
    limitsize =
      FALSE
  )
  
  
  ggsave(
    
    filename =
      file.path(
        figure_folder,
        paste0(
          base_name,
          "_600dpi.tiff"
        )
      ),
    
    plot =
      plot_object,
    
    width =
      width,
    
    height =
      height,
    
    units =
      "in",
    
    dpi =
      600,
    
    compression =
      "lzw",
    
    bg =
      "white",
    
    limitsize =
      FALSE
  )
}


save_panel(
  panel_A,
  "Figure2A_Evidence_Supported_CLASH_Network",
  11,
  8
)

save_panel(
  panel_B,
  "Figure2B_Vancomycin_Responsive_Target_Subset",
  11,
  7
)

save_panel(
  panel_C,
  "Figure2C_Benchmark_and_nrdF_Context",
  11,
  7
)


# ------------------------------------------------------------
# 24. COMBINED FIGURE
# ------------------------------------------------------------

combined_figure <- (
  
  panel_A /
    
    panel_B /
    
    panel_C
  
) +
  
  patchwork::plot_layout(
    
    heights =
      c(
        1.15,
        1,
        1
      ),
    
    guides =
      "collect"
  ) +
  
  patchwork::plot_annotation(
    
    tag_levels =
      "A",
    
    theme =
      theme(
        
        plot.background =
          element_rect(
            fill = "white",
            colour = NA
          ),
        
        plot.tag =
          element_text(
            face = "bold",
            size = 17
          )
      )
  ) &
  
  theme(
    
    legend.position =
      "bottom",
    
    plot.background =
      element_rect(
        fill = "white",
        colour = NA
      )
  )


print(
  combined_figure
)


# ------------------------------------------------------------
# 25. EXPORT COMBINED FIGURE
# ------------------------------------------------------------

ggsave(
  
  file.path(
    figure_folder,
    "Figure2_ABC_Evidence_Layered_sRNA_Network.pdf"
  ),
  
  combined_figure,
  
  width = 12,
  height = 21,
  
  units = "in",
  
  device = cairo_pdf,
  
  bg = "white",
  
  limitsize = FALSE
)


ggsave(
  
  file.path(
    figure_folder,
    "Figure2_ABC_Evidence_Layered_sRNA_Network_600dpi.png"
  ),
  
  combined_figure,
  
  width = 12,
  height = 21,
  
  units = "in",
  
  dpi = 600,
  
  bg = "white",
  
  limitsize = FALSE
)


ggsave(
  
  file.path(
    figure_folder,
    "Figure2_ABC_Evidence_Layered_sRNA_Network_600dpi.tiff"
  ),
  
  combined_figure,
  
  width = 12,
  height = 21,
  
  units = "in",
  
  dpi = 600,
  
  compression = "lzw",
  
  bg = "white",
  
  limitsize = FALSE
)


# ------------------------------------------------------------
# 26. SAVE PANEL TABLES
# ------------------------------------------------------------

write_csv(
  
  layout_A$edges,
  
  file.path(
    table_folder,
    "Figure2A_edges.csv"
  )
)

write_csv(
  
  layout_A$nodes,
  
  file.path(
    table_folder,
    "Figure2A_nodes.csv"
  )
)


write_csv(
  
  layout_B$edges,
  
  file.path(
    table_folder,
    "Figure2B_edges.csv"
  )
)

write_csv(
  
  layout_B$nodes,
  
  file.path(
    table_folder,
    "Figure2B_nodes.csv"
  )
)


write_csv(
  
  layout_C$edges,
  
  file.path(
    table_folder,
    "Figure2C_edges.csv"
  )
)

write_csv(
  
  layout_C$nodes,
  
  file.path(
    table_folder,
    "Figure2C_nodes.csv"
  )
)


# ------------------------------------------------------------
# 27. FIGURE AUDIT SUMMARY
# ------------------------------------------------------------

figure_summary <- tibble(
  
  metric = c(
    
    "Total CLASH-supported network edges",
    
    "Higher/intermediate CLASH support edges",
    
    "Vancomycin-responsive target edges",
    
    "Benchmark/nrdF context edges",
    
    "Total unique regulatory RNAs",
    
    "Total unique mRNA targets",
    
    "Responsive-subset regulatory RNAs",
    
    "Responsive-subset mRNA targets",
    
    "nrdF candidate interactions"
  ),
  
  value = c(
    
    nrow(network),
    
    nrow(panel_A_edges),
    
    nrow(panel_B_edges),
    
    nrow(panel_C_edges),
    
    n_distinct(
      network$source_node
    ),
    
    n_distinct(
      network$target_node
    ),
    
    n_distinct(
      panel_B_edges$source_node
    ),
    
    n_distinct(
      panel_B_edges$target_node
    ),
    
    nrow(nrdf)
  )
)


write_csv(
  
  figure_summary,
  
  file.path(
    audit_folder,
    "Script15_Revised_Figure_Summary.csv"
  )
)


# ------------------------------------------------------------
# 28. CRITICAL INTERACTION AUDIT
# ------------------------------------------------------------

critical_audit <- network |>
  filter(
    
    str_detect(
      str_to_lower(
        paste(
          source_display,
          source_node,
          source_raw_names
        )
      ),
      "rsaoi|sau[-_ ]?6477|spra2|rsaj"
    ) |
      
      str_detect(
        str_to_lower(
          paste(
            target_display,
            target_node,
            target_raw_names
          )
        ),
        "\\bnrdf\\b"
      )
  ) |>
  
  select(
    
    source_display,
    
    target_display,
    
    source_class_display,
    
    total_hybrid_count,
    
    maximum_number_of_experiments,
    
    best_adjusted_p_value,
    
    best_connection_score,
    
    clash_evidence_score,
    
    clash_support_category,
    
    target_log2_fold_change,
    
    target_deseq_adjusted_p,
    
    target_is_significant_deg,
    
    target_status,
    
    benchmark_annotation,
    
    evidence_profile
  )


write_csv(
  
  critical_audit,
  
  file.path(
    audit_folder,
    "Critical_Interaction_Figure_Audit.csv"
  )
)


# ------------------------------------------------------------
# 29. MANUSCRIPT-SAFE FIGURE NOTES
# ------------------------------------------------------------

figure_notes <- c(
  
  "SCRIPT 15 REVISED — FIGURE INTERPRETATION",
  
  "",
  
  paste0(
    "1. CLASH-supported regulatory RNA-mRNA interactions ",
    "define the interaction network."
  ),
  
  paste0(
    "2. Target differential expression is represented as ",
    "an independent evidence layer and is not required for ",
    "retention of a CLASH-supported interaction."
  ),
  
  paste0(
    "3. Panel A displays interactions with higher or ",
    "intermediate CLASH support."
  ),
  
  paste0(
    "4. Panel B displays the subset whose mRNA targets meet ",
    "the prespecified vancomycin differential-expression criterion."
  ),
  
  paste0(
    "5. Panel C provides descriptive context for recovered ",
    "benchmark sRNAs and the putative 3UTR-associated ",
    "RNA-nrdF interactions."
  ),
  
  paste0(
    "6. Edge width reflects CLASH hybrid support."
  ),
  
  paste0(
    "7. Edge line type represents descriptive CLASH-support ",
    "categories rather than experimentally validated ",
    "regulatory confidence."
  ),
  
  paste0(
    "8. Arrow direction represents analytical assignment ",
    "from regulatory RNA to mRNA target and does not establish ",
    "activation, repression, or causal direction."
  ),
  
  paste0(
    "9. Putative 3UTR-associated RNAs are candidate ",
    "regulatory RNA species; their independent production ",
    "or processing is not established by CLASH alone."
  ),
  
  paste0(
    "10. Target mRNA differential abundance provides ",
    "orthogonal transcriptomic context and should not be ",
    "interpreted as proof of direct regulation."
  )
)


writeLines(
  
  figure_notes,
  
  file.path(
    audit_folder,
    "Script15_Revised_Figure_Notes.txt"
  )
)


# ------------------------------------------------------------
# 30. SESSION INFO
# ------------------------------------------------------------

capture.output(
  
  sessionInfo(),
  
  file =
    file.path(
      audit_folder,
      "sessionInfo.txt"
    )
)


# ------------------------------------------------------------
# 31. FINAL REPORT
# ------------------------------------------------------------

cat(
  "\n========================================\n"
)

cat(
  "SCRIPT 15 REVISED COMPLETED SUCCESSFULLY\n"
)

cat(
  "========================================\n"
)

cat(
  "\nOutput folder:\n",
  output_folder,
  "\n"
)

cat(
  "\nFIGURE SUMMARY\n"
)

print(
  figure_summary,
  n = Inf
)

cat(
  "\nCRITICAL INTERACTIONS\n"
)

print(
  critical_audit,
  n = Inf
)

cat(
  "\nMost important outputs:\n"
)

cat(
  "1. figures/Figure2A_Evidence_Supported_CLASH_Network_600dpi.png\n"
)

cat(
  "2. figures/Figure2B_Vancomycin_Responsive_Target_Subset_600dpi.png\n"
)

cat(
  "3. figures/Figure2C_Benchmark_and_nrdF_Context_600dpi.png\n"
)

cat(
  "4. figures/Figure2_ABC_Evidence_Layered_sRNA_Network_600dpi.png\n"
)

cat(
  "5. audit/Script15_Revised_Figure_Summary.csv\n"
)

cat(
  "6. audit/Critical_Interaction_Figure_Audit.csv\n"
)

cat(
  "\nNEXT STEP:\n"
)

cat(
  paste0(
    "Inspect the figure summary and critical interaction audit ",
    "before modifying downstream functional-annotation scripts.\n"
  )
)