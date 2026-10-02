# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages
BiocManager::install(
  c("DESeq2", "SummarizedExperiment"),
  ask = FALSE,
  update = FALSE
)
library(DESeq2)
library(SummarizedExperiment)

packageVersion("DESeq2")
packageVersion("SummarizedExperiment")
# ============================================================
# BacRegRNA Project
# Publication-ready Figure 1
#
# Panel A: PCA
# Panel B: Sample-distance heatmap
# Panel C: Volcano plot
# Panel D: Top-24 annotated DEG heatmap
#
# Main improvements:
# - Automatic installation of missing CRAN packages
# - Automatic installation of DESeq2 and SummarizedExperiment
#   through Bioconductor
# - Protected separation between Panels A/B and C/D
# - Additional space for heatmap labels and legends
# - Export as PDF, TIFF 600 dpi, and PNG 600 dpi
# ============================================================


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

vst_file <- file.path(
  project_folder,
  "09_RNAseq_analysis",
  "GSE254530_VST_object.rds"
)

annotated_results_file <- file.path(
  project_folder,
  "10_annotation",
  "tables",
  "GSE254530_DESeq2_annotated_all_genes.csv"
)

annotated_significant_file <- file.path(
  project_folder,
  "10_annotation",
  "tables",
  "GSE254530_DESeq2_annotated_significant_genes.csv"
)

top30_file <- file.path(
  project_folder,
  "10_annotation",
  "tables",
  "GSE254530_top30_annotated_heatmap_genes.csv"
)

output_folder <- file.path(
  project_folder,
  "11_publication_figures_final"
)

panel_folder <- file.path(
  output_folder,
  "individual_panels"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  panel_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Install and load required packages
# ------------------------------------------------------------

cran_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "patchwork",
  "ggplotify",
  "scales"
)

bioconductor_packages <- c(
  "DESeq2",
  "SummarizedExperiment"
)


# Install missing CRAN packages
missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_cran) > 0) {
  
  message(
    "Installing missing CRAN packages: ",
    paste(missing_cran, collapse = ", ")
  )
  
  install.packages(
    missing_cran,
    dependencies = TRUE
  )
}


# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  
  message("Installing BiocManager...")
  
  install.packages(
    "BiocManager"
  )
}


# Install missing Bioconductor packages
missing_bioc <- bioconductor_packages[
  !vapply(
    bioconductor_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_bioc) > 0) {
  
  message(
    "Installing missing Bioconductor packages: ",
    paste(missing_bioc, collapse = ", ")
  )
  
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}


# Confirm all packages are available
all_required_packages <- c(
  cran_packages,
  bioconductor_packages
)

still_missing <- all_required_packages[
  !vapply(
    all_required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(still_missing) > 0) {
  
  stop(
    paste0(
      "Installation did not complete for: ",
      paste(still_missing, collapse = ", "),
      "\nRestart RStudio and run the script again."
    )
  )
}


# Load packages
suppressPackageStartupMessages({
  
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(patchwork)
  library(ggplotify)
  library(scales)
  library(DESeq2)
  library(SummarizedExperiment)
  
})


# ------------------------------------------------------------
# 3. Check required input files
# ------------------------------------------------------------

required_files <- c(
  vst_file,
  annotated_results_file,
  annotated_significant_file,
  top30_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    paste0(
      "The following input files were not found:\n\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


# ------------------------------------------------------------
# 4. Read data
# ------------------------------------------------------------

vsd <- readRDS(
  vst_file
)

vst_matrix <- SummarizedExperiment::assay(
  vsd
)

annotated_results <- read_csv(
  annotated_results_file,
  show_col_types = FALSE
)

annotated_significant <- read_csv(
  annotated_significant_file,
  show_col_types = FALSE
)

top30_genes <- read_csv(
  top30_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 5. Validate important columns
# ------------------------------------------------------------

required_results_columns <- c(
  "gene_id",
  "log2_fold_change",
  "adjusted_p_value",
  "significance"
)

missing_results_columns <- setdiff(
  required_results_columns,
  names(annotated_results)
)

if (length(missing_results_columns) > 0) {
  
  stop(
    paste0(
      "Missing columns in annotated results file: ",
      paste(
        missing_results_columns,
        collapse = ", "
      )
    )
  )
}

required_top30_columns <- c(
  "gene_id",
  "log2_fold_change",
  "adjusted_p_value",
  "significance"
)

missing_top30_columns <- setdiff(
  required_top30_columns,
  names(top30_genes)
)

if (length(missing_top30_columns) > 0) {
  
  stop(
    paste0(
      "Missing columns in the top-30 heatmap file: ",
      paste(
        missing_top30_columns,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 6. Short sample labels
# ------------------------------------------------------------

short_sample_label <- function(x) {
  
  case_when(
    
    str_detect(
      x,
      regex(
        "Cont1",
        ignore_case = TRUE
      )
    ) ~ "C1",
    
    str_detect(
      x,
      regex(
        "Cont2",
        ignore_case = TRUE
      )
    ) ~ "C2",
    
    str_detect(
      x,
      regex(
        "Cont3",
        ignore_case = TRUE
      )
    ) ~ "C3",
    
    str_detect(
      x,
      regex(
        "Van1",
        ignore_case = TRUE
      )
    ) ~ "V1",
    
    str_detect(
      x,
      regex(
        "Van2",
        ignore_case = TRUE
      )
    ) ~ "V2",
    
    str_detect(
      x,
      regex(
        "Van3",
        ignore_case = TRUE
      )
    ) ~ "V3",
    
    TRUE ~ x
  )
}

original_names <- colnames(
  vst_matrix
)

short_names <- vapply(
  original_names,
  short_sample_label,
  FUN.VALUE = character(1)
)

if (anyDuplicated(short_names) > 0) {
  
  stop(
    "Duplicated short sample labels were detected."
  )
}

sample_key <- tibble(
  original_sample = original_names,
  sample = short_names,
  condition = if_else(
    str_detect(
      short_names,
      "^C"
    ),
    "Untreated",
    "Treated"
  )
)

sample_key$condition <- factor(
  sample_key$condition,
  levels = c(
    "Untreated",
    "Treated"
  )
)

colnames(vst_matrix) <- short_names


# ------------------------------------------------------------
# 7. Colours and shapes
# ------------------------------------------------------------

condition_colours <- c(
  "Untreated" = "#0072B2",
  "Treated" = "#D55E00"
)

regulation_colours <- c(
  "Downregulated" = "#0072B2",
  "Not significant" = "grey80",
  "Upregulated" = "#D55E00"
)

condition_shapes <- c(
  "Untreated" = 16,
  "Treated" = 17
)


# ------------------------------------------------------------
# 8. Publication theme
# ------------------------------------------------------------

publication_theme <- theme_classic(
  base_size = 10,
  base_family = "sans"
) +
  theme(
    
    axis.title = element_text(
      size = 10,
      colour = "black"
    ),
    
    axis.text = element_text(
      size = 9,
      colour = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.35
    ),
    
    axis.ticks = element_line(
      linewidth = 0.30
    ),
    
    legend.text = element_text(
      size = 8.5
    ),
    
    legend.title = element_text(
      size = 8.5,
      face = "bold"
    ),
    
    legend.key.height = grid::unit(
      0.32,
      "cm"
    ),
    
    legend.key.width = grid::unit(
      0.42,
      "cm"
    ),
    
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )


# ------------------------------------------------------------
# 9. Gene-label cleaning function
# ------------------------------------------------------------

clean_gene_label <- function(
    gene_id,
    gene_symbol,
    product_short,
    product
) {
  
  gene_symbol <- na_if(
    str_squish(
      as.character(
        gene_symbol
      )
    ),
    ""
  )
  
  product_short <- na_if(
    str_squish(
      as.character(
        product_short
      )
    ),
    ""
  )
  
  product <- na_if(
    str_squish(
      as.character(
        product
      )
    ),
    ""
  )
  
  selected_product <- coalesce(
    product_short,
    product
  )
  
  selected_product <- case_when(
    
    str_detect(
      selected_product,
      regex(
        "^LCP family glycopolymer transferase",
        ignore_case = TRUE
      )
    ) ~ "LCP glycopolymer transferase",
    
    str_detect(
      selected_product,
      regex(
        "beta-class phenol-soluble modulin",
        ignore_case = TRUE
      )
    ) ~ "β-type phenol-soluble modulin",
    
    str_detect(
      selected_product,
      regex(
        "ABC transporter ATP-binding protein",
        ignore_case = TRUE
      )
    ) ~ "ABC transporter ATPase",
    
    str_detect(
      selected_product,
      regex(
        "ABC transporter permease",
        ignore_case = TRUE
      )
    ) ~ "ABC transporter permease",
    
    str_detect(
      selected_product,
      regex(
        "maltodextrin ABC transporter",
        ignore_case = TRUE
      )
    ) ~ "Maltodextrin transporter",
    
    str_detect(
      selected_product,
      regex(
        "thermonuclease",
        ignore_case = TRUE
      )
    ) ~ "Thermonuclease",
    
    str_detect(
      selected_product,
      regex(
        "^hypothetical protein",
        ignore_case = TRUE
      )
    ) ~ gene_id,
    
    str_detect(
      selected_product,
      regex(
        "^membrane protein$",
        ignore_case = TRUE
      )
    ) ~ gene_id,
    
    TRUE ~ selected_product
  )
  
  label <- coalesce(
    gene_symbol,
    selected_product,
    gene_id
  )
  
  if_else(
    str_length(label) > 31,
    paste0(
      str_sub(
        label,
        1,
        28
      ),
      "..."
    ),
    label
  )
}


# ============================================================
# PANEL A: PCA
# ============================================================

pca_data <- DESeq2::plotPCA(
  vsd,
  intgroup = "condition",
  ntop = 500,
  returnData = TRUE
)

percent_variance <- round(
  100 * attr(
    pca_data,
    "percentVar"
  ),
  1
)

pca_data <- pca_data |>
  rownames_to_column(
    "original_sample"
  ) |>
  mutate(
    
    sample = vapply(
      original_sample,
      short_sample_label,
      FUN.VALUE = character(1)
    )
    
  ) |>
  select(
    sample,
    PC1,
    PC2
  ) |>
  left_join(
    
    sample_key |>
      select(
        sample,
        condition
      ),
    
    by = "sample"
  )

panel_A <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    colour = condition,
    shape = condition,
    label = sample
  )
) +
  
  geom_point(
    size = 3.8,
    stroke = 0.55
  ) +
  
  geom_text_repel(
    size = 3.2,
    box.padding = 0.85,
    point.padding = 0.55,
    min.segment.length = 0,
    segment.size = 0.25,
    max.overlaps = Inf,
    seed = 123
  ) +
  
  scale_colour_manual(
    values = condition_colours
  ) +
  
  scale_shape_manual(
    values = condition_shapes
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.15,
        0.18
      )
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.18,
        0.18
      )
    )
  ) +
  
  labs(
    x = paste0(
      "PC1 (",
      percent_variance[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      percent_variance[2],
      "%)"
    ),
    colour = NULL,
    shape = NULL
  ) +
  
  publication_theme +
  
  theme(
    legend.position = "top",
    legend.justification = "left",
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )


# ============================================================
# PANEL B: SAMPLE-DISTANCE HEATMAP
# ============================================================

sample_distance_matrix <- as.matrix(
  dist(
    t(vst_matrix)
  )
)

rownames(sample_distance_matrix) <- short_names
colnames(sample_distance_matrix) <- short_names

sample_annotation <- data.frame(
  Condition = sample_key$condition
)

rownames(sample_annotation) <- sample_key$sample

distance_heatmap <- pheatmap(
  
  sample_distance_matrix,
  
  annotation_col = sample_annotation,
  
  annotation_colors = list(
    Condition = condition_colours
  ),
  
  annotation_names_col = FALSE,
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  border_color = "white",
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  fontsize = 8.5,
  fontsize_row = 8.5,
  fontsize_col = 8.5,
  
  angle_col = 0,
  
  legend = TRUE,
  annotation_legend = TRUE,
  
  treeheight_row = 26,
  treeheight_col = 26,
  
  cellwidth = 35,
  cellheight = 35,
  
  silent = TRUE
)

panel_B_raw <- ggplotify::as.ggplot(
  distance_heatmap$gtable
)

panel_B <- patchwork::wrap_elements(
  full = panel_B_raw,
  clip = TRUE
)


# ============================================================
# PANEL C: VOLCANO PLOT
# ============================================================

volcano_data <- annotated_results |>
  mutate(
    
    adjusted_p_plot = case_when(
      
      is.na(adjusted_p_value) ~ 1,
      
      adjusted_p_value <= 0 ~
        .Machine$double.xmin,
      
      TRUE ~ adjusted_p_value
    ),
    
    minus_log10_adjusted_p =
      -log10(
        adjusted_p_plot
      ),
    
    volcano_class = case_when(
      
      significance == "Upregulated" ~
        "Upregulated",
      
      significance == "Downregulated" ~
        "Downregulated",
      
      TRUE ~
        "Not significant"
    ),
    
    clean_label = clean_gene_label(
      
      gene_id,
      
      if (
        "gene_symbol" %in%
        names(annotated_results)
      ) {
        gene_symbol
      } else {
        NA_character_
      },
      
      if (
        "product_short" %in%
        names(annotated_results)
      ) {
        product_short
      } else {
        NA_character_
      },
      
      if (
        "product" %in%
        names(annotated_results)
      ) {
        product
      } else {
        NA_character_
      }
    )
  )

volcano_data$volcano_class <- factor(
  volcano_data$volcano_class,
  levels = c(
    "Downregulated",
    "Not significant",
    "Upregulated"
  )
)

volcano_candidates <- volcano_data |>
  filter(
    
    volcano_class %in% c(
      "Upregulated",
      "Downregulated"
    ),
    
    !is.na(
      adjusted_p_value
    )
  ) |>
  mutate(
    
    symbol_available =
      if (
        "gene_symbol" %in%
        names(volcano_data)
      ) {
        
        !is.na(gene_symbol) &
          gene_symbol != ""
        
      } else {
        
        FALSE
        
      }
  )

volcano_labels <- bind_rows(
  
  volcano_candidates |>
    filter(
      volcano_class ==
        "Upregulated"
    ) |>
    arrange(
      desc(
        symbol_available
      ),
      adjusted_p_value,
      desc(
        abs(
          log2_fold_change
        )
      )
    ) |>
    slice_head(
      n = 4
    ),
  
  volcano_candidates |>
    filter(
      volcano_class ==
        "Downregulated"
    ) |>
    arrange(
      desc(
        symbol_available
      ),
      adjusted_p_value,
      desc(
        abs(
          log2_fold_change
        )
      )
    ) |>
    slice_head(
      n = 4
    )
) |>
  distinct(
    gene_id,
    .keep_all = TRUE
  )

maximum_lfc <- ceiling(
  max(
    abs(
      volcano_data$log2_fold_change
    ),
    na.rm = TRUE
  )
)

maximum_lfc <- max(
  maximum_lfc,
  5
)

panel_C <- ggplot(
  volcano_data,
  aes(
    x = log2_fold_change,
    y = minus_log10_adjusted_p,
    colour = volcano_class
  )
) +
  
  geom_point(
    size = 1.25,
    alpha = 0.72
  ) +
  
  geom_vline(
    xintercept = c(
      -1,
      1
    ),
    linetype = "dashed",
    linewidth = 0.28,
    colour = "grey45"
  ) +
  
  geom_hline(
    yintercept = -log10(
      0.05
    ),
    linetype = "dashed",
    linewidth = 0.28,
    colour = "grey45"
  ) +
  
  geom_text_repel(
    data = volcano_labels,
    aes(
      x = log2_fold_change,
      y = minus_log10_adjusted_p,
      label = clean_label
    ),
    inherit.aes = FALSE,
    size = 2.9,
    box.padding = 0.70,
    point.padding = 0.40,
    min.segment.length = 0,
    segment.size = 0.25,
    max.overlaps = Inf,
    seed = 123
  ) +
  
  scale_colour_manual(
    values = regulation_colours,
    drop = FALSE
  ) +
  
  scale_x_continuous(
    limits = c(
      -maximum_lfc,
      maximum_lfc
    ),
    breaks = pretty_breaks(
      n = 5
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.02,
        0.10
      )
    )
  ) +
  
  labs(
    x = expression(
      Log[2]~fold~change
    ),
    y = expression(
      -Log[10]~adjusted~italic(P)
    ),
    colour = NULL
  ) +
  
  publication_theme +
  
  theme(
    legend.position = "top",
    legend.justification = "left",
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )


# ============================================================
# PANEL D: TOP-24 DEG HEATMAP
# ============================================================

top24_genes <- bind_rows(
  
  top30_genes |>
    filter(
      significance ==
        "Upregulated"
    ) |>
    arrange(
      adjusted_p_value,
      desc(
        abs(
          log2_fold_change
        )
      )
    ) |>
    slice_head(
      n = 12
    ),
  
  top30_genes |>
    filter(
      significance ==
        "Downregulated"
    ) |>
    arrange(
      adjusted_p_value,
      desc(
        abs(
          log2_fold_change
        )
      )
    ) |>
    slice_head(
      n = 12
    )
) |>
  
  distinct(
    gene_id,
    .keep_all = TRUE
  ) |>
  
  filter(
    gene_id %in%
      rownames(
        vst_matrix
      )
  )

if (nrow(top24_genes) == 0) {
  
  stop(
    "No heatmap genes were matched to the VST matrix."
  )
}

top24_genes <- top24_genes |>
  mutate(
    
    heatmap_label = clean_gene_label(
      
      gene_id,
      
      if (
        "gene_symbol" %in%
        names(top24_genes)
      ) {
        gene_symbol
      } else {
        NA_character_
      },
      
      if (
        "product_short" %in%
        names(top24_genes)
      ) {
        product_short
      } else {
        NA_character_
      },
      
      if (
        "product" %in%
        names(top24_genes)
      ) {
        product
      } else {
        NA_character_
      }
    )
  ) |>
  
  group_by(
    heatmap_label
  ) |>
  
  mutate(
    
    heatmap_label = if_else(
      
      n() > 1,
      
      paste0(
        heatmap_label,
        " (",
        str_extract(
          gene_id,
          "[0-9]+$"
        ),
        ")"
      ),
      
      heatmap_label
    )
  ) |>
  
  ungroup()

heatmap_matrix_scaled <- vst_matrix[
  top24_genes$gene_id,
  ,
  drop = FALSE
]

heatmap_matrix_scaled <- t(
  scale(
    t(
      heatmap_matrix_scaled
    )
  )
)

heatmap_matrix_scaled[
  is.na(
    heatmap_matrix_scaled
  )
] <- 0


# ------------------------------------------------------------
# Cluster genes within each expression-direction group
# ------------------------------------------------------------

cluster_subset <- function(
    matrix_data,
    indices
) {
  
  if (length(indices) <= 1) {
    
    return(
      indices
    )
  }
  
  subset_matrix <- matrix_data[
    indices,
    ,
    drop = FALSE
  ]
  
  hc <- hclust(
    dist(
      subset_matrix
    ),
    method = "complete"
  )
  
  indices[
    hc$order
  ]
}


up_indices <- which(
  top24_genes$significance ==
    "Upregulated"
)

down_indices <- which(
  top24_genes$significance ==
    "Downregulated"
)

final_order <- c(
  
  cluster_subset(
    heatmap_matrix_scaled,
    up_indices
  ),
  
  cluster_subset(
    heatmap_matrix_scaled,
    down_indices
  )
)

heatmap_matrix <- heatmap_matrix_scaled[
  final_order,
  ,
  drop = FALSE
]

top24_genes <- top24_genes[
  final_order,
  ,
  drop = FALSE
]

rownames(
  heatmap_matrix
) <- top24_genes$heatmap_label

colnames(
  heatmap_matrix
) <- short_names

row_annotation <- data.frame(
  Regulation =
    top24_genes$significance
)

rownames(
  row_annotation
) <- top24_genes$heatmap_label

number_up <- sum(
  top24_genes$significance ==
    "Upregulated"
)

row_gap <- NULL

if (
  number_up > 0 &&
  number_up < nrow(
    top24_genes
  )
) {
  
  row_gap <- number_up
}


deg_heatmap <- pheatmap(
  
  heatmap_matrix,
  
  annotation_col =
    sample_annotation,
  
  annotation_row =
    row_annotation,
  
  annotation_colors = list(
    
    Condition =
      condition_colours,
    
    Regulation = c(
      "Upregulated" =
        "#D55E00",
      "Downregulated" =
        "#0072B2"
    )
  ),
  
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  
  gaps_row = row_gap,
  
  border_color = "white",
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  fontsize = 7.8,
  fontsize_row = 7.6,
  fontsize_col = 8.5,
  
  angle_col = 0,
  
  scale = "none",
  
  legend = TRUE,
  annotation_legend = TRUE,
  
  legend_breaks = c(
    -2,
    -1,
    0,
    1,
    2
  ),
  
  treeheight_col = 24,
  treeheight_row = 0,
  
  cellwidth = 27,
  cellheight = 17,
  
  silent = TRUE
)

panel_D_raw <- ggplotify::as.ggplot(
  deg_heatmap$gtable
)

panel_D <- patchwork::wrap_elements(
  full = panel_D_raw,
  clip = TRUE
)


# ------------------------------------------------------------
# 10. Export individual panels
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    panel_folder,
    "Figure1A_PCA_final.pdf"
  ),
  plot = panel_A,
  width = 5.5,
  height = 4.6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    panel_folder,
    "Figure1B_distance_final.pdf"
  ),
  plot = panel_B,
  width = 6.1,
  height = 4.8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    panel_folder,
    "Figure1C_volcano_final.pdf"
  ),
  plot = panel_C,
  width = 5.8,
  height = 5.3,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    panel_folder,
    "Figure1D_heatmap_final.pdf"
  ),
  plot = panel_D,
  width = 7.8,
  height = 7.3,
  units = "in",
  device = cairo_pdf
)


# ============================================================
# 11. Combine panels
# ============================================================

# A protected spacer column is inserted between the left
# and right panels:
#
# A | spacer | B
# C | spacer | D

layout_design <- "
ASB
CSD
"

combined_figure <- wrap_plots(
  
  A = panel_A,
  B = panel_B,
  C = panel_C,
  D = panel_D,
  S = plot_spacer(),
  
  design = layout_design,
  
  widths = c(
    0.95,
    0.08,
    1.30
  ),
  
  heights = c(
    0.84,
    1.28
  )
) +
  
  plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    
    plot.tag = element_text(
      face = "bold",
      size = 18,
      colour = "black"
    ),
    
    plot.tag.position = c(
      0.006,
      0.995
    ),
    
    plot.margin = margin(
      t = 3,
      r = 3,
      b = 3,
      l = 3
    )
  )


# ------------------------------------------------------------
# 12. Export final combined figure
# ------------------------------------------------------------

final_width <- 15.2
final_height <- 11.8

ggsave(
  filename = file.path(
    output_folder,
    "Figure_1_ABCD_final.pdf"
  ),
  plot = combined_figure,
  width = final_width,
  height = final_height,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    output_folder,
    "Figure_1_ABCD_final_600dpi.tiff"
  ),
  plot = combined_figure,
  width = final_width,
  height = final_height,
  units = "in",
  dpi = 600,
  compression = "lzw",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    output_folder,
    "Figure_1_ABCD_final_600dpi.png"
  ),
  plot = combined_figure,
  width = final_width,
  height = final_height,
  units = "in",
  dpi = 600,
  limitsize = FALSE
)


# ------------------------------------------------------------
# 13. Save supporting tables
# ------------------------------------------------------------

write_csv(
  pca_data,
  file.path(
    output_folder,
    "Figure1A_PCA_coordinates.csv"
  )
)

write_csv(
  volcano_labels,
  file.path(
    output_folder,
    "Figure1C_labels.csv"
  )
)

write_csv(
  top24_genes,
  file.path(
    output_folder,
    "Figure1D_top24_genes.csv"
  )
)


# ------------------------------------------------------------
# 14. Save session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_folder,
    "Figure1_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 15. Final report
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FINAL FIGURE COMPLETED\n")
cat("============================================\n")

cat(
  "PC1:",
  percent_variance[1],
  "%\n"
)

cat(
  "PC2:",
  percent_variance[2],
  "%\n"
)

cat(
  "Volcano labels:",
  nrow(
    volcano_labels
  ),
  "\n"
)

cat(
  "Heatmap genes:",
  nrow(
    top24_genes
  ),
  "\n"
)

cat(
  "Final figure width:",
  final_width,
  "inches\n"
)

cat(
  "Final figure height:",
  final_height,
  "inches\n"
)

cat(
  "\nFiles saved inside:\n",
  output_folder,
  "\n"
)

cat("\nGenerated files:\n")
cat("- Figure_1_ABCD_final.pdf\n")
cat("- Figure_1_ABCD_final_600dpi.tiff\n")
cat("- Figure_1_ABCD_final_600dpi.png\n")
cat("- Individual panel PDFs\n")
cat("- Supporting CSV tables\n")
cat("- Figure1_sessionInfo.txt\n")

cat("\n============================================\n")