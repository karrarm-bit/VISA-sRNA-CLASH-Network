# ============================================================
# 28_Create_Figure4_Integrated_Functional_Impact_Map.R
#
# Purpose:
# Create a manuscript-ready functional impact map in which:
# - Each unique target gene is displayed once.
# - Multiple regulatory RNAs targeting the same gene are merged.
# - Target expression direction and log2 fold change are shown.
# - Gene-specific putative biological implications are displayed.
# - CLASH-supported binding is distinguished from validated
#   regulatory direction.
#
# Input:
# D:/Bac-sRNA/16_functional_annotation/tables/
# Manuscript_functional_annotation_table.csv
#
# Outputs:
# Figure4_integrated_functional_impact_map.png
# Figure4_integrated_functional_impact_map.tiff
# Figure4_integrated_functional_impact_map.pdf
# Figure4_integrated_functional_impact_map.svg
# Figure4_integrated_functional_impact_map_data.csv
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "readr",
  "stringr",
  "ggplot2",
  "svglite"
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


# ------------------------------------------------------------
# 2. Load packages
# ------------------------------------------------------------

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(svglite)


# ------------------------------------------------------------
# 3. Input and output paths
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

input_folder <- file.path(
  project_folder,
  "16_functional_annotation",
  "tables"
)

input_file <- file.path(
  input_folder,
  "Manuscript_functional_annotation_table.csv"
)

output_folder <- file.path(
  project_folder,
  "16_functional_annotation",
  "Figure4_integrated_functional_impact_map"
)

if (!dir.exists(output_folder)) {
  dir.create(
    output_folder,
    recursive = TRUE
  )
}


output_png <- file.path(
  output_folder,
  "Figure4_integrated_functional_impact_map.png"
)

output_tiff <- file.path(
  output_folder,
  "Figure4_integrated_functional_impact_map.tiff"
)

output_pdf <- file.path(
  output_folder,
  "Figure4_integrated_functional_impact_map.pdf"
)

output_svg <- file.path(
  output_folder,
  "Figure4_integrated_functional_impact_map.svg"
)

output_csv <- file.path(
  output_folder,
  "Figure4_integrated_functional_impact_map_data.csv"
)


# ------------------------------------------------------------
# 4. Detect input file
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  
  possible_files <- list.files(
    path = input_folder,
    pattern = "Manuscript.*functional.*annotation.*\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(possible_files) == 0) {
    stop(
      paste0(
        "The functional annotation file was not found in:\n",
        input_folder
      )
    )
  }
  
  input_file <- possible_files[1]
  
  message(
    "Detected input file:\n",
    input_file
  )
}


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


rename_first_match <- function(
    data,
    new_name,
    candidates
) {
  
  detected <- candidates[
    candidates %in% names(data)
  ]
  
  if (
    length(detected) > 0 &&
    !new_name %in% names(data)
  ) {
    names(data)[
      names(data) == detected[1]
    ] <- new_name
  }
  
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


first_non_missing <- function(x) {
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA)
  }
  
  x[1]
}


# ------------------------------------------------------------
# 6. Read input data
# ------------------------------------------------------------

raw_data <- read_csv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
) |>
  clean_column_names()


# ------------------------------------------------------------
# 7. Standardise column names
# ------------------------------------------------------------

raw_data <- raw_data |>
  
  rename_first_match(
    "regulatory_rna",
    c(
      "regulatory_rna",
      "srna",
      "source_rna",
      "regulator"
    )
  ) |>
  
  rename_first_match(
    "rna_class",
    c(
      "rna_class",
      "regulatory_rna_class",
      "source_class"
    )
  ) |>
  
  rename_first_match(
    "target_gene",
    c(
      "target_gene",
      "gene_symbol",
      "gene",
      "target"
    )
  ) |>
  
  rename_first_match(
    "target_locus",
    c(
      "target_locus",
      "gene_id",
      "locus_tag"
    )
  ) |>
  
  rename_first_match(
    "target_product",
    c(
      "target_product",
      "gene_product",
      "product",
      "description"
    )
  ) |>
  
  rename_first_match(
    "functional_category",
    c(
      "functional_category",
      "functional_system",
      "functional_group",
      "category",
      "system"
    )
  ) |>
  
  rename_first_match(
    "log2_fold_change",
    c(
      "log2_fold_change",
      "log2foldchange",
      "log2fc",
      "target_log2fc"
    )
  ) |>
  
  rename_first_match(
    "adjusted_p_value",
    c(
      "adjusted_p_value",
      "deseq2_adjusted_p",
      "padj",
      "fdr"
    )
  ) |>
  
  rename_first_match(
    "confidence_class",
    c(
      "confidence_class",
      "interaction_confidence",
      "evidence_class"
    )
  )


# ------------------------------------------------------------
# 8. Validate essential columns
# ------------------------------------------------------------

required_columns <- c(
  "regulatory_rna",
  "target_gene",
  "log2_fold_change"
)

missing_columns <- setdiff(
  required_columns,
  names(raw_data)
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
        names(raw_data),
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 9. Add optional columns when absent
# ------------------------------------------------------------

if (!"rna_class" %in% names(raw_data)) {
  raw_data$rna_class <- "Regulatory RNA"
}

if (!"target_locus" %in% names(raw_data)) {
  raw_data$target_locus <- NA_character_
}

if (!"target_product" %in% names(raw_data)) {
  raw_data$target_product <- NA_character_
}

if (!"functional_category" %in% names(raw_data)) {
  raw_data$functional_category <- NA_character_
}

if (!"adjusted_p_value" %in% names(raw_data)) {
  raw_data$adjusted_p_value <- NA_real_
}

if (!"confidence_class" %in% names(raw_data)) {
  raw_data$confidence_class <- "CLASH-supported"
}


# ------------------------------------------------------------
# 10. Clean input data
# ------------------------------------------------------------

clean_data <- raw_data |>
  mutate(
    
    regulatory_rna =
      clean_text(regulatory_rna),
    
    rna_class =
      clean_text(rna_class),
    
    target_gene =
      clean_text(target_gene),
    
    target_locus =
      clean_text(target_locus),
    
    target_product =
      clean_text(target_product),
    
    functional_category =
      clean_text(functional_category),
    
    confidence_class =
      clean_text(confidence_class),
    
    log2_fold_change =
      safe_numeric(log2_fold_change),
    
    adjusted_p_value =
      safe_numeric(adjusted_p_value)
  ) |>
  
  filter(
    !is.na(regulatory_rna),
    !is.na(target_gene),
    !is.na(log2_fold_change)
  )


# ------------------------------------------------------------
# 11. Create regulatory RNA display labels
# ------------------------------------------------------------

clean_data <- clean_data |>
  mutate(
    
    regulatory_rna_display = paste0(
      regulatory_rna,
      " (",
      ifelse(
        is.na(rna_class),
        "Regulatory RNA",
        rna_class
      ),
      ")"
    )
  )


# ------------------------------------------------------------
# 12. Collapse repeated interactions to unique targets
#
# Multiple regulatory RNAs targeting the same gene are merged.
# ------------------------------------------------------------

target_data <- clean_data |>
  group_by(
    target_gene
  ) |>
  summarise(
    
    target_locus =
      first_non_missing(target_locus),
    
    target_product =
      first_non_missing(target_product),
    
    original_functional_category =
      first_non_missing(functional_category),
    
    log2_fold_change =
      first_non_missing(log2_fold_change),
    
    adjusted_p_value =
      first_non_missing(adjusted_p_value),
    
    regulatory_rnas = paste(
      unique(regulatory_rna_display),
      collapse = "\n"
    ),
    
    regulatory_rna_names = paste(
      unique(regulatory_rna),
      collapse = "; "
    ),
    
    number_of_regulatory_rnas =
      n_distinct(regulatory_rna),
    
    confidence_classes = paste(
      unique(confidence_class),
      collapse = "; "
    ),
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 13. Correct functional categories using gene-specific rules
# ------------------------------------------------------------

target_data <- target_data |>
  mutate(
    
    functional_system = case_when(
      
      str_to_lower(target_gene) == "guab" ~
        "Purine nucleotide biosynthesis",
      
      str_to_lower(target_gene) == "nrdf" ~
        "DNA-precursor synthesis",
      
      str_to_lower(target_gene) == "moda" ~
        "Molybdate transport",
      
      str_to_lower(target_gene) == "qoxb" ~
        "Respiratory electron transport",
      
      str_to_lower(target_gene) == "spa" ~
        "Surface virulence and host interaction",
      
      str_to_lower(target_gene) == "dhak" ~
        "Carbon metabolism",
      
      str_to_lower(target_gene) == "gltb" ~
        "Nitrogen and amino-acid metabolism",
      
      !is.na(original_functional_category) ~
        original_functional_category,
      
      TRUE ~
        "Functional system not assigned"
    )
  )


# ------------------------------------------------------------
# 14. Classify observed target response
# ------------------------------------------------------------

target_data <- target_data |>
  mutate(
    
    target_direction = case_when(
      
      log2_fold_change < 0 ~
        "Downregulated",
      
      log2_fold_change > 0 ~
        "Upregulated",
      
      TRUE ~
        "No change"
    ),
    
    direction_symbol = case_when(
      
      target_direction == "Downregulated" ~
        "\u2193",
      
      target_direction == "Upregulated" ~
        "\u2191",
      
      TRUE ~
        "\u2192"
    ),
    
    regulatory_interpretation = case_when(
      
      target_direction == "Downregulated" ~
        "Directionally compatible with putative repression",
      
      target_direction == "Upregulated" ~
        "Mechanistically ambiguous",
      
      TRUE ~
        "No directional transcriptomic response"
    )
  )


# ------------------------------------------------------------
# 15. Assign gene-specific biological implications
#
# Statements are deliberately phrased as putative implications.
# ------------------------------------------------------------

target_data <- target_data |>
  mutate(
    
    biological_implication = case_when(
      
      str_to_lower(target_gene) == "guab" &
        target_direction == "Downregulated" ~
        "Potential reduction in purine nucleotide biosynthesis",
      
      str_to_lower(target_gene) == "guab" &
        target_direction == "Upregulated" ~
        "Potential enhancement of purine nucleotide biosynthesis",
      
      str_to_lower(target_gene) == "nrdf" &
        target_direction == "Downregulated" ~
        "Potential reduction in DNA-precursor synthesis",
      
      str_to_lower(target_gene) == "nrdf" &
        target_direction == "Upregulated" ~
        "Potential enhancement of DNA-precursor synthesis",
      
      str_to_lower(target_gene) == "moda" &
        target_direction == "Downregulated" ~
        "Potential reduction in molybdate uptake",
      
      str_to_lower(target_gene) == "moda" &
        target_direction == "Upregulated" ~
        "Potential enhancement of molybdate uptake",
      
      str_to_lower(target_gene) == "qoxb" &
        target_direction == "Downregulated" ~
        "Potential reduction in respiratory electron transport",
      
      str_to_lower(target_gene) == "qoxb" &
        target_direction == "Upregulated" ~
        "Potential enhancement of respiratory adaptation",
      
      str_to_lower(target_gene) == "spa" &
        target_direction == "Downregulated" ~
        "Potential reduction in surface-associated virulence",
      
      str_to_lower(target_gene) == "spa" &
        target_direction == "Upregulated" ~
        "Potential enhancement of surface-associated virulence",
      
      str_to_lower(target_gene) == "dhak" &
        target_direction == "Downregulated" ~
        "Potential reduction in carbon metabolic activity",
      
      str_to_lower(target_gene) == "dhak" &
        target_direction == "Upregulated" ~
        "Potential carbon-metabolic adaptation",
      
      str_to_lower(target_gene) == "gltb" &
        target_direction == "Downregulated" ~
        "Potential reduction in nitrogen and amino-acid metabolism",
      
      str_to_lower(target_gene) == "gltb" &
        target_direction == "Upregulated" ~
        "Potential nitrogen and amino-acid metabolic adaptation",
      
      TRUE ~
        "Functional consequence remains uncertain"
    )
  )


# ------------------------------------------------------------
# 16. Create figure labels
# ------------------------------------------------------------

target_data <- target_data |>
  mutate(
    
    regulatory_label = regulatory_rnas,
    
    target_label = paste0(
      target_gene,
      "  ",
      direction_symbol,
      "\nlog2FC = ",
      sprintf(
        "%.2f",
        log2_fold_change
      )
    ),
    
    implication_label = paste0(
      functional_system,
      "\n",
      str_wrap(
        biological_implication,
        width = 42
      )
    ),
    
    interpretation_short = case_when(
      
      target_direction == "Downregulated" ~
        "Compatible with putative repression",
      
      target_direction == "Upregulated" ~
        "Direction remains mechanistically ambiguous",
      
      TRUE ~
        "No directional response"
    )
  )


# ------------------------------------------------------------
# 17. Define target ordering
# ------------------------------------------------------------

preferred_order <- c(
  "spa",
  "modA",
  "qoxB",
  "nrdF",
  "guaB",
  "dhaK",
  "gltB"
)

target_data <- target_data |>
  mutate(
    
    ordering_index = match(
      target_gene,
      preferred_order
    ),
    
    ordering_index = ifelse(
      is.na(ordering_index),
      999,
      ordering_index
    )
  ) |>
  
  arrange(
    ordering_index,
    target_gene
  ) |>
  
  mutate(
    row_position =
      rev(seq_len(n()))
  )


# ------------------------------------------------------------
# 18. Export processed data
# ------------------------------------------------------------

write_csv(
  target_data |>
    select(
      target_gene,
      target_locus,
      target_product,
      regulatory_rna_names,
      number_of_regulatory_rnas,
      confidence_classes,
      log2_fold_change,
      adjusted_p_value,
      target_direction,
      regulatory_interpretation,
      functional_system,
      biological_implication
    ),
  output_csv,
  na = ""
)


# ------------------------------------------------------------
# 19. Set plot positions
# ------------------------------------------------------------

x_regulatory <- 0.0
x_target <- 2.1
x_implication <- 4.65

target_data <- target_data |>
  mutate(
    x_regulatory = x_regulatory,
    x_target = x_target,
    x_implication = x_implication
  )


# ------------------------------------------------------------
# 20. Create Figure 4
# ------------------------------------------------------------

figure4 <- ggplot(
  target_data
) +
  
  # CLASH-supported regulatory RNA-to-target interaction
  geom_segment(
    aes(
      x = x_regulatory + 0.60,
      xend = x_target - 0.47,
      y = row_position,
      yend = row_position
    ),
    linewidth = 0.8,
    linetype = "dashed",
    colour = "#666666",
    arrow = arrow(
      length = grid::unit(
        0.12,
        "inches"
      ),
      type = "closed"
    )
  ) +
  
  # Target-to-functional implication link
  geom_segment(
    aes(
      x = x_target + 0.49,
      xend = x_implication - 0.72,
      y = row_position,
      yend = row_position
    ),
    linewidth = 0.8,
    colour = "#8A8A8A",
    arrow = arrow(
      length = grid::unit(
        0.12,
        "inches"
      ),
      type = "closed"
    )
  ) +
  
  # Regulatory RNA boxes
  geom_label(
    aes(
      x = x_regulatory,
      y = row_position,
      label = regulatory_label
    ),
    fill = "#EAF2F8",
    colour = "#1F4E78",
    fontface = "bold",
    size = 3.15,
    lineheight = 0.93,
    label.size = 0.4,
    label.padding = grid::unit(
      0.25,
      "lines"
    )
  ) +
  
  # Target response boxes
  geom_label(
    aes(
      x = x_target,
      y = row_position,
      label = target_label,
      fill = target_direction
    ),
    colour = "black",
    fontface = "bold",
    size = 3.65,
    lineheight = 0.95,
    label.size = 0.4,
    label.padding = grid::unit(
      0.28,
      "lines"
    )
  ) +
  
  # Biological implication boxes
  geom_label(
    aes(
      x = x_implication,
      y = row_position,
      label = implication_label,
      fill = functional_system
    ),
    colour = "#202020",
    size = 3.15,
    lineheight = 0.93,
    label.size = 0.4,
    label.padding = grid::unit(
      0.25,
      "lines"
    )
  ) +
  
  # Target-response colours
  scale_fill_manual(
    aesthetics = "fill",
    values = c(
      
      "Downregulated" = "#F4B7B7",
      "Upregulated" = "#B8DEC7",
      "No change" = "#D9D9D9",
      
      "Purine nucleotide biosynthesis" = "#F8E6A0",
      "DNA-precursor synthesis" = "#F3D17A",
      "Molybdate transport" = "#B9D7F0",
      "Respiratory electron transport" = "#D5C5F2",
      "Surface virulence and host interaction" = "#F3C1D9",
      "Carbon metabolism" = "#C7E8D3",
      "Nitrogen and amino-acid metabolism" = "#B8E0D6",
      "Functional system not assigned" = "#E6E6E6"
    ),
    guide = "none"
  ) +
  
  # Column headings
  annotate(
    "text",
    x = x_regulatory,
    y = max(target_data$row_position) + 1.05,
    label = "Regulatory RNA(s)",
    fontface = "bold",
    size = 4.8
  ) +
  
  annotate(
    "text",
    x = x_target,
    y = max(target_data$row_position) + 1.05,
    label = "Observed target response",
    fontface = "bold",
    size = 4.8
  ) +
  
  annotate(
    "text",
    x = x_implication,
    y = max(target_data$row_position) + 1.05,
    label = "Putative biological implication",
    fontface = "bold",
    size = 4.8
  ) +
  
  # CLASH annotation
  annotate(
    "text",
    x = (
      x_regulatory +
        x_target
    ) / 2,
    y = max(target_data$row_position) + 0.53,
    label = "CLASH-supported interaction",
    fontface = "italic",
    colour = "#555555",
    size = 3.4
  ) +
  
  # Repression-compatible label
  annotate(
    "label",
    x = 2.1,
    y = 0.13,
    label = paste0(
      "Downregulated targets: directionally compatible with ",
      "putative repression"
    ),
    fill = "#F9E1E1",
    colour = "#7F1D1D",
    size = 3.3,
    fontface = "italic",
    label.size = 0.3
  ) +
  
  # Ambiguous upregulation label
  annotate(
    "label",
    x = 4.45,
    y = 0.13,
    label = paste0(
      "Upregulated targets: regulatory direction remains ",
      "mechanistically ambiguous"
    ),
    fill = "#E0F0E6",
    colour = "#245C3B",
    size = 3.3,
    fontface = "italic",
    label.size = 0.3
  ) +
  
  # Core caution statement
  annotate(
    "text",
    x = 2.65,
    y = -0.38,
    label = paste0(
      "CLASH-supported binding does not establish the direction ",
      "of post-transcriptional regulation."
    ),
    fontface = "italic",
    colour = "#555555",
    size = 3.55
  ) +
  
  scale_x_continuous(
    limits = c(
      -1.05,
      5.85
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  scale_y_continuous(
    limits = c(
      -0.70,
      max(target_data$row_position) + 1.40
    ),
    breaks = NULL,
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_void(
    base_size = 12
  ) +
  
  theme(
    
    text = element_text(
      family = "Arial"
    ),
    
    plot.margin = margin(
      28,
      40,
      42,
      40
    )
  )


# ------------------------------------------------------------
# 21. Display Figure 4
# ------------------------------------------------------------

print(figure4)


# ------------------------------------------------------------
# 22. Save high-resolution PNG
# ------------------------------------------------------------

ggsave(
  filename = output_png,
  plot = figure4,
  width = 15,
  height = 9.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 23. Save journal-quality TIFF
# ------------------------------------------------------------

ggsave(
  filename = output_tiff,
  plot = figure4,
  width = 15,
  height = 9.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 24. Save vector PDF
# ------------------------------------------------------------

ggsave(
  filename = output_pdf,
  plot = figure4,
  width = 15,
  height = 9.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# ------------------------------------------------------------
# 25. Save editable SVG
# ------------------------------------------------------------

ggsave(
  filename = output_svg,
  plot = figure4,
  width = 15,
  height = 9.5,
  units = "in",
  device = svglite,
  bg = "white"
)


# ------------------------------------------------------------
# 26. Final report
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "Integrated Figure 4 created successfully\n",
  "============================================\n\n",
  
  "Unique target genes displayed: ",
  nrow(target_data),
  "\n",
  
  "Regulatory RNA-target interactions represented: ",
  sum(target_data$number_of_regulatory_rnas),
  "\n\n",
  
  "PNG:\n",
  normalizePath(
    output_png,
    winslash = "/"
  ),
  "\n\n",
  
  "TIFF:\n",
  normalizePath(
    output_tiff,
    winslash = "/"
  ),
  "\n\n",
  
  "PDF:\n",
  normalizePath(
    output_pdf,
    winslash = "/"
  ),
  "\n\n",
  
  "Editable SVG:\n",
  normalizePath(
    output_svg,
    winslash = "/"
  ),
  "\n\n",
  
  "Processed data:\n",
  normalizePath(
    output_csv,
    winslash = "/"
  ),
  "\n\n"
)

print(
  target_data |>
    select(
      regulatory_rna_names,
      target_gene,
      log2_fold_change,
      target_direction,
      functional_system,
      biological_implication,
      regulatory_interpretation
    )
)