# ============================================================
# BacRegRNA Project
# Script 09 FINAL FIXED
# GSE254530 RNA-seq differential-expression analysis
#
# Comparison:
#   Vancomycin-treated versus untreated
#
# Samples:
#   3 untreated
#   3 vancomycin-treated
# ============================================================


# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

input_folder <- file.path(
  project_folder,
  "07_processed_data",
  "downloaded_files",
  "GSE254530"
)

design_file <- file.path(
  project_folder,
  "06_curated_design",
  "BacRegRNA_final_curated_sample_design.csv"
)

output_folder <- file.path(
  project_folder,
  "09_RNAseq_analysis"
)

figure_folder <- file.path(
  output_folder,
  "figures"
)

table_folder <- file.path(
  output_folder,
  "tables"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
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

setwd(project_folder)


# ------------------------------------------------------------
# 2. Install required CRAN packages
# ------------------------------------------------------------

cran_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "purrr",
  "tibble",
  "tidyr",
  "ggplot2",
  "pheatmap",
  "ggrepel"
)

missing_cran <- cran_packages[
  !cran_packages %in% rownames(
    installed.packages()
  )
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}


# ------------------------------------------------------------
# 3. Install Bioconductor packages
# ------------------------------------------------------------

if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "DESeq2",
  "EnhancedVolcano"
)

missing_bioc <- bioc_packages[
  !vapply(
    bioc_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_bioc) > 0) {
  
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

# apeglm is optional
if (!requireNamespace(
  "apeglm",
  quietly = TRUE
)) {
  
  try(
    BiocManager::install(
      "apeglm",
      ask = FALSE,
      update = FALSE
    ),
    silent = TRUE
  )
}


# ------------------------------------------------------------
# 4. Load packages
# ------------------------------------------------------------

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(DESeq2)
library(EnhancedVolcano)


# ------------------------------------------------------------
# 5. Verify input folders and files
# ------------------------------------------------------------

if (!dir.exists(input_folder)) {
  
  stop(
    "Input folder was not found:\n",
    input_folder
  )
}

if (!file.exists(design_file)) {
  
  stop(
    "Curated sample-design file was not found:\n",
    design_file
  )
}


# ------------------------------------------------------------
# 6. Locate six RNA-seq count files
# ------------------------------------------------------------

count_files <- list.files(
  input_folder,
  pattern = "^GSM[0-9]+_.*RNA\\.txt\\.gz$",
  full.names = TRUE
)

count_files <- count_files[
  !str_detect(
    basename(count_files),
    regex(
      "DESeq2|DEseq2",
      ignore_case = TRUE
    )
  )
]

if (length(count_files) != 6) {
  
  stop(
    "Expected 6 RNA count files, but found ",
    length(count_files),
    ".\n\nDetected files:\n",
    paste(
      basename(count_files),
      collapse = "\n"
    )
  )
}

cat("\nRNA count files detected:\n")

print(
  basename(count_files)
)


# ------------------------------------------------------------
# 7. Function to read one count file
# ------------------------------------------------------------

read_count_file <- function(path) {
  
  sample_name <- basename(path) |>
    str_remove("\\.txt\\.gz$")
  
  data <- read_tsv(
    file = path,
    col_names = c(
      "gene_id",
      sample_name
    ),
    col_types = cols(
      gene_id = col_character(),
      .default = col_double()
    ),
    show_col_types = FALSE,
    progress = FALSE,
    trim_ws = TRUE
  )
  
  data <- data |>
    filter(
      !is.na(gene_id),
      gene_id != ""
    )
  
  if (ncol(data) != 2) {
    
    stop(
      "Unexpected number of columns in:\n",
      path
    )
  }
  
  if (anyDuplicated(data$gene_id) > 0) {
    
    data <- data |>
      group_by(gene_id) |>
      summarise(
        across(
          everything(),
          ~ sum(.x, na.rm = TRUE)
        ),
        .groups = "drop"
      )
  }
  
  data
}


# ------------------------------------------------------------
# 8. Read all count files
# ------------------------------------------------------------

count_tables <- purrr::map(
  count_files,
  read_count_file
)


# ------------------------------------------------------------
# 9. Build raw count table
# Explicit purrr::reduce avoids Bioconductor conflict
# ------------------------------------------------------------

count_data <- purrr::reduce(
  count_tables,
  dplyr::full_join,
  by = "gene_id"
) |>
  dplyr::mutate(
    dplyr::across(
      -gene_id,
      ~ tidyr::replace_na(.x, 0)
    )
  )


# ------------------------------------------------------------
# 10. Convert to integer count matrix
# ------------------------------------------------------------

count_matrix <- count_data |>
  tibble::column_to_rownames(
    "gene_id"
  ) |>
  as.matrix()

if (any(count_matrix < 0, na.rm = TRUE)) {
  
  stop(
    "Negative values were detected in the count matrix."
  )
}

count_matrix <- round(
  count_matrix
)

storage.mode(count_matrix) <- "integer"

cat(
  "\nCount matrix dimensions:",
  nrow(count_matrix),
  "genes x",
  ncol(count_matrix),
  "samples\n"
)

cat("\nCount-matrix columns:\n")

print(
  colnames(count_matrix)
)


# ------------------------------------------------------------
# 11. Read final curated sample design
# ------------------------------------------------------------

full_design <- read_csv(
  design_file,
  show_col_types = FALSE
)

sample_design <- full_design |>
  filter(
    dataset == "GSE254530",
    assay_curated == "RNA-seq",
    vancomycin_status %in% c(
      "Untreated",
      "Treated"
    )
  ) |>
  mutate(
    replicate_curated = suppressWarnings(
      as.numeric(replicate_curated)
    )
  )

if (nrow(sample_design) != 6) {
  
  stop(
    "Expected 6 curated RNA-seq samples, but found ",
    nrow(sample_design)
  )
}


# ------------------------------------------------------------
# 12. Helper for matching one count column
# ------------------------------------------------------------

find_one_column <- function(
    pattern,
    available_columns
) {
  
  matches <- available_columns[
    str_detect(
      available_columns,
      regex(
        pattern,
        ignore_case = TRUE
      )
    )
  ]
  
  if (length(matches) != 1) {
    return(NA_character_)
  }
  
  matches[1]
}


# ------------------------------------------------------------
# 13. Match metadata samples to count columns
# ------------------------------------------------------------

sample_design <- sample_design |>
  mutate(
    
    file_sample_name = case_when(
      
      vancomycin_status == "Untreated" &
        replicate_curated == 1 ~
        find_one_column(
          "Cont1_30mins_RNA$",
          colnames(count_matrix)
        ),
      
      vancomycin_status == "Untreated" &
        replicate_curated == 2 ~
        find_one_column(
          "Cont2_30mins_RNA$",
          colnames(count_matrix)
        ),
      
      vancomycin_status == "Untreated" &
        replicate_curated == 3 ~
        find_one_column(
          "Cont3_30mins_RNA$",
          colnames(count_matrix)
        ),
      
      vancomycin_status == "Treated" &
        replicate_curated == 1 ~
        find_one_column(
          "Van1_30mins_RNA$",
          colnames(count_matrix)
        ),
      
      vancomycin_status == "Treated" &
        replicate_curated == 2 ~
        find_one_column(
          "Van2_30mins_RNA$",
          colnames(count_matrix)
        ),
      
      vancomycin_status == "Treated" &
        replicate_curated == 3 ~
        find_one_column(
          "Van3_30mins_RNA$",
          colnames(count_matrix)
        ),
      
      TRUE ~ NA_character_
    )
  )

if (any(is.na(
  sample_design$file_sample_name
))) {
  
  cat(
    "\nSamples that could not be matched:\n"
  )
  
  print(
    sample_design |>
      filter(
        is.na(file_sample_name)
      ) |>
      select(
        sample_title,
        vancomycin_status,
        replicate_curated
      )
  )
  
  stop(
    "Some samples could not be matched to count files."
  )
}

if (anyDuplicated(
  sample_design$file_sample_name
) > 0) {
  
  stop(
    "Duplicate count-file assignments were detected."
  )
}


# ------------------------------------------------------------
# 14. Create condition factor and reorder samples
# ------------------------------------------------------------

sample_design <- sample_design |>
  mutate(
    condition = factor(
      vancomycin_status,
      levels = c(
        "Untreated",
        "Treated"
      )
    )
  ) |>
  arrange(
    condition,
    replicate_curated
  )

count_matrix <- count_matrix[
  ,
  sample_design$file_sample_name,
  drop = FALSE
]

rownames(sample_design) <-
  sample_design$file_sample_name

if (!identical(
  colnames(count_matrix),
  rownames(sample_design)
)) {
  
  stop(
    "Count-matrix columns do not match metadata rows."
  )
}

cat("\nFinal sample order:\n")

print(
  sample_design |>
    select(
      file_sample_name,
      condition,
      replicate_curated
    )
)


# ------------------------------------------------------------
# 15. Save raw count matrix and sample design
# ------------------------------------------------------------

write_csv(
  count_data,
  file.path(
    table_folder,
    "GSE254530_raw_count_matrix.csv"
  )
)

sample_design_export <- sample_design |>
  tibble::rownames_to_column(
    "matrix_column"
  ) |>
  select(
    matrix_column,
    dataset,
    sample_id,
    gsm_accession,
    sample_title,
    file_sample_name,
    condition,
    replicate_curated,
    concentration_curated,
    exposure_time_curated
  )

write_csv(
  sample_design_export,
  file.path(
    table_folder,
    "GSE254530_sample_design.csv"
  )
)


# ------------------------------------------------------------
# 16. Filter low-count genes
# Rule: >=10 reads in at least 3 samples
# ------------------------------------------------------------

keep_genes <- rowSums(
  count_matrix >= 10
) >= 3

filtered_count_matrix <- count_matrix[
  keep_genes,
  ,
  drop = FALSE
]

if (nrow(filtered_count_matrix) == 0) {
  
  stop(
    "No genes remained after filtering."
  )
}

filter_summary <- tibble(
  total_genes = nrow(count_matrix),
  retained_genes =
    nrow(filtered_count_matrix),
  removed_genes =
    nrow(count_matrix) -
    nrow(filtered_count_matrix),
  retained_percent = round(
    100 *
      nrow(filtered_count_matrix) /
      nrow(count_matrix),
    2
  ),
  filter_rule =
    "Count >=10 in at least 3 of 6 samples"
)

write_csv(
  filter_summary,
  file.path(
    table_folder,
    "GSE254530_filter_summary.csv"
  )
)

cat(
  "\nGenes retained after filtering:",
  nrow(filtered_count_matrix),
  "of",
  nrow(count_matrix),
  "\n"
)


# ------------------------------------------------------------
# 17. Create DESeq2 object
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = filtered_count_matrix,
  colData = sample_design,
  design = ~ condition
)

dds$condition <- relevel(
  dds$condition,
  ref = "Untreated"
)


# ------------------------------------------------------------
# 18. Run DESeq2
# ------------------------------------------------------------

dds <- DESeq(
  dds,
  quiet = FALSE
)

cat("\nDESeq2 coefficient names:\n")

print(
  resultsNames(dds)
)


# ------------------------------------------------------------
# 19. Standard DESeq2 results
# P-values and padj are taken from here
# ------------------------------------------------------------

standard_results <- results(
  dds,
  contrast = c(
    "condition",
    "Treated",
    "Untreated"
  ),
  alpha = 0.05
)


# ------------------------------------------------------------
# 20. LFC shrinkage
# Only effect-size estimates are taken from this object
# ------------------------------------------------------------

coefficient_names <- resultsNames(dds)

coefficient_name <- coefficient_names[
  str_detect(
    coefficient_names,
    "condition_Treated_vs_Untreated"
  )
][1]

if (
  requireNamespace(
    "apeglm",
    quietly = TRUE
  ) &&
  !is.na(coefficient_name)
) {
  
  shrunk_results <- tryCatch(
    
    lfcShrink(
      dds,
      coef = coefficient_name,
      type = "apeglm"
    ),
    
    error = function(e) {
      
      message(
        "apeglm shrinkage failed. ",
        "Standard log2FC values will be used. Reason: ",
        conditionMessage(e)
      )
      
      standard_results
    }
  )
  
} else {
  
  message(
    "apeglm was unavailable. ",
    "Standard DESeq2 log2FC values will be used."
  )
  
  shrunk_results <- standard_results
}


# ------------------------------------------------------------
# 21. Convert standard results safely
# This section fixes: object 'baseMean' not found
# ------------------------------------------------------------

standard_df <- as.data.frame(
  standard_results
)

standard_df$gene_id <- rownames(
  standard_df
)

required_standard_columns <- c(
  "gene_id",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)

missing_standard_columns <- setdiff(
  required_standard_columns,
  names(standard_df)
)

if (length(
  missing_standard_columns
) > 0) {
  
  cat(
    "\nColumns actually present in standard results:\n"
  )
  
  print(
    names(standard_df)
  )
  
  stop(
    "Missing expected columns in standard DESeq2 results: ",
    paste(
      missing_standard_columns,
      collapse = ", "
    )
  )
}

standard_table <- tibble::as_tibble(
  standard_df
) |>
  dplyr::transmute(
    
    gene_id =
      as.character(.data$gene_id),
    
    base_mean =
      as.numeric(.data$baseMean),
    
    unshrunk_log2_fold_change =
      as.numeric(.data$log2FoldChange),
    
    standard_error =
      as.numeric(.data$lfcSE),
    
    test_statistic =
      as.numeric(.data$stat),
    
    p_value =
      as.numeric(.data$pvalue),
    
    adjusted_p_value =
      as.numeric(.data$padj)
  )


# ------------------------------------------------------------
# 22. Convert shrunken results safely
# ------------------------------------------------------------

shrunk_df <- as.data.frame(
  shrunk_results
)

shrunk_df$gene_id <- rownames(
  shrunk_df
)

if (!"log2FoldChange" %in%
    names(shrunk_df)) {
  
  stop(
    "Shrunken results do not contain log2FoldChange."
  )
}

if ("lfcSE" %in% names(shrunk_df)) {
  
  shrunk_table <- tibble::as_tibble(
    shrunk_df
  ) |>
    dplyr::transmute(
      
      gene_id =
        as.character(.data$gene_id),
      
      log2_fold_change =
        as.numeric(.data$log2FoldChange),
      
      shrunken_standard_error =
        as.numeric(.data$lfcSE)
    )
  
} else {
  
  shrunk_table <- tibble::as_tibble(
    shrunk_df
  ) |>
    dplyr::transmute(
      
      gene_id =
        as.character(.data$gene_id),
      
      log2_fold_change =
        as.numeric(.data$log2FoldChange),
      
      shrunken_standard_error =
        NA_real_
    )
}


# ------------------------------------------------------------
# 23. Combine statistics and shrunken effect sizes
# ------------------------------------------------------------

results_table <- standard_table |>
  dplyr::left_join(
    shrunk_table,
    by = "gene_id"
  ) |>
  dplyr::mutate(
    
    log2_fold_change =
      dplyr::coalesce(
        log2_fold_change,
        unshrunk_log2_fold_change
      ),
    
    significance =
      dplyr::case_when(
        
        !is.na(adjusted_p_value) &
          adjusted_p_value < 0.05 &
          log2_fold_change >= 1 ~
          "Upregulated",
        
        !is.na(adjusted_p_value) &
          adjusted_p_value < 0.05 &
          log2_fold_change <= -1 ~
          "Downregulated",
        
        TRUE ~
          "Not significant"
      )
  ) |>
  dplyr::arrange(
    is.na(adjusted_p_value),
    adjusted_p_value,
    dplyr::desc(
      abs(log2_fold_change)
    )
  )


# ------------------------------------------------------------
# 24. Save complete and significant DE results
# ------------------------------------------------------------

write_csv(
  results_table,
  file.path(
    table_folder,
    "GSE254530_DESeq2_results_all_genes.csv"
  )
)

significant_results <- results_table |>
  filter(
    significance %in% c(
      "Upregulated",
      "Downregulated"
    )
  )

write_csv(
  significant_results,
  file.path(
    table_folder,
    "GSE254530_DESeq2_significant_genes.csv"
  )
)

de_summary <- results_table |>
  count(
    significance,
    name = "number_of_genes"
  ) |>
  tidyr::complete(
    significance = c(
      "Upregulated",
      "Downregulated",
      "Not significant"
    ),
    fill = list(
      number_of_genes = 0
    )
  )

write_csv(
  de_summary,
  file.path(
    table_folder,
    "GSE254530_DE_summary.csv"
  )
)


# ------------------------------------------------------------
# 25. Save normalized counts
# ------------------------------------------------------------

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

normalized_counts_table <-
  as.data.frame(
    normalized_counts
  ) |>
  tibble::rownames_to_column(
    "gene_id"
  ) |>
  as_tibble()

write_csv(
  normalized_counts_table,
  file.path(
    table_folder,
    "GSE254530_normalized_counts.csv"
  )
)


# ------------------------------------------------------------
# 26. Variance-stabilizing transformation
# ------------------------------------------------------------

vsd <- vst(
  dds,
  blind = TRUE
)

vst_matrix <- assay(
  vsd
)

write_csv(
  as.data.frame(vst_matrix) |>
    tibble::rownames_to_column(
      "gene_id"
    ),
  file.path(
    table_folder,
    "GSE254530_VST_expression_matrix.csv"
  )
)


# ------------------------------------------------------------
# 27. PCA
# ------------------------------------------------------------

pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percent_variance <- round(
  100 *
    attr(
      pca_data,
      "percentVar"
    ),
  1
)

pca_data_export <- pca_data |>
  tibble::rownames_to_column(
    "sample_row"
  )

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    shape = condition,
    label = name
  )
) +
  geom_point(
    size = 4
  ) +
  ggrepel::geom_text_repel(
    max.overlaps = Inf
  ) +
  labs(
    title = "GSE254530 RNA-seq PCA",
    subtitle =
      "Vancomycin-treated versus untreated samples",
    x = paste0(
      "PC1: ",
      percent_variance[1],
      "% variance"
    ),
    y = paste0(
      "PC2: ",
      percent_variance[2],
      "% variance"
    ),
    shape = "Condition"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    figure_folder,
    "GSE254530_PCA.png"
  ),
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

write_csv(
  pca_data_export,
  file.path(
    table_folder,
    "GSE254530_PCA_coordinates.csv"
  )
)


# ------------------------------------------------------------
# 28. Sample-distance heatmap
# ------------------------------------------------------------

sample_distance_matrix <- as.matrix(
  dist(
    t(vst_matrix)
  )
)

annotation_data <- data.frame(
  Condition =
    sample_design$condition
)

rownames(annotation_data) <-
  rownames(sample_design)

png(
  filename = file.path(
    figure_folder,
    "GSE254530_sample_distance_heatmap.png"
  ),
  width = 1800,
  height = 1600,
  res = 250
)

pheatmap(
  sample_distance_matrix,
  annotation_col = annotation_data,
  annotation_row = annotation_data,
  main =
    "GSE254530 sample-to-sample distance"
)

dev.off()


# ------------------------------------------------------------
# 29. Volcano plot
# ------------------------------------------------------------

volcano_data <- results_table |>
  mutate(
    
    volcano_adjusted_p = case_when(
      
      is.na(adjusted_p_value) ~ 1,
      
      adjusted_p_value <= 0 ~
        .Machine$double.xmin,
      
      TRUE ~ adjusted_p_value
    )
  )

label_genes <- significant_results |>
  filter(
    !is.na(adjusted_p_value)
  ) |>
  arrange(
    adjusted_p_value
  ) |>
  slice_head(
    n = min(
      20,
      n()
    )
  ) |>
  pull(
    gene_id
  )

volcano_plot <- EnhancedVolcano(
  volcano_data,
  lab = volcano_data$gene_id,
  x = "log2_fold_change",
  y = "volcano_adjusted_p",
  pCutoff = 0.05,
  FCcutoff = 1,
  title =
    "Vancomycin-treated versus untreated",
  subtitle =
    "GSE254530 RNA-seq",
  caption =
    "Adjusted P < 0.05 and |log2FC| >= 1",
  pointSize = 2,
  labSize = 3,
  selectLab = label_genes
)

ggsave(
  filename = file.path(
    figure_folder,
    "GSE254530_volcano_plot.png"
  ),
  plot = volcano_plot,
  width = 9,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 30. Heatmap of top 50 significant genes
# ------------------------------------------------------------

top_genes <- significant_results |>
  filter(
    !is.na(adjusted_p_value)
  ) |>
  arrange(
    adjusted_p_value
  ) |>
  slice_head(
    n = min(
      50,
      n()
    )
  ) |>
  pull(
    gene_id
  )

if (length(top_genes) >= 2) {
  
  top_heatmap_matrix <- vst_matrix[
    top_genes,
    ,
    drop = FALSE
  ]
  
  top_heatmap_matrix <- t(
    scale(
      t(top_heatmap_matrix)
    )
  )
  
  top_heatmap_matrix[
    is.na(top_heatmap_matrix)
  ] <- 0
  
  png(
    filename = file.path(
      figure_folder,
      "GSE254530_top50_DEG_heatmap.png"
    ),
    width = 1800,
    height = 2200,
    res = 250
  )
  
  pheatmap(
    top_heatmap_matrix,
    annotation_col = annotation_data,
    show_rownames = TRUE,
    fontsize_row = 6,
    main =
      "Top vancomycin-responsive genes"
  )
  
  dev.off()
}


# ------------------------------------------------------------
# 31. Locate authors' DESeq2 file
# ------------------------------------------------------------

author_file <- list.files(
  input_folder,
  pattern =
    "GSE254530_RNAseq_VAN30_DEseq2\\.txt\\.gz$",
  full.names = TRUE
)

author_comparison <- tibble()
comparison_summary <- tibble()


# ------------------------------------------------------------
# 32. Compare our results with authors' results
# ------------------------------------------------------------

if (length(author_file) == 1) {
  
  author_results <- read_tsv(
    author_file,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  names(author_results) <- names(
    author_results
  ) |>
    str_to_lower() |>
    str_replace_all(
      "[^a-z0-9]+",
      "_"
    ) |>
    str_replace_all(
      "^_|_$",
      ""
    )
  
  cat(
    "\nPublished result columns:\n"
  )
  
  print(
    names(author_results)
  )
  
  author_gene_col <- names(
    author_results
  )[
    str_detect(
      names(author_results),
      "^gene_id$"
    )
  ][1]
  
  author_lfc_col <- names(
    author_results
  )[
    str_detect(
      names(author_results),
      "^log2_fc$|log2.*fc"
    )
  ][1]
  
  author_padj_col <- names(
    author_results
  )[
    str_detect(
      names(author_results),
      "^p_adj$|^padj$|adjust"
    )
  ][1]
  
  valid_author_columns <- all(
    !is.na(
      c(
        author_gene_col,
        author_lfc_col,
        author_padj_col
      )
    )
  )
  
  if (valid_author_columns) {
    
    author_clean <- author_results |>
      transmute(
        
        gene_id =
          as.character(
            .data[[author_gene_col]]
          ),
        
        author_log2_fold_change =
          suppressWarnings(
            as.numeric(
              .data[[author_lfc_col]]
            )
          ),
        
        author_adjusted_p_value =
          suppressWarnings(
            as.numeric(
              .data[[author_padj_col]]
            )
          )
      )
    
    author_comparison <-
      results_table |>
      select(
        
        gene_id,
        
        our_log2_fold_change =
          log2_fold_change,
        
        our_adjusted_p_value =
          adjusted_p_value
      ) |>
      inner_join(
        author_clean,
        by = "gene_id"
      ) |>
      mutate(
        
        same_direction =
          case_when(
            
            sign(
              our_log2_fold_change
            ) ==
              sign(
                author_log2_fold_change
              ) ~ "Yes",
            
            TRUE ~ "No"
          ),
        
        both_significant =
          case_when(
            
            !is.na(
              our_adjusted_p_value
            ) &
              
              !is.na(
                author_adjusted_p_value
              ) &
              
              our_adjusted_p_value < 0.05 &
              
              author_adjusted_p_value < 0.05 ~
              "Yes",
            
            TRUE ~ "No"
          )
      )
    
    write_csv(
      author_comparison,
      file.path(
        table_folder,
        "GSE254530_comparison_with_author_results.csv"
      )
    )
    
    correlation_value <- cor(
      author_comparison$
        our_log2_fold_change,
      author_comparison$
        author_log2_fold_change,
      use = "complete.obs",
      method = "pearson"
    )
    
    comparison_summary <- tibble(
      
      genes_compared =
        nrow(author_comparison),
      
      pearson_log2FC_correlation =
        correlation_value,
      
      same_direction_percent =
        round(
          mean(
            author_comparison$
              same_direction == "Yes",
            na.rm = TRUE
          ) * 100,
          2
        ),
      
      both_significant_genes =
        sum(
          author_comparison$
            both_significant == "Yes",
          na.rm = TRUE
        )
    )
    
    write_csv(
      comparison_summary,
      file.path(
        table_folder,
        "GSE254530_author_comparison_summary.csv"
      )
    )
    
    comparison_plot <- ggplot(
      author_comparison,
      aes(
        x =
          author_log2_fold_change,
        y =
          our_log2_fold_change
      )
    ) +
      geom_point(
        alpha = 0.5
      ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed"
      ) +
      labs(
        title =
          "Comparison with published DESeq2 results",
        subtitle = paste0(
          "Pearson r = ",
          round(
            correlation_value,
            3
          )
        ),
        x =
          "Published log2 fold change",
        y =
          "Reanalysis log2 fold change"
      ) +
      theme_classic()
    
    ggsave(
      filename = file.path(
        figure_folder,
        "GSE254530_author_result_comparison.png"
      ),
      plot = comparison_plot,
      width = 7,
      height = 6,
      dpi = 300
    )
    
  } else {
    
    warning(
      "Could not identify required columns ",
      "in the published DESeq2 table."
    )
  }
  
} else {
  
  warning(
    "Exactly one published DESeq2-result file ",
    "was not found."
  )
}


# ------------------------------------------------------------
# 33. Save R objects
# ------------------------------------------------------------

saveRDS(
  dds,
  file.path(
    output_folder,
    "GSE254530_DESeq2_object.rds"
  )
)

saveRDS(
  vsd,
  file.path(
    output_folder,
    "GSE254530_VST_object.rds"
  )
)

saveRDS(
  results_table,
  file.path(
    output_folder,
    "GSE254530_DESeq2_results_table.rds"
  )
)


# ------------------------------------------------------------
# 34. Save session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    output_folder,
    "sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 35. Final report
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 09 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat(
  "Genes before filtering:",
  nrow(count_matrix),
  "\n"
)

cat(
  "Genes after filtering:",
  nrow(filtered_count_matrix),
  "\n"
)

cat(
  "Upregulated genes:",
  sum(
    results_table$
      significance ==
      "Upregulated"
  ),
  "\n"
)

cat(
  "Downregulated genes:",
  sum(
    results_table$
      significance ==
      "Downregulated"
  ),
  "\n"
)

cat(
  "Not significant:",
  sum(
    results_table$
      significance ==
      "Not significant"
  ),
  "\n"
)

if (nrow(comparison_summary) == 1) {
  
  cat(
    "Genes compared with authors:",
    comparison_summary$
      genes_compared,
    "\n"
  )
  
  cat(
    "Log2FC correlation with authors:",
    round(
      comparison_summary$
        pearson_log2FC_correlation,
      3
    ),
    "\n"
  )
  
  cat(
    "Same-direction genes:",
    comparison_summary$
      same_direction_percent,
    "%\n"
  )
}

cat(
  "\nOutputs saved inside:\n",
  output_folder,
  "\n"
)
# ============================================================
# STEP 09 COMPLETION BLOCK
# Run after the previous script stopped at de_summary
# ============================================================

# Verify that required objects still exist
required_objects <- c(
  "results_table",
  "significant_results",
  "dds",
  "sample_design",
  "table_folder",
  "figure_folder",
  "output_folder"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Required objects are missing: ",
    paste(missing_objects, collapse = ", "),
    "\nRun the main Step 09 script again up to results_table."
  )
}

# ------------------------------------------------------------
# 1. Build DE summary without count() conflict
# ------------------------------------------------------------

de_summary <- results_table |>
  dplyr::group_by(significance) |>
  dplyr::summarise(
    number_of_genes = dplyr::n(),
    .groups = "drop"
  ) |>
  tidyr::complete(
    significance = c(
      "Upregulated",
      "Downregulated",
      "Not significant"
    ),
    fill = list(
      number_of_genes = 0L
    )
  ) |>
  dplyr::mutate(
    display_order = match(
      significance,
      c(
        "Upregulated",
        "Downregulated",
        "Not significant"
      )
    )
  ) |>
  dplyr::arrange(display_order) |>
  dplyr::select(-display_order)

readr::write_csv(
  de_summary,
  file.path(
    table_folder,
    "GSE254530_DE_summary.csv"
  )
)

print(de_summary)

# ------------------------------------------------------------
# 2. Save normalized counts
# ------------------------------------------------------------

normalized_counts <- DESeq2::counts(
  dds,
  normalized = TRUE
)

normalized_counts_table <- as.data.frame(
  normalized_counts
) |>
  tibble::rownames_to_column("gene_id") |>
  tibble::as_tibble()

readr::write_csv(
  normalized_counts_table,
  file.path(
    table_folder,
    "GSE254530_normalized_counts.csv"
  )
)

# ------------------------------------------------------------
# 3. VST transformation
# ------------------------------------------------------------

vsd <- DESeq2::vst(
  dds,
  blind = TRUE
)

vst_matrix <- SummarizedExperiment::assay(vsd)

readr::write_csv(
  as.data.frame(vst_matrix) |>
    tibble::rownames_to_column("gene_id"),
  file.path(
    table_folder,
    "GSE254530_VST_expression_matrix.csv"
  )
)

# ------------------------------------------------------------
# 4. PCA figure
# ------------------------------------------------------------

pca_data <- DESeq2::plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percent_variance <- round(
  100 * attr(pca_data, "percentVar"),
  1
)

pca_export <- pca_data |>
  tibble::rownames_to_column("sample_row")

pca_plot <- ggplot2::ggplot(
  pca_data,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    shape = condition,
    label = name
  )
) +
  ggplot2::geom_point(size = 4) +
  ggrepel::geom_text_repel(
    max.overlaps = Inf
  ) +
  ggplot2::labs(
    title = "GSE254530 RNA-seq PCA",
    subtitle = "Vancomycin-treated versus untreated",
    x = paste0(
      "PC1: ",
      percent_variance[1],
      "% variance"
    ),
    y = paste0(
      "PC2: ",
      percent_variance[2],
      "% variance"
    ),
    shape = "Condition"
  ) +
  ggplot2::theme_classic()

ggplot2::ggsave(
  filename = file.path(
    figure_folder,
    "GSE254530_PCA.png"
  ),
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

readr::write_csv(
  pca_export,
  file.path(
    table_folder,
    "GSE254530_PCA_coordinates.csv"
  )
)

# ------------------------------------------------------------
# 5. Sample-distance heatmap
# ------------------------------------------------------------

sample_distance_matrix <- as.matrix(
  stats::dist(
    t(vst_matrix)
  )
)

annotation_data <- data.frame(
  Condition = sample_design$condition
)

rownames(annotation_data) <- rownames(sample_design)

grDevices::png(
  filename = file.path(
    figure_folder,
    "GSE254530_sample_distance_heatmap.png"
  ),
  width = 1800,
  height = 1600,
  res = 250
)

pheatmap::pheatmap(
  sample_distance_matrix,
  annotation_col = annotation_data,
  annotation_row = annotation_data,
  main = "GSE254530 sample-to-sample distance"
)

grDevices::dev.off()

# ------------------------------------------------------------
# 6. Volcano plot
# ------------------------------------------------------------

volcano_data <- results_table |>
  dplyr::mutate(
    volcano_adjusted_p = dplyr::case_when(
      is.na(adjusted_p_value) ~ 1,
      adjusted_p_value <= 0 ~ .Machine$double.xmin,
      TRUE ~ adjusted_p_value
    )
  )

label_genes <- significant_results |>
  dplyr::filter(
    !is.na(adjusted_p_value)
  ) |>
  dplyr::arrange(
    adjusted_p_value
  ) |>
  dplyr::slice_head(
    n = min(20, nrow(significant_results))
  ) |>
  dplyr::pull(gene_id)

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  volcano_data,
  lab = volcano_data$gene_id,
  x = "log2_fold_change",
  y = "volcano_adjusted_p",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Vancomycin-treated versus untreated",
  subtitle = "GSE254530 RNA-seq",
  caption = "Adjusted P < 0.05 and |log2FC| >= 1",
  pointSize = 2,
  labSize = 3,
  selectLab = label_genes
)

ggplot2::ggsave(
  filename = file.path(
    figure_folder,
    "GSE254530_volcano_plot.png"
  ),
  plot = volcano_plot,
  width = 9,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 7. Heatmap of top 50 significant genes
# ------------------------------------------------------------

top_genes <- significant_results |>
  dplyr::filter(
    !is.na(adjusted_p_value)
  ) |>
  dplyr::arrange(
    adjusted_p_value
  ) |>
  dplyr::slice_head(
    n = min(
      50,
      nrow(significant_results)
    )
  ) |>
  dplyr::pull(gene_id)

if (length(top_genes) >= 2) {
  
  top_heatmap_matrix <- vst_matrix[
    top_genes,
    ,
    drop = FALSE
  ]
  
  top_heatmap_matrix <- t(
    scale(
      t(top_heatmap_matrix)
    )
  )
  
  top_heatmap_matrix[
    is.na(top_heatmap_matrix)
  ] <- 0
  
  grDevices::png(
    filename = file.path(
      figure_folder,
      "GSE254530_top50_DEG_heatmap.png"
    ),
    width = 1800,
    height = 2200,
    res = 250
  )
  
  pheatmap::pheatmap(
    top_heatmap_matrix,
    annotation_col = annotation_data,
    show_rownames = TRUE,
    fontsize_row = 6,
    main = "Top vancomycin-responsive genes"
  )
  
  grDevices::dev.off()
}

# ------------------------------------------------------------
# 8. Save R objects
# ------------------------------------------------------------

saveRDS(
  dds,
  file.path(
    output_folder,
    "GSE254530_DESeq2_object.rds"
  )
)

saveRDS(
  vsd,
  file.path(
    output_folder,
    "GSE254530_VST_object.rds"
  )
)

saveRDS(
  results_table,
  file.path(
    output_folder,
    "GSE254530_DESeq2_results_table.rds"
  )
)

# ------------------------------------------------------------
# 9. Save session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    output_folder,
    "sessionInfo.txt"
  )
)

# ------------------------------------------------------------
# 10. Confirm created figures
# ------------------------------------------------------------

created_figures <- list.files(
  figure_folder,
  pattern = "\\.png$",
  full.names = TRUE
)

figure_report <- tibble::tibble(
  figure = basename(created_figures),
  exists = file.exists(created_figures),
  size_KB = round(
    file.info(created_figures)$size / 1024,
    2
  )
)

readr::write_csv(
  figure_report,
  file.path(
    table_folder,
    "GSE254530_figure_creation_report.csv"
  )
)

cat("\n============================================\n")
cat("STEP 09 COMPLETION FINISHED\n")
cat("============================================\n")

cat(
  "Upregulated genes:",
  sum(
    results_table$significance ==
      "Upregulated"
  ),
  "\n"
)

cat(
  "Downregulated genes:",
  sum(
    results_table$significance ==
      "Downregulated"
  ),
  "\n"
)

cat(
  "Figures created:",
  nrow(figure_report),
  "\n\n"
)

print(figure_report)

cat(
  "\nFigures saved inside:\n",
  figure_folder,
  "\n"
)