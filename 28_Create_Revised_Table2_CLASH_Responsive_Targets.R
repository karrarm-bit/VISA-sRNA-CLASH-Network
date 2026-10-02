# ============================================================
# SCRIPT 28
# Create Revised Table 2:
# CLASH-supported regulatory RNA-mRNA interactions involving
# vancomycin-responsive target genes
#
# Revised evidence-layered framework
# ============================================================

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "tidyr",
  "openxlsx",
  "officer",
  "flextable"
)

new_packages <- required_packages[
  !(required_packages %in% rownames(installed.packages()))
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(openxlsx)
  library(officer)
  library(flextable)
})


# ------------------------------------------------------------
# 1. Project directories
# ------------------------------------------------------------

project_dir <- "D:/Bac-sRNA"

input_file <- file.path(
  project_dir,
  "14_revised_evidence_layered_network",
  "tables",
  "Evidence_Layered_CLASH_Network_ALL.csv"
)

outdir <- file.path(
  project_dir,
  "28_Revised_Table2"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Check main input
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop(
    paste0(
      "\nERROR: Main revised network file not found:\n",
      input_file,
      "\n"
    )
  )
}

cat("\nReading revised evidence-layered CLASH network...\n")

dat <- read_csv(
  input_file,
  show_col_types = FALSE
)

cat("Rows in complete revised network: ", nrow(dat), "\n")
cat("Columns: ", ncol(dat), "\n\n")

print(names(dat))


# ------------------------------------------------------------
# 3. Validate required columns
# ------------------------------------------------------------

required_cols <- c(
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
  "transcriptomic_support_score"
)

missing_cols <- setdiff(required_cols, names(dat))

if (length(missing_cols) > 0) {
  
  cat("\nMissing required columns:\n")
  print(missing_cols)
  
  stop(
    "\nERROR: Required columns are missing from revised network.\n"
  )
}

cat("\nAll required revised-network columns detected.\n")


# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------

clean_chr <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}


format_number <- function(x, digits = 3) {
  
  x <- suppressWarnings(as.numeric(x))
  
  ifelse(
    is.na(x),
    "",
    formatC(
      x,
      format = "f",
      digits = digits
    )
  )
}


format_pvalue <- function(x) {
  
  x <- suppressWarnings(as.numeric(x))
  
  ifelse(
    is.na(x),
    "",
    ifelse(
      x == 0,
      "<1 × 10^-300",
      format(
        x,
        scientific = TRUE,
        digits = 3
      )
    )
  )
}


normalize_id <- function(x) {
  
  x <- clean_chr(x)
  
  x <- str_replace_all(
    x,
    "\\s+",
    ""
  )
  
  toupper(x)
}


# ------------------------------------------------------------
# 5. Standardize main data
# ------------------------------------------------------------

dat <- dat %>%
  mutate(
    source_label = clean_chr(source_label),
    source_RNA_class = clean_chr(source_RNA_class),
    
    target_node = clean_chr(target_node),
    target_label = clean_chr(target_label),
    target_raw_names = clean_chr(target_raw_names),
    
    total_hybrid_count =
      suppressWarnings(as.numeric(total_hybrid_count)),
    
    maximum_number_of_experiments =
      suppressWarnings(as.numeric(maximum_number_of_experiments)),
    
    best_adjusted_p_value =
      suppressWarnings(as.numeric(best_adjusted_p_value)),
    
    best_connection_score =
      suppressWarnings(as.numeric(best_connection_score)),
    
    clash_evidence_score =
      suppressWarnings(as.numeric(clash_evidence_score)),
    
    target_log2_fold_change =
      suppressWarnings(as.numeric(target_log2_fold_change)),
    
    target_deseq_adjusted_p =
      suppressWarnings(as.numeric(target_deseq_adjusted_p))
  )


# ------------------------------------------------------------
# 6. Select significant target-DEG interactions
#
# IMPORTANT:
# Target differential expression is NOT required for retention
# in the complete CLASH network.
#
# Table 2 is specifically the responsive-target subset.
# ------------------------------------------------------------

responsive <- dat %>%
  filter(
    target_is_significant_deg %in% TRUE
  )

cat(
  "\nCLASH interactions involving significant target DEGs: ",
  nrow(responsive),
  "\n"
)

cat(
  "Unique responsive target genes: ",
  n_distinct(responsive$target_label),
  "\n"
)


# ------------------------------------------------------------
# 7. Check expected counts
# ------------------------------------------------------------

if (nrow(responsive) != 10) {
  
  warning(
    paste0(
      "Expected 10 responsive interactions, but detected ",
      nrow(responsive),
      ". Check upstream data."
    )
  )
}

if (n_distinct(responsive$target_label) != 7) {
  
  warning(
    paste0(
      "Expected 7 unique responsive target genes, but detected ",
      n_distinct(responsive$target_label),
      "."
    )
  )
}


# ------------------------------------------------------------
# 8. Conservative RNA-class renaming
#
# Critical revision:
# Do NOT call RS04120 / RS05585 "3'UTR-derived regulatory sRNA".
# Term-seq did not establish independent processing.
# ------------------------------------------------------------

responsive <- responsive %>%
  mutate(
    
    Regulatory_RNA = case_when(
      
      str_detect(
        source_label,
        regex("RS04120", ignore_case = TRUE)
      ) ~ "3UTR-RS04120",
      
      str_detect(
        source_label,
        regex("RS05585", ignore_case = TRUE)
      ) ~ "3UTR-RS05585",
      
      TRUE ~ source_label
    ),
    
    Regulatory_RNA_class = case_when(
      
      str_detect(
        source_label,
        regex("RS04120|RS05585", ignore_case = TRUE)
      ) ~ "Putative 3′UTR-associated RNA candidate",
      
      source_RNA_class == "Named sRNA" ~
        "Named sRNA",
      
      source_RNA_class ==
        "Putative 3UTR-associated RNA" ~
        "Putative 3′UTR-associated RNA candidate",
      
      source_RNA_class ==
        "Putative 5UTR-associated RNA" ~
        "Putative 5′UTR-associated RNA candidate",
      
      source_RNA_class ==
        "Intergenic RNA candidate" ~
        "Intergenic RNA candidate",
      
      source_RNA_class ==
        "Numbered sRNA candidate" ~
        "Predicted sRNA candidate",
      
      TRUE ~ source_RNA_class
    )
  )


# ------------------------------------------------------------
# 9. Search automatically for a strain-specific annotation file
#
# We need:
#   Target locus
#   Predicted protein product
#
# Search likely CSV/XLSX files in project.
# ------------------------------------------------------------

cat("\nSearching for strain-specific annotation resources...\n")

all_files <- list.files(
  project_dir,
  recursive = TRUE,
  full.names = TRUE
)

annotation_candidates <- all_files[
  str_detect(
    basename(all_files),
    regex(
      "master|knowledge|annotation|annotated|JKD6008",
      ignore_case = TRUE
    )
  ) &
    str_detect(
      all_files,
      regex("\\.(csv|xlsx)$", ignore_case = TRUE)
    )
]

cat(
  "Potential annotation files found: ",
  length(annotation_candidates),
  "\n"
)


# ------------------------------------------------------------
# 10. Function to read candidate annotation file
# ------------------------------------------------------------

safe_read_table <- function(f) {
  
  tryCatch({
    
    if (
      str_detect(
        f,
        regex("\\.csv$", ignore_case = TRUE)
      )
    ) {
      
      read_csv(
        f,
        show_col_types = FALSE
      )
      
    } else if (
      str_detect(
        f,
        regex("\\.xlsx$", ignore_case = TRUE)
      )
    ) {
      
      openxlsx::read.xlsx(f)
      
    } else {
      
      NULL
    }
    
  }, error = function(e) {
    
    NULL
  })
}


# ------------------------------------------------------------
# 11. Flexible annotation-column detector
# ------------------------------------------------------------

find_existing_col <- function(df, candidates) {
  
  hit <- candidates[
    candidates %in% names(df)
  ]
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  hit[1]
}


# ------------------------------------------------------------
# 12. Build annotation lookup from available project resources
# ------------------------------------------------------------

annotation_lookup <- tibble()

for (f in annotation_candidates) {
  
  tmp <- safe_read_table(f)
  
  if (is.null(tmp)) next
  if (nrow(tmp) == 0) next
  
  locus_col <- find_existing_col(
    tmp,
    c(
      "current_locus_tag",
      "locus_tag",
      "Locus_tag",
      "locus",
      "Target locus",
      "target_locus",
      "gene_locus"
    )
  )
  
  gene_col <- find_existing_col(
    tmp,
    c(
      "gene_symbol",
      "gene",
      "Gene",
      "symbol",
      "gene_name",
      "master_label"
    )
  )
  
  protein_col <- find_existing_col(
    tmp,
    c(
      "protein_id",
      "refseq_protein",
      "protein_accession",
      "RefSeq_protein",
      "protein"
    )
  )
  
  product_col <- find_existing_col(
    tmp,
    c(
      "product",
      "gene_product",
      "protein_product",
      "product_description",
      "Predicted protein product",
      "description"
    )
  )
  
  # File must contain at least one identifier
  # plus useful annotation.
  has_identifier <-
    !is.na(locus_col) ||
    !is.na(gene_col) ||
    !is.na(protein_col)
  
  has_annotation <-
    !is.na(locus_col) ||
    !is.na(product_col)
  
  if (!has_identifier || !has_annotation) next
  
  
  n <- nrow(tmp)
  
  temp_lookup <- tibble(
    locus = if (!is.na(locus_col)) {
      clean_chr(tmp[[locus_col]])
    } else {
      rep(NA_character_, n)
    },
    
    gene = if (!is.na(gene_col)) {
      clean_chr(tmp[[gene_col]])
    } else {
      rep(NA_character_, n)
    },
    
    protein = if (!is.na(protein_col)) {
      clean_chr(tmp[[protein_col]])
    } else {
      rep(NA_character_, n)
    },
    
    product = if (!is.na(product_col)) {
      clean_chr(tmp[[product_col]])
    } else {
      rep(NA_character_, n)
    },
    
    annotation_source = basename(f)
  )
  
  annotation_lookup <- bind_rows(
    annotation_lookup,
    temp_lookup
  )
}


# ------------------------------------------------------------
# 13. Remove empty annotation rows
# ------------------------------------------------------------

if (nrow(annotation_lookup) > 0) {
  
  annotation_lookup <- annotation_lookup %>%
    filter(
      !(
        is.na(locus) &
          is.na(gene) &
          is.na(protein)
      )
    ) %>%
    distinct()
  
  cat(
    "Annotation lookup rows assembled: ",
    nrow(annotation_lookup),
    "\n"
  )
  
} else {
  
  cat(
    "\nWARNING: No suitable annotation table was automatically detected.\n"
  )
}


# ------------------------------------------------------------
# 14. Build normalized annotation keys
# ------------------------------------------------------------

if (nrow(annotation_lookup) > 0) {
  
  annotation_long <- bind_rows(
    
    annotation_lookup %>%
      filter(!is.na(locus)) %>%
      transmute(
        key = normalize_id(locus),
        locus,
        gene,
        protein,
        product,
        annotation_source
      ),
    
    annotation_lookup %>%
      filter(!is.na(gene)) %>%
      transmute(
        key = normalize_id(gene),
        locus,
        gene,
        protein,
        product,
        annotation_source
      ),
    
    annotation_lookup %>%
      filter(!is.na(protein)) %>%
      transmute(
        key = normalize_id(protein),
        locus,
        gene,
        protein,
        product,
        annotation_source
      )
    
  ) %>%
    filter(!is.na(key)) %>%
    distinct(key, .keep_all = TRUE)
  
} else {
  
  annotation_long <- tibble(
    key = character(),
    locus = character(),
    gene = character(),
    protein = character(),
    product = character(),
    annotation_source = character()
  )
}


# ------------------------------------------------------------
# 15. Create target identifiers for matching
# ------------------------------------------------------------

responsive <- responsive %>%
  mutate(
    key_target_label =
      normalize_id(target_label),
    
    key_target_node =
      normalize_id(target_node)
  )


# ------------------------------------------------------------
# 16. First annotation match: target label
# ------------------------------------------------------------

ann_gene <- annotation_long %>%
  select(
    key,
    locus_gene = locus,
    product_gene = product,
    annotation_source_gene = annotation_source
  )

responsive <- responsive %>%
  left_join(
    ann_gene,
    by = c(
      "key_target_label" = "key"
    )
  )


# ------------------------------------------------------------
# 17. Second annotation match: target node
# ------------------------------------------------------------

ann_node <- annotation_long %>%
  select(
    key,
    locus_node = locus,
    product_node = product,
    annotation_source_node = annotation_source
  )

responsive <- responsive %>%
  left_join(
    ann_node,
    by = c(
      "key_target_node" = "key"
    )
  )


# ------------------------------------------------------------
# 18. Combine annotation matches
# ------------------------------------------------------------

responsive <- responsive %>%
  mutate(
    
    Target_locus = coalesce(
      locus_gene,
      locus_node
    ),
    
    Predicted_protein_product = coalesce(
      product_gene,
      product_node
    ),
    
    Annotation_source = coalesce(
      annotation_source_gene,
      annotation_source_node
    )
  )


# ------------------------------------------------------------
# 19. Safe fallback for known responsive targets
#
# These values correspond to the strain-specific annotations
# already used in the previous manuscript table.
#
# They are used ONLY if automatic annotation lookup did not
# resolve the field.
# ------------------------------------------------------------

responsive <- responsive %>%
  mutate(
    
    Target_locus = case_when(
      
      !is.na(Target_locus) ~ Target_locus,
      
      target_label == "spa" ~
        "SAA6008_00089",
      
      target_label == "modA" ~
        "SAA6008_02316",
      
      target_label == "gltB" ~
        "SAA6008_00475",
      
      target_label == "guaB" ~
        "SAA6008_00389",
      
      target_label == "dhaK" ~
        "SAA6008_00666",
      
      target_label == "nrdF" ~
        "SAA6008_00747",
      
      target_label == "qoxB" ~
        "SAA6008_01015",
      
      TRUE ~ NA_character_
    ),
    
    Predicted_protein_product = case_when(
      
      !is.na(Predicted_protein_product) ~
        Predicted_protein_product,
      
      target_label == "spa" ~
        "staphylococcal protein A",
      
      target_label == "modA" ~
        "molybdate ABC transporter substrate-binding protein",
      
      target_label == "gltB" ~
        "glutamate synthase large subunit",
      
      target_label == "guaB" ~
        "IMP dehydrogenase",
      
      target_label == "dhaK" ~
        "dihydroxyacetone kinase subunit DhaK",
      
      target_label == "nrdF" ~
        "class 1b ribonucleoside-diphosphate reductase subunit beta",
      
      target_label == "qoxB" ~
        "cytochrome aa3 quinol oxidase subunit I",
      
      TRUE ~ NA_character_
    )
  )


# ------------------------------------------------------------
# 20. Define target response
# ------------------------------------------------------------

responsive <- responsive %>%
  mutate(
    
    Target_response = case_when(
      
      target_log2_fold_change > 0 ~
        "Induced",
      
      target_log2_fold_change < 0 ~
        "Repressed",
      
      TRUE ~
        "No directional change"
    )
  )


# ------------------------------------------------------------
# 21. IMPORTANT:
# Use revised CLASH support score/category ONLY.
#
# Do NOT use:
# legacy_comparable_total_score
# legacy_confidence_category
# ------------------------------------------------------------

responsive <- responsive %>%
  mutate(
    
    CLASH_support_score =
      clash_evidence_score,
    
    CLASH_support_category =
      clash_support_category
  )


# ------------------------------------------------------------
# 22. Sort table
#
# CLASH evidence first.
# Transcriptomic response remains a separate evidence layer.
# ------------------------------------------------------------

responsive <- responsive %>%
  arrange(
    desc(CLASH_support_score),
    desc(total_hybrid_count),
    best_adjusted_p_value
  )


# ------------------------------------------------------------
# 23. Build manuscript Table 2
# ------------------------------------------------------------

table2 <- responsive %>%
  transmute(
    
    `Regulatory RNA` =
      Regulatory_RNA,
    
    `Regulatory RNA class` =
      Regulatory_RNA_class,
    
    `Target gene` =
      target_label,
    
    `Target locus` =
      Target_locus,
    
    `Predicted protein product` =
      Predicted_protein_product,
    
    `Target log2FC` =
      round(
        target_log2_fold_change,
        3
      ),
    
    `Target response` =
      Target_response,
    
    `Target adjusted P` =
      target_deseq_adjusted_p,
    
    `CLASH hybrid count` =
      as.integer(
        total_hybrid_count
      ),
    
    `CLASH adjusted P` =
      best_adjusted_p_value,
    
    `Connection score` =
      round(
        best_connection_score,
        4
      ),
    
    `Experiments` =
      as.integer(
        maximum_number_of_experiments
      ),
    
    `CLASH-support score` =
      as.integer(
        CLASH_support_score
      ),
    
    `CLASH-support category` =
      CLASH_support_category
  )


# ------------------------------------------------------------
# 24. Check candidate terminology
# ------------------------------------------------------------

cat("\nChecking revised 3'UTR terminology...\n")

print(
  table2 %>%
    filter(
      str_detect(
        `Regulatory RNA`,
        regex(
          "RS04120|RS05585",
          ignore_case = TRUE
        )
      )
    )
)


# ------------------------------------------------------------
# 25. Ensure old confidence terminology is absent
# ------------------------------------------------------------

old_terms <- c(
  "High confidence",
  "Moderate confidence",
  "Supported candidate",
  "Exploratory interaction"
)

old_term_detected <- any(
  unlist(table2) %in% old_terms,
  na.rm = TRUE
)

if (old_term_detected) {
  
  warning(
    "Old confidence terminology was detected in revised Table 2."
  )
  
} else {
  
  cat(
    "\nPASS: Old High/Moderate/Supported confidence labels are absent.\n"
  )
}


# ------------------------------------------------------------
# 26. Export raw CSV
# ------------------------------------------------------------

csv_file <- file.path(
  outdir,
  "Table2_Revised_CLASH_Responsive_Interactions.csv"
)

write_csv(
  table2,
  csv_file
)


# ------------------------------------------------------------
# 27. Excel version
# ------------------------------------------------------------

xlsx_file <- file.path(
  outdir,
  "Table2_Revised_CLASH_Responsive_Interactions.xlsx"
)

wb <- createWorkbook()

addWorksheet(
  wb,
  "Revised Table 2"
)

writeData(
  wb,
  sheet = "Revised Table 2",
  x = table2,
  startRow = 1,
  startCol = 1
)

header_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 10,
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "TopBottom",
  borderStyle = "thin"
)

body_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 9,
  valign = "center",
  wrapText = TRUE
)

addStyle(
  wb,
  sheet = "Revised Table 2",
  style = header_style,
  rows = 1,
  cols = 1:ncol(table2),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Revised Table 2",
  style = body_style,
  rows = 2:(nrow(table2) + 1),
  cols = 1:ncol(table2),
  gridExpand = TRUE
)

freezePane(
  wb,
  sheet = "Revised Table 2",
  firstRow = TRUE
)

setColWidths(
  wb,
  sheet = "Revised Table 2",
  cols = 1:ncol(table2),
  widths = "auto"
)

saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)


# ------------------------------------------------------------
# 28. Prepare formatted Word table
# ------------------------------------------------------------

table2_word <- table2 %>%
  mutate(
    
    `Target log2FC` =
      format_number(
        `Target log2FC`,
        3
      ),
    
    `Target adjusted P` =
      format_pvalue(
        `Target adjusted P`
      ),
    
    `CLASH adjusted P` =
      format_pvalue(
        `CLASH adjusted P`
      ),
    
    `Connection score` =
      format_number(
        `Connection score`,
        4
      )
  )


# ------------------------------------------------------------
# 29. Build flextable
# ------------------------------------------------------------

ft <- flextable(
  table2_word
)

ft <- theme_booktabs(ft)

ft <- font(
  ft,
  fontname = "Times New Roman",
  part = "all"
)

ft <- fontsize(
  ft,
  size = 8,
  part = "body"
)

ft <- fontsize(
  ft,
  size = 8,
  part = "header"
)

ft <- bold(
  ft,
  bold = TRUE,
  part = "header"
)

ft <- align(
  ft,
  align = "center",
  part = "header"
)

ft <- align(
  ft,
  j = c(
    "Regulatory RNA",
    "Regulatory RNA class",
    "Predicted protein product"
  ),
  align = "left",
  part = "body"
)

ft <- align(
  ft,
  j = setdiff(
    names(table2_word),
    c(
      "Regulatory RNA",
      "Regulatory RNA class",
      "Predicted protein product"
    )
  ),
  align = "center",
  part = "body"
)

ft <- valign(
  ft,
  valign = "center",
  part = "all"
)

ft <- autofit(ft)

ft <- set_table_properties(
  ft,
  layout = "autofit",
  width = 1
)


# ------------------------------------------------------------
# 30. Table caption
# ------------------------------------------------------------

caption_text <- paste0(
  "Table 2. CLASH-supported regulatory RNA–mRNA interactions ",
  "involving vancomycin-responsive target genes."
)


# ------------------------------------------------------------
# 31. Revised manuscript-safe note
# ------------------------------------------------------------

note_text <- paste0(
  
  "Note: Table 2 shows the subset of CLASH-supported regulatory RNA–mRNA ",
  "interactions whose target genes met the predefined differential-expression ",
  "criteria in GSE254530. Target differential expression was treated as an ",
  "independent transcriptomic evidence layer and was not required for retention ",
  "of an interaction in the complete RNase III-CLASH network. CLASH evidence ",
  "was obtained from GSE254532. The CLASH-support score summarizes physical ",
  "interaction evidence from supporting observations, experimental recurrence, ",
  "hybrid abundance, connection score, and CLASH-adjusted P value; target ",
  "transcriptional response was not included in this score. CLASH-support ",
  "categories are descriptive evidence strata rather than functional-validation ",
  "classes. RS04120 and RS05585 are conservatively described as putative ",
  "3′UTR-associated RNA candidates because genomic context is compatible with ",
  "3′UTR association, whereas Term-seq did not demonstrate a reproducible sharp ",
  "candidate 3′ boundary. Accordingly, independent RNA processing and direct ",
  "regulation of nrdF remain unproven."
)


# ------------------------------------------------------------
# 32. Create Word document
# ------------------------------------------------------------

doc <- read_docx()

doc <- body_add_par(
  doc,
  caption_text,
  style = "Normal"
)

doc <- body_add_flextable(
  doc,
  value = ft
)

doc <- body_add_par(
  doc,
  "",
  style = "Normal"
)

doc <- body_add_par(
  doc,
  note_text,
  style = "Normal"
)


# ------------------------------------------------------------
# 33. Landscape orientation
# ------------------------------------------------------------

doc <- body_end_section_landscape(
  doc
)


# ------------------------------------------------------------
# 34. Save Word document
# ------------------------------------------------------------

word_file <- file.path(
  outdir,
  "Table2_Revised_CLASH_Responsive_Interactions.docx"
)

print(
  doc,
  target = word_file
)


# ------------------------------------------------------------
# 35. Export annotation audit
# ------------------------------------------------------------

annotation_audit <- responsive %>%
  select(
    Regulatory_RNA,
    target_label,
    target_node,
    Target_locus,
    Predicted_protein_product,
    Annotation_source
  )

write_csv(
  annotation_audit,
  file.path(
    outdir,
    "Table2_Target_Annotation_Audit.csv"
  )
)


# ------------------------------------------------------------
# 36. Export evidence audit
# ------------------------------------------------------------

evidence_audit <- responsive %>%
  select(
    Regulatory_RNA,
    Regulatory_RNA_class,
    target_label,
    total_hybrid_count,
    maximum_number_of_experiments,
    best_adjusted_p_value,
    best_connection_score,
    clash_evidence_score,
    clash_support_category,
    target_log2_fold_change,
    target_deseq_adjusted_p,
    target_is_significant_deg,
    target_transcriptomic_layer
  )

write_csv(
  evidence_audit,
  file.path(
    outdir,
    "Table2_Evidence_Layer_Audit.csv"
  )
)


# ------------------------------------------------------------
# 37. Verification
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("SCRIPT 28 COMPLETED\n")
cat("============================================================\n\n")

cat(
  "Complete CLASH network = ",
  nrow(dat),
  "\n"
)

cat(
  "Table 2 interactions = ",
  nrow(table2),
  "\n"
)

cat(
  "Unique responsive target genes = ",
  n_distinct(responsive$target_label),
  "\n"
)

cat(
  "RS04120 rows = ",
  sum(
    str_detect(
      table2$`Regulatory RNA`,
      regex(
        "RS04120",
        ignore_case = TRUE
      )
    )
  ),
  "\n"
)

cat(
  "RS05585 rows = ",
  sum(
    str_detect(
      table2$`Regulatory RNA`,
      regex(
        "RS05585",
        ignore_case = TRUE
      )
    )
  ),
  "\n"
)


cat("\nCLASH-support categories:\n")

print(
  table(
    table2$`CLASH-support category`,
    useNA = "ifany"
  )
)


cat("\nTarget genes:\n")

print(
  table(
    table2$`Target gene`
  )
)


cat("\nFinal revised Table 2:\n")

print(
  table2,
  n = Inf
)


# ------------------------------------------------------------
# 38. Save verification summary
# ------------------------------------------------------------

verification <- tibble(
  Check = c(
    "Complete CLASH network interactions",
    "Responsive-target Table 2 interactions",
    "Unique responsive target genes",
    "RS04120 rows",
    "RS05585 rows",
    "Old confidence terminology detected"
  ),
  
  Value = c(
    nrow(dat),
    nrow(table2),
    n_distinct(responsive$target_label),
    
    sum(
      str_detect(
        table2$`Regulatory RNA`,
        regex("RS04120", ignore_case = TRUE)
      )
    ),
    
    sum(
      str_detect(
        table2$`Regulatory RNA`,
        regex("RS05585", ignore_case = TRUE)
      )
    ),
    
    old_term_detected
  )
)

write_csv(
  verification,
  file.path(
    outdir,
    "Script28_Output_Verification.csv"
  )
)


# ------------------------------------------------------------
# 39. Save session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "Script28_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 40. Final output paths
# ------------------------------------------------------------

cat("\nOutputs:\n")

cat(
  "WORD:\n",
  word_file,
  "\n\n"
)

cat(
  "EXCEL:\n",
  xlsx_file,
  "\n\n"
)

cat(
  "CSV:\n",
  csv_file,
  "\n\n"
)

cat(
  "Annotation audit:\n",
  file.path(
    outdir,
    "Table2_Target_Annotation_Audit.csv"
  ),
  "\n\n"
)

cat(
  "Evidence audit:\n",
  file.path(
    outdir,
    "Table2_Evidence_Layer_Audit.csv"
  ),
  "\n"
)

cat("\n============================================================\n")