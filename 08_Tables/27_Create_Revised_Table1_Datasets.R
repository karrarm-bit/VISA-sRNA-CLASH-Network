# ============================================================
# SCRIPT 27
# CREATE REVISED TABLE 1
# Characteristics and analytical roles of GEO datasets
# ============================================================

# 1. Packages -------------------------------------------------

packages <- c(
  "officer",
  "flextable",
  "openxlsx"
)

missing_packages <- packages[
  !packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(officer)
library(flextable)
library(openxlsx)


# 2. Output directory -----------------------------------------

outdir <- "D:/Bac-sRNA/27_Revised_Table1"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Table data -----------------------------------------------

table1 <- data.frame(
  
  `GEO series` = c(
    "GSE254530",
    "GSE254532",
    "GSE254531",
    "GSE158830"
  ),
  
  `Organism/strain` = c(
    "S. aureus JKD6008",
    "S. aureus JKD6008",
    "S. aureus JKD6008",
    "S. aureus JKD6008 and JKD6009"
  ),
  
  Assay = c(
    "RNA-seq",
    "RNase III-CLASH",
    "Ribo-seq",
    "RNA-seq, RNase III-CLASH, dRNA-seq, and Term-seq"
  ),
  
  `Experimental groups` = c(
    "Untreated vs treated",
    "Untreated, treated, and untagged technical controls",
    "Untreated baseline samples",
    "Multiple genotypes, untreated and treated conditions, plus technical controls"
  ),
  
  `Vancomycin exposure` = c(
    "8 µg mL⁻¹, 30 min",
    "8 µg mL⁻¹, 30 min in treated samples",
    "None",
    "3 µg mL⁻¹, 10 min for treated RNA-seq samples"
  ),
  
  `Sample size` = c(
    "6",
    "10",
    "2",
    "39"
  ),
  
  `Analytical role` = c(
    "Primary differential-expression analysis",
    "Reconstruction of the RNase III-CLASH interaction network",
    "Baseline translational context only",
    "Supporting transcript annotation and Term-seq-based candidate boundary assessment"
  ),
  
  Reference = c(
    "Wu, Pang et al., 2024",
    "Wu, Pang et al., 2024",
    "Wu, Pang et al., 2024",
    "Mediati, Wong et al., 2022"
  ),
  
  check.names = FALSE
)


# 4. Table caption and note -----------------------------------

caption_text <- paste0(
  "Table 1. Characteristics and analytical roles of the GEO ",
  "datasets included in the present study."
)

note_text <- paste0(
  "Note: Vancomycin-exposure details apply only to the corresponding ",
  "treated samples. For GSE158830, the reported exposure of 3 µg mL⁻¹ ",
  "for 10 min applies specifically to the treated RNA-seq samples. ",
  "RNA-seq, RNA sequencing; Ribo-seq, ribosome profiling; CLASH, ",
  "cross-linking, ligation and sequencing of hybrids; dRNA-seq, ",
  "differential RNA sequencing. The datasets were integrated as ",
  "complementary evidence layers according to their experimental roles ",
  "and were not treated as directly interchangeable. GSE254531 contained ",
  "untreated Ribo-seq samples only and was therefore considered as ",
  "baseline translational context rather than evidence of ",
  "vancomycin-dependent translational regulation."
)


# 5. Create formatted flextable -------------------------------

ft <- flextable(table1)

ft <- set_header_labels(
  ft,
  `GEO series` = "GEO series",
  `Organism/strain` = "Organism/strain",
  Assay = "Assay",
  `Experimental groups` = "Experimental groups",
  `Vancomycin exposure` = "Vancomycin exposure",
  `Sample size` = "Sample size",
  `Analytical role` = "Analytical role",
  Reference = "Ref."
)

# Font
ft <- font(
  ft,
  fontname = "Times New Roman",
  part = "all"
)

ft <- fontsize(
  ft,
  size = 9,
  part = "body"
)

ft <- fontsize(
  ft,
  size = 9,
  part = "header"
)

# Header formatting
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

ft <- valign(
  ft,
  valign = "center",
  part = "all"
)

# Body alignment
ft <- align(
  ft,
  j = c(
    "GEO series",
    "Sample size"
  ),
  align = "center",
  part = "body"
)

ft <- align(
  ft,
  j = c(
    "Organism/strain",
    "Assay",
    "Experimental groups",
    "Vancomycin exposure",
    "Analytical role",
    "Reference"
  ),
  align = "left",
  part = "body"
)

# Borders similar to manuscript table
border_black <- fp_border(
  color = "black",
  width = 0.75
)

border_light <- fp_border(
  color = "black",
  width = 0.35
)

ft <- border_remove(ft)

ft <- hline_top(
  ft,
  border = border_black,
  part = "header"
)

ft <- hline_bottom(
  ft,
  border = border_black,
  part = "header"
)

ft <- hline(
  ft,
  border = border_light,
  part = "body"
)

ft <- hline_bottom(
  ft,
  border = border_black,
  part = "body"
)

# Padding
ft <- padding(
  ft,
  padding.top = 3,
  padding.bottom = 3,
  padding.left = 3,
  padding.right = 3,
  part = "all"
)

# Widths tuned for landscape page
ft <- width(ft, j = "GEO series", width = 0.85)
ft <- width(ft, j = "Organism/strain", width = 1.25)
ft <- width(ft, j = "Assay", width = 1.35)
ft <- width(ft, j = "Experimental groups", width = 1.70)
ft <- width(ft, j = "Vancomycin exposure", width = 1.45)
ft <- width(ft, j = "Sample size", width = 0.65)
ft <- width(ft, j = "Analytical role", width = 2.15)
ft <- width(ft, j = "Reference", width = 1.15)

ft <- set_table_properties(
  ft,
  layout = "fixed",
  width = 1
)


# 6. Create Word document -------------------------------------

doc <- read_docx()

# Landscape orientation
doc <- body_end_section_landscape(doc)

# Caption
doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      caption_text,
      prop = fp_text(
        font.family = "Times New Roman",
        font.size = 10,
        bold = TRUE,
        italic = TRUE
      )
    )
  )
)

# Add table
doc <- body_add_flextable(
  doc,
  value = ft
)

# Add note
doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      note_text,
      prop = fp_text(
        font.family = "Times New Roman",
        font.size = 8,
        italic = TRUE
      )
    )
  )
)

word_file <- file.path(
  outdir,
  "Table1_Revised_GEO_Datasets.docx"
)

print(
  doc,
  target = word_file
)


# 7. Create Excel version -------------------------------------

xlsx_file <- file.path(
  outdir,
  "Table1_Revised_GEO_Datasets.xlsx"
)

wb <- createWorkbook()

addWorksheet(
  wb,
  "Table 1"
)

writeData(
  wb,
  sheet = "Table 1",
  x = table1,
  startRow = 3,
  startCol = 1
)

writeData(
  wb,
  sheet = "Table 1",
  x = caption_text,
  startRow = 1,
  startCol = 1
)

writeData(
  wb,
  sheet = "Table 1",
  x = note_text,
  startRow = 9,
  startCol = 1
)

mergeCells(
  wb,
  sheet = "Table 1",
  cols = 1:8,
  rows = 1
)

mergeCells(
  wb,
  sheet = "Table 1",
  cols = 1:8,
  rows = 9
)

header_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 10,
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "TopBottom"
)

body_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 10,
  valign = "center",
  wrapText = TRUE,
  border = "Bottom"
)

caption_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 10,
  textDecoration = "bold",
  italic = TRUE
)

note_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 9,
  italic = TRUE,
  wrapText = TRUE,
  valign = "top"
)

addStyle(
  wb,
  "Table 1",
  header_style,
  rows = 3,
  cols = 1:8,
  gridExpand = TRUE
)

addStyle(
  wb,
  "Table 1",
  body_style,
  rows = 4:7,
  cols = 1:8,
  gridExpand = TRUE
)

addStyle(
  wb,
  "Table 1",
  caption_style,
  rows = 1,
  cols = 1
)

addStyle(
  wb,
  "Table 1",
  note_style,
  rows = 9,
  cols = 1
)

setColWidths(
  wb,
  "Table 1",
  cols = 1:8,
  widths = c(
    14, 23, 25, 34,
    27, 12, 42, 23
  )
)

setRowHeights(
  wb,
  "Table 1",
  rows = 3:7,
  heights = 45
)

setRowHeights(
  wb,
  "Table 1",
  rows = 9,
  heights = 75
)

saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)


# 8. Export plain CSV for reproducibility ---------------------

write.csv(
  table1,
  file.path(
    outdir,
    "Table1_Revised_GEO_Datasets.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# 9. Console verification -------------------------------------

cat("\n")
cat("====================================================\n")
cat("REVISED TABLE 1 CREATED SUCCESSFULLY\n")
cat("====================================================\n\n")

cat("Word:\n")
cat(word_file, "\n\n")

cat("Excel:\n")
cat(xlsx_file, "\n\n")

cat("CSV:\n")
cat(
  file.path(
    outdir,
    "Table1_Revised_GEO_Datasets.csv"
  ),
  "\n\n"
)

cat("Rows =", nrow(table1), "\n")
cat("Columns =", ncol(table1), "\n")

cat("\nAnalytical roles:\n")
print(
  table1[, c(
    "GEO series",
    "Analytical role"
  )]
)

cat("\n====================================================\n")