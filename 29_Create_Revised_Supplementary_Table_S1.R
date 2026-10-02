# ============================================================
# SCRIPT 29
# Revised Supplementary Table S1
#
# CLASH-support scoring framework used to summarize evidence
# for regulatory RNA-mRNA interactions
#
# Outputs:
#   1. Word (.docx)
#   2. Excel (.xlsx)
#   3. CSV files
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "tibble",
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
  library(dplyr)
  library(tibble)
  library(openxlsx)
  library(officer)
  library(flextable)
})


# ------------------------------------------------------------
# 1. Output directory
# ------------------------------------------------------------

project_dir <- "D:/Bac-sRNA"

outdir <- file.path(
  project_dir,
  "29_Revised_Supplementary_Table_S1"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Main revised scoring framework
# ------------------------------------------------------------

scoring_table <- tribble(
  
  ~`Evidence component`,
  ~`Scoring condition`,
  ~`Points assigned`,
  
  "Supporting observations",
  "≥3 supporting CLASH observations",
  "2",
  
  "",
  "2 supporting CLASH observations",
  "1",
  
  "",
  "<2 supporting CLASH observations",
  "0",
  
  
  "Experimental recurrence",
  "Supported in ≥3 independent experiments",
  "2",
  
  "",
  "Supported in 2 independent experiments",
  "1",
  
  "",
  "Supported in <2 independent experiments",
  "0",
  
  
  "CLASH hybrid abundance",
  "≥50 total CLASH hybrids",
  "3",
  
  "",
  "20–49 total CLASH hybrids",
  "2",
  
  "",
  "5–19 total CLASH hybrids",
  "1",
  
  "",
  "<5 total CLASH hybrids",
  "0",
  
  
  "Connection score",
  "Connection score ≥0.50",
  "3",
  
  "",
  "Connection score ≥0.20 and <0.50",
  "2",
  
  "",
  "Connection score ≥0.05 and <0.20",
  "1",
  
  "",
  "Connection score <0.05 or unavailable",
  "0",
  
  
  "CLASH-adjusted P value",
  "Adjusted P ≤1 × 10⁻¹⁰",
  "3",
  
  "",
  "Adjusted P >1 × 10⁻¹⁰ and ≤0.001",
  "2",
  
  "",
  "Adjusted P >0.001 and <0.05",
  "1",
  
  "",
  "Adjusted P ≥0.05 or unavailable",
  "0"
)


# ------------------------------------------------------------
# 3. CLASH-support category definitions
# ------------------------------------------------------------

category_table <- tribble(
  
  ~`Total CLASH-support score`,
  ~`CLASH-support category`,
  
  "≥9",
  "Higher CLASH support",
  
  "6–8",
  "Intermediate CLASH support",
  
  "3–5",
  "Limited CLASH support",
  
  "<3",
  "Minimal CLASH support"
)


# ------------------------------------------------------------
# 4. Revised note
# ------------------------------------------------------------

note_text <- paste0(
  
  "Note: Each interaction was evaluated using a study-specific ",
  "CLASH-support framework incorporating the number of supporting ",
  "CLASH observations, experimental recurrence, total hybrid abundance, ",
  "connection score, and CLASH-adjusted P value. Component scores were ",
  "summed to generate a maximum CLASH-support score of 13 points. ",
  
  "Target transcriptional response and regulatory-RNA annotation status ",
  "were evaluated as separate evidence layers and were not included in ",
  "the CLASH-support score. ",
  
  "Ribo-seq data from GSE254531 were inspected descriptively as a ",
  "baseline translational resource and were not incorporated into the ",
  "CLASH-support score. ",
  
  "The resulting categories (Higher, Intermediate, Limited, and Minimal ",
  "CLASH support) describe relative support within the RNase III-CLASH ",
  "dataset and should not be interpreted as levels of biological or ",
  "functional validation."
)


# ------------------------------------------------------------
# 5. Export CSV files
# ------------------------------------------------------------

write.csv(
  scoring_table,
  file.path(
    outdir,
    "Supplementary_Table_S1_CLASH_Scoring_Framework.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  category_table,
  file.path(
    outdir,
    "Supplementary_Table_S1_CLASH_Support_Categories.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 6. Create Excel workbook
# ------------------------------------------------------------

wb <- createWorkbook()

addWorksheet(
  wb,
  "Supplementary Table S1"
)


# ------------------------------------------------------------
# 7. Excel title
# ------------------------------------------------------------

excel_title <- paste0(
  "Supplementary Table S1. ",
  "CLASH-support scoring framework used to summarize evidence ",
  "for regulatory RNA–mRNA interactions."
)

writeData(
  wb,
  sheet = "Supplementary Table S1",
  x = excel_title,
  startRow = 1,
  startCol = 1
)

mergeCells(
  wb,
  sheet = "Supplementary Table S1",
  cols = 1:3,
  rows = 1
)


# ------------------------------------------------------------
# 8. Write main scoring table
# ------------------------------------------------------------

writeData(
  wb,
  sheet = "Supplementary Table S1",
  x = scoring_table,
  startRow = 3,
  startCol = 1,
  headerStyle = createStyle(
    fontName = "Times New Roman",
    fontSize = 11,
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    wrapText = TRUE,
    border = "TopBottom",
    borderStyle = "thin"
  )
)


# ------------------------------------------------------------
# 9. Excel styles
# ------------------------------------------------------------

title_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 11,
  textDecoration = c("bold", "italic"),
  halign = "left",
  valign = "center",
  wrapText = TRUE
)

body_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 10,
  valign = "center",
  wrapText = TRUE
)

component_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 10,
  textDecoration = "bold",
  valign = "center",
  wrapText = TRUE
)

center_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 10,
  halign = "center",
  valign = "center",
  wrapText = TRUE
)

section_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 10,
  textDecoration = "bold",
  halign = "left",
  valign = "center"
)

note_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 9,
  textDecoration = "italic",
  valign = "top",
  wrapText = TRUE
)


# ------------------------------------------------------------
# 10. Apply Excel title style
# ------------------------------------------------------------

addStyle(
  wb,
  sheet = "Supplementary Table S1",
  style = title_style,
  rows = 1,
  cols = 1:3,
  gridExpand = TRUE
)


# ------------------------------------------------------------
# 11. Main table body style
# ------------------------------------------------------------

main_start_row <- 4
main_end_row <- 3 + nrow(scoring_table)

addStyle(
  wb,
  sheet = "Supplementary Table S1",
  style = body_style,
  rows = main_start_row:main_end_row,
  cols = 1:2,
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Supplementary Table S1",
  style = center_style,
  rows = main_start_row:main_end_row,
  cols = 3,
  gridExpand = TRUE
)


# ------------------------------------------------------------
# 12. Bold evidence-component rows
# ------------------------------------------------------------

component_rows <- which(
  scoring_table$`Evidence component` != ""
) + 3

addStyle(
  wb,
  sheet = "Supplementary Table S1",
  style = component_style,
  rows = component_rows,
  cols = 1,
  gridExpand = TRUE,
  stack = TRUE
)


# ------------------------------------------------------------
# 13. Category section
# ------------------------------------------------------------

category_heading_row <- main_end_row + 3

writeData(
  wb,
  sheet = "Supplementary Table S1",
  x = "CLASH-support category definitions",
  startRow = category_heading_row,
  startCol = 1
)

mergeCells(
  wb,
  sheet = "Supplementary Table S1",
  cols = 1:2,
  rows = category_heading_row
)

addStyle(
  wb,
  sheet = "Supplementary Table S1",
  style = section_style,
  rows = category_heading_row,
  cols = 1:2,
  gridExpand = TRUE
)


# ------------------------------------------------------------
# 14. Write category table
# ------------------------------------------------------------

category_table_start <- category_heading_row + 2

writeData(
  wb,
  sheet = "Supplementary Table S1",
  x = category_table,
  startRow = category_table_start,
  startCol = 1,
  headerStyle = createStyle(
    fontName = "Times New Roman",
    fontSize = 10,
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    wrapText = TRUE,
    border = "TopBottom",
    borderStyle = "thin"
  )
)

category_body_start <- category_table_start + 1

category_body_end <-
  category_table_start + nrow(category_table)

addStyle(
  wb,
  sheet = "Supplementary Table S1",
  style = center_style,
  rows = category_body_start:category_body_end,
  cols = 1:2,
  gridExpand = TRUE
)


# ------------------------------------------------------------
# 15. Add note to Excel
# ------------------------------------------------------------

note_row <- category_body_end + 3

writeData(
  wb,
  sheet = "Supplementary Table S1",
  x = note_text,
  startRow = note_row,
  startCol = 1
)

mergeCells(
  wb,
  sheet = "Supplementary Table S1",
  cols = 1:3,
  rows = note_row
)

addStyle(
  wb,
  sheet = "Supplementary Table S1",
  style = note_style,
  rows = note_row,
  cols = 1:3,
  gridExpand = TRUE
)

setRowHeights(
  wb,
  sheet = "Supplementary Table S1",
  rows = note_row,
  heights = 100
)


# ------------------------------------------------------------
# 16. Excel dimensions
# ------------------------------------------------------------

setColWidths(
  wb,
  sheet = "Supplementary Table S1",
  cols = 1,
  widths = 28
)

setColWidths(
  wb,
  sheet = "Supplementary Table S1",
  cols = 2,
  widths = 55
)

setColWidths(
  wb,
  sheet = "Supplementary Table S1",
  cols = 3,
  widths = 18
)

freezePane(
  wb,
  sheet = "Supplementary Table S1",
  firstActiveRow = 4
)


# ------------------------------------------------------------
# 17. Save Excel
# ------------------------------------------------------------

excel_file <- file.path(
  outdir,
  "Supplementary_Table_S1_Revised.xlsx"
)

saveWorkbook(
  wb,
  excel_file,
  overwrite = TRUE
)


# ============================================================
# WORD VERSION
# ============================================================


# ------------------------------------------------------------
# 18. Main scoring flextable
# ------------------------------------------------------------

ft1 <- flextable(scoring_table)

ft1 <- theme_booktabs(ft1)

ft1 <- font(
  ft1,
  fontname = "Times New Roman",
  part = "all"
)

ft1 <- fontsize(
  ft1,
  size = 9.5,
  part = "body"
)

ft1 <- fontsize(
  ft1,
  size = 9.5,
  part = "header"
)

ft1 <- bold(
  ft1,
  bold = TRUE,
  part = "header"
)

ft1 <- align(
  ft1,
  j = 1:3,
  align = "center",
  part = "header"
)

ft1 <- align(
  ft1,
  j = 1:2,
  align = "left",
  part = "body"
)

ft1 <- align(
  ft1,
  j = 3,
  align = "center",
  part = "body"
)

ft1 <- valign(
  ft1,
  valign = "center",
  part = "all"
)

# Bold first row of each evidence component
component_indices <- which(
  scoring_table$`Evidence component` != ""
)

for (i in component_indices) {
  
  ft1 <- bold(
    ft1,
    i = i,
    j = 1,
    bold = TRUE,
    part = "body"
  )
}

ft1 <- width(
  ft1,
  j = 1,
  width = 2.0
)

ft1 <- width(
  ft1,
  j = 2,
  width = 4.6
)

ft1 <- width(
  ft1,
  j = 3,
  width = 1.2
)

ft1 <- set_table_properties(
  ft1,
  layout = "fixed"
)


# ------------------------------------------------------------
# 19. Category flextable
# ------------------------------------------------------------

ft2 <- flextable(category_table)

ft2 <- theme_booktabs(ft2)

ft2 <- font(
  ft2,
  fontname = "Times New Roman",
  part = "all"
)

ft2 <- fontsize(
  ft2,
  size = 9.5,
  part = "all"
)

ft2 <- bold(
  ft2,
  bold = TRUE,
  part = "header"
)

ft2 <- align(
  ft2,
  align = "center",
  part = "all"
)

ft2 <- valign(
  ft2,
  valign = "center",
  part = "all"
)

ft2 <- width(
  ft2,
  j = 1,
  width = 2.8
)

ft2 <- width(
  ft2,
  j = 2,
  width = 3.8
)

ft2 <- set_table_properties(
  ft2,
  layout = "fixed"
)


# ------------------------------------------------------------
# 20. Word title
# ------------------------------------------------------------

word_title <- paste0(
  "Supplementary Table S1. ",
  "CLASH-support scoring framework used to summarize evidence ",
  "for regulatory RNA–mRNA interactions."
)


# ------------------------------------------------------------
# 21. Create Word document
# ------------------------------------------------------------

doc <- read_docx()

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      word_title,
      prop = fp_text(
        font.family = "Times New Roman",
        font.size = 10,
        bold = TRUE,
        italic = TRUE
      )
    )
  )
)

doc <- body_add_par(
  doc,
  "",
  style = "Normal"
)

doc <- body_add_flextable(
  doc,
  value = ft1
)

doc <- body_add_par(
  doc,
  "",
  style = "Normal"
)


# ------------------------------------------------------------
# 22. Category heading
# ------------------------------------------------------------

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      "CLASH-support category definitions",
      prop = fp_text(
        font.family = "Times New Roman",
        font.size = 10,
        bold = TRUE
      )
    )
  )
)

doc <- body_add_par(
  doc,
  "",
  style = "Normal"
)

doc <- body_add_flextable(
  doc,
  value = ft2
)

doc <- body_add_par(
  doc,
  "",
  style = "Normal"
)


# ------------------------------------------------------------
# 23. Add revised note
# ------------------------------------------------------------

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      note_text,
      prop = fp_text(
        font.family = "Times New Roman",
        font.size = 9,
        italic = TRUE
      )
    )
  )
)


# ------------------------------------------------------------
# 24. Save Word
# ------------------------------------------------------------

word_file <- file.path(
  outdir,
  "Supplementary_Table_S1_Revised.docx"
)

print(
  doc,
  target = word_file
)


# ------------------------------------------------------------
# 25. Verification file
# ------------------------------------------------------------

verification <- tibble(
  Item = c(
    "Maximum possible CLASH-support score",
    "Target transcriptional response included in score",
    "Regulatory-RNA annotation included in score",
    "Ribo-seq included in score",
    "Higher CLASH support threshold",
    "Intermediate CLASH support range",
    "Limited CLASH support range",
    "Minimal CLASH support threshold"
  ),
  
  Revised_value = c(
    "13",
    "No",
    "No",
    "No",
    "≥9",
    "6–8",
    "3–5",
    "<3"
  )
)

write.csv(
  verification,
  file.path(
    outdir,
    "Supplementary_Table_S1_Verification.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 26. Final console report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("REVISED SUPPLEMENTARY TABLE S1 COMPLETED\n")
cat("============================================================\n\n")

cat("Maximum CLASH-support score: 13\n")
cat("Target transcriptional response in score: NO\n")
cat("Regulatory-RNA annotation in score: NO\n")
cat("Ribo-seq in score: NO\n\n")

cat("Categories:\n")
cat("Higher CLASH support       = >=9\n")
cat("Intermediate CLASH support = 6-8\n")
cat("Limited CLASH support      = 3-5\n")
cat("Minimal CLASH support      = <3\n\n")

cat("WORD:\n")
cat(word_file, "\n\n")

cat("EXCEL:\n")
cat(excel_file, "\n\n")

cat("OUTPUT DIRECTORY:\n")
cat(outdir, "\n")

cat("\n============================================================\n")