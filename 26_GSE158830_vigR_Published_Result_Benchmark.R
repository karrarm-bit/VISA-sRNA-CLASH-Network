# ============================================================
# SCRIPT 26
# GSE158830 / vigR PUBLISHED-RESULT BENCHMARK
#
# Study:
# Mediati et al. - Nature Communications
# "RNase III-CLASH of multi-drug resistant Staphylococcus aureus
# reveals a regulatory mRNA 3′UTR required for intermediate
# vancomycin resistance"
#
# PURPOSE
# -------
# 1. Audit the locally available GSE158830 processed files.
# 2. Search specifically for vigR, isaA and folD.
# 3. Determine which published features can actually be tested
#    with the local data.
# 4. Separate:
#       Recovered
#       Partially recovered
#       Not recovered
#       Not tested
#
# IMPORTANT
# ---------
# This script does NOT treat absence from an RNA-seq table as
# evidence against an experimentally validated RNA-RNA interaction.
#
# It also does NOT claim reproduction of the complete original
# study, because the original work included CLASH, RNA-seq,
# dRNA-seq, Term-seq and experimental validation.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

packages <- c(
  "readxl",
  "readr",
  "dplyr",
  "stringr",
  "tidyr",
  "purrr",
  "openxlsx",
  "tibble"
)

missing_packages <- packages[
  !packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(readxl)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(openxlsx)
library(tibble)


# ============================================================
# 2. PROJECT PATHS
# ============================================================

project <- "D:/Bac-sRNA"

gse_dir <- file.path(
  project,
  "07_processed_data",
  "downloaded_files",
  "GSE158830"
)

outdir <- file.path(
  project,
  "26_GSE158830_vigR_benchmark"
)

tabdir <- file.path(
  outdir,
  "tables"
)

auditdir <- file.path(
  outdir,
  "audit"
)

dir.create(
  tabdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  auditdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. EXPECTED LOCAL FILES
# ============================================================

mut_vs_wt_file <- file.path(
  gse_dir,
  "GSE158830_2021.03.16_MUTvsWT.DEseq2_combined_output.xlsx"
)

mut_vs_repair_file <- file.path(
  gse_dir,
  "GSE158830_DEseq_output_vigR_MutvsREPAIR.xlsx"
)


input_files <- tibble(
  
  file_role = c(
    "vigR mutant versus WT DESeq2 output",
    "vigR mutant versus repair DESeq2 output"
  ),
  
  path = c(
    mut_vs_wt_file,
    mut_vs_repair_file
  ),
  
  exists = c(
    file.exists(mut_vs_wt_file),
    file.exists(mut_vs_repair_file)
  )
)


write_csv(
  input_files,
  file.path(
    auditdir,
    "Script26_Input_File_Audit.csv"
  )
)


cat("\nInput-file audit:\n")
print(input_files)


if (!all(input_files$exists)) {
  
  stop(
    paste0(
      "\nOne or more expected GSE158830 processed files are missing.\n",
      "See: ",
      file.path(
        auditdir,
        "Script26_Input_File_Audit.csv"
      )
    )
  )
}


# ============================================================
# 4. WORKBOOK INVENTORY
# ============================================================

inventory_workbook <- function(path, label) {
  
  sheets <- excel_sheets(path)
  
  tibble(
    workbook = label,
    path = path,
    sheet = sheets
  )
}


sheet_inventory <- bind_rows(
  
  inventory_workbook(
    mut_vs_wt_file,
    "MUT_vs_WT"
  ),
  
  inventory_workbook(
    mut_vs_repair_file,
    "MUT_vs_REPAIR"
  )
)


write_csv(
  sheet_inventory,
  file.path(
    auditdir,
    "GSE158830_Workbook_Sheet_Inventory.csv"
  )
)


cat("\nWorkbook sheets:\n")
print(sheet_inventory)


# ============================================================
# 5. SAFE SHEET READER
# ============================================================

safe_read_sheet <- function(path, sheet) {
  
  tryCatch(
    
    {
      x <- read_excel(
        path,
        sheet = sheet
      )
      
      x <- as.data.frame(
        x,
        stringsAsFactors = FALSE
      )
      
      x
    },
    
    error = function(e) {
      
      message(
        "Could not read: ",
        basename(path),
        " / ",
        sheet
      )
      
      data.frame()
    }
  )
}


# ============================================================
# 6. READ EVERY SHEET
# ============================================================

all_sheets <- sheet_inventory %>%
  
  mutate(
    
    data = map2(
      path,
      sheet,
      safe_read_sheet
    ),
    
    n_rows = map_int(
      data,
      nrow
    ),
    
    n_columns = map_int(
      data,
      ncol
    )
  )


sheet_summary <- all_sheets %>%
  
  select(
    workbook,
    sheet,
    n_rows,
    n_columns
  )


write_csv(
  sheet_summary,
  file.path(
    auditdir,
    "GSE158830_Sheet_Summary.csv"
  )
)


# ============================================================
# 7. COLUMN INVENTORY
# ============================================================

column_inventory <- all_sheets %>%
  
  mutate(
    
    column_names = map(
      data,
      names
    )
    
  ) %>%
  
  select(
    workbook,
    sheet,
    column_names
  ) %>%
  
  unnest_longer(
    column_names,
    values_to = "column_name"
  )


write_csv(
  column_inventory,
  file.path(
    auditdir,
    "GSE158830_Column_Inventory.csv"
  )
)


# ============================================================
# 8. GENERIC ROW SEARCH FUNCTION
#
# Search every column of every sheet.
# This avoids assuming which column contains gene names.
# ============================================================

search_dataframe <- function(
    df,
    pattern
) {
  
  if (nrow(df) == 0 || ncol(df) == 0) {
    return(
      tibble()
    )
  }
  
  char_df <- df %>%
    mutate(
      across(
        everything(),
        as.character
      )
    )
  
  hit_matrix <- sapply(
    
    char_df,
    
    function(x) {
      
      str_detect(
        str_to_lower(
          replace_na(
            x,
            ""
          )
        ),
        
        regex(
          pattern,
          ignore_case = TRUE
        )
      )
    }
  )
  
  
  if (is.null(dim(hit_matrix))) {
    
    hit_rows <- which(
      hit_matrix
    )
    
  } else {
    
    hit_rows <- which(
      apply(
        hit_matrix,
        1,
        any
      )
    )
  }
  
  
  if (length(hit_rows) == 0) {
    
    return(
      tibble()
    )
  }
  
  
  as_tibble(
    df[hit_rows, , drop = FALSE]
  ) %>%
    
    mutate(
      .source_row = hit_rows,
      .before = 1
    )
}


# ============================================================
# 9. SEARCH TARGETS ACROSS ALL SHEETS
# ============================================================

benchmark_terms <- c(
  "vigR",
  "isaA",
  "folD"
)


search_one_term <- function(term) {
  
  results <- list()
  
  counter <- 1
  
  
  for (i in seq_len(nrow(all_sheets))) {
    
    dat <- all_sheets$data[[i]]
    
    hits <- search_dataframe(
      dat,
      paste0(
        "\\b",
        term,
        "\\b"
      )
    )
    
    
    if (nrow(hits) > 0) {
      
      hits <- hits %>%
        
        mutate(
          searched_term = term,
          workbook = all_sheets$workbook[i],
          sheet = all_sheets$sheet[i],
          .before = 1
        )
      
      
      results[[counter]] <- hits
      
      counter <- counter + 1
    }
  }
  
  
  if (length(results) == 0) {
    
    return(
      tibble(
        searched_term = character(),
        workbook = character(),
        sheet = character()
      )
    )
  }
  
  
  bind_rows(results)
}


vigR_hits <- search_one_term("vigR")
isaA_hits <- search_one_term("isaA")
folD_hits <- search_one_term("folD")


# ============================================================
# 10. EXPORT RAW SEARCH HITS
# ============================================================

write_csv(
  vigR_hits,
  file.path(
    tabdir,
    "vigR_All_Local_GSE158830_Hits.csv"
  )
)


write_csv(
  isaA_hits,
  file.path(
    tabdir,
    "isaA_All_Local_GSE158830_Hits.csv"
  )
)


write_csv(
  folD_hits,
  file.path(
    tabdir,
    "folD_All_Local_GSE158830_Hits.csv"
  )
)


# ============================================================
# 11. SEARCH SUMMARY
# ============================================================

term_summary <- tibble(
  
  feature = c(
    "vigR",
    "isaA",
    "folD"
  ),
  
  local_rows_recovered = c(
    nrow(vigR_hits),
    nrow(isaA_hits),
    nrow(folD_hits)
  ),
  
  detected_in_local_processed_files = c(
    nrow(vigR_hits) > 0,
    nrow(isaA_hits) > 0,
    nrow(folD_hits) > 0
  )
)


write_csv(
  term_summary,
  file.path(
    tabdir,
    "GSE158830_Key_Feature_Recovery_Summary.csv"
  )
)


# ============================================================
# 12. PER-WORKBOOK RECOVERY
# ============================================================

make_workbook_recovery <- function(
    hits,
    feature
) {
  
  if (nrow(hits) == 0) {
    
    return(
      tibble(
        feature = feature,
        workbook = character(),
        sheet = character(),
        n_rows = integer()
      )
    )
  }
  
  
  hits %>%
    
    count(
      workbook,
      sheet,
      name = "n_rows"
    ) %>%
    
    mutate(
      feature = feature,
      .before = 1
    )
}


workbook_recovery <- bind_rows(
  
  make_workbook_recovery(
    vigR_hits,
    "vigR"
  ),
  
  make_workbook_recovery(
    isaA_hits,
    "isaA"
  ),
  
  make_workbook_recovery(
    folD_hits,
    "folD"
  )
)


write_csv(
  workbook_recovery,
  file.path(
    tabdir,
    "GSE158830_Key_Features_By_Workbook.csv"
  )
)


# ============================================================
# 13. DEFINE WHAT CAN ACTUALLY BE TESTED
#
# Published central features:
#
# A. vigR is central to the original study.
# B. isaA is a published direct regulatory target.
# C. folD was also experimentally examined as a vigR-binding RNA.
# D. vigR 3′UTR contributes to VISA / vancomycin phenotype.
#
# However:
# Our local processed files here are DESeq2 outputs.
# Therefore:
# - presence / expression effects can be benchmarked.
# - direct RNA-RNA interaction cannot be re-proven here.
# - cell-wall phenotype cannot be tested here.
# ============================================================


feature_present <- function(x) {
  
  ifelse(
    x > 0,
    "Recovered",
    "Not recovered"
  )
}


vigR_status <- feature_present(
  nrow(vigR_hits)
)

isaA_status <- feature_present(
  nrow(isaA_hits)
)

folD_status <- feature_present(
  nrow(folD_hits)
)


# ============================================================
# 14. FORMAL PUBLISHED-RESULT BENCHMARK
# ============================================================

benchmark <- tibble(
  
  benchmark_id = c(
    "GSE158830-PB01",
    "GSE158830-PB02",
    "GSE158830-PB03",
    "GSE158830-PB04",
    "GSE158830-PB05"
  ),
  
  published_feature = c(
    
    paste0(
      "vigR is a central transcript / regulatory RNA feature ",
      "of the original GSE158830 study"
    ),
    
    paste0(
      "isaA is linked to vigR regulation in the original study"
    ),
    
    paste0(
      "folD was examined as a vigR-interacting transcript ",
      "in the original study"
    ),
    
    paste0(
      "Direct vigR-isaA RNA-RNA interaction"
    ),
    
    paste0(
      "Requirement of the vigR 3′UTR for the vancomycin-",
      "intermediate resistance / cell-wall phenotype"
    )
  ),
  
  evidence_available_in_current_local_audit = c(
    
    "Processed GSE158830 DESeq2 workbooks",
    
    "Processed GSE158830 DESeq2 workbooks",
    
    "Processed GSE158830 DESeq2 workbooks",
    
    paste0(
      "No direct interaction experiment is reproduced by ",
      "the two local DESeq2 workbooks"
    ),
    
    paste0(
      "No cell-wall or antibiotic susceptibility experiment ",
      "is reproduced by the local computational files"
    )
  ),
  
  status = c(
    vigR_status,
    isaA_status,
    folD_status,
    "Not tested",
    "Not tested"
  ),
  
  quantitative_result = c(
    
    paste0(
      nrow(vigR_hits),
      " matching row(s) across local processed workbooks"
    ),
    
    paste0(
      nrow(isaA_hits),
      " matching row(s) across local processed workbooks"
    ),
    
    paste0(
      nrow(folD_hits),
      " matching row(s) across local processed workbooks"
    ),
    
    "Not applicable",
    
    "Not applicable"
  ),
  
  interpretation = c(
    
    paste0(
      "Presence of vigR-related records in the locally available ",
      "processed source data provides an internal dataset-level ",
      "benchmark, but does not by itself reproduce the full ",
      "mechanistic conclusions of the original study."
    ),
    
    paste0(
      "Recovery of isaA in the processed expression data is ",
      "consistent with its representation in the published ",
      "vigR regulatory context. Direct regulation is not inferred ",
      "from transcriptomic recovery alone."
    ),
    
    paste0(
      "Recovery of folD in the processed expression data confirms ",
      "that this published vigR-associated transcript is represented ",
      "in the local source outputs. This is not treated as independent ",
      "validation of RNA-RNA binding."
    ),
    
    paste0(
      "The original study used interaction-specific and experimental ",
      "evidence. The current DESeq2-file audit does not independently ",
      "reproduce that direct interaction."
    ),
    
    paste0(
      "The published phenotype was experimentally established in the ",
      "source study and lies outside the scope of the present ",
      "computational reproduction audit."
    )
  ),
  
  manuscript_role = c(
    "Dataset benchmark",
    "Published-feature benchmark",
    "Published-feature benchmark",
    "Not tested in current audit",
    "Not tested in current audit"
  )
)


write_csv(
  benchmark,
  file.path(
    tabdir,
    "GSE158830_Published_Result_Benchmark_FULL.csv"
  )
)


# ============================================================
# 15. MANUSCRIPT-READY VERSION
# ============================================================

manuscript_table <- benchmark %>%
  
  select(
    
    `Published feature` =
      published_feature,
    
    `Evidence examined` =
      evidence_available_in_current_local_audit,
    
    `Recovery status` =
      status,
    
    `Quantitative result` =
      quantitative_result,
    
    `Interpretation` =
      interpretation
  )


write_csv(
  manuscript_table,
  file.path(
    tabdir,
    "Supplementary_Table_GSE158830_Benchmark.csv"
  )
)


# ============================================================
# 16. CRITICAL INTERPRETATION TABLE
# ============================================================

interpretation_rules <- tibble(
  
  category = c(
    "Recovered",
    "Partially recovered",
    "Not recovered",
    "Not tested"
  ),
  
  definition = c(
    
    paste0(
      "The specified feature was detected in the exact local ",
      "data layer examined."
    ),
    
    paste0(
      "Some but not all components of the specified published ",
      "feature were recovered."
    ),
    
    paste0(
      "The specified feature was not detected in the exact local ",
      "data layer examined. This does not establish biological absence."
    ),
    
    paste0(
      "The local data layer does not directly test the published ",
      "feature or mechanism."
    )
  )
)


write_csv(
  interpretation_rules,
  file.path(
    tabdir,
    "Benchmark_Status_Definitions.csv"
  )
)


# ============================================================
# 17. EXCEL WORKBOOK
# ============================================================

xlsx_file <- file.path(
  tabdir,
  "GSE158830_vigR_Published_Result_Benchmark.xlsx"
)


wb <- createWorkbook()


addWorksheet(
  wb,
  "Benchmark"
)

writeData(
  wb,
  "Benchmark",
  benchmark
)


addWorksheet(
  wb,
  "Manuscript table"
)

writeData(
  wb,
  "Manuscript table",
  manuscript_table
)


addWorksheet(
  wb,
  "Feature summary"
)

writeData(
  wb,
  "Feature summary",
  term_summary
)


addWorksheet(
  wb,
  "Feature by workbook"
)

writeData(
  wb,
  "Feature by workbook",
  workbook_recovery
)


addWorksheet(
  wb,
  "Sheet inventory"
)

writeData(
  wb,
  "Sheet inventory",
  sheet_summary
)


addWorksheet(
  wb,
  "Status definitions"
)

writeData(
  wb,
  "Status definitions",
  interpretation_rules
)


header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE
)


for (sheet in names(wb)) {
  
  tmp <- readWorkbook(
    wb,
    sheet = sheet
  )
  
  ncols <- ncol(tmp)
  
  
  if (ncols > 0) {
    
    addStyle(
      wb,
      sheet = sheet,
      style = header_style,
      rows = 1,
      cols = 1:ncols,
      gridExpand = TRUE
    )
    
    freezePane(
      wb,
      sheet = sheet,
      firstRow = TRUE
    )
    
    setColWidths(
      wb,
      sheet = sheet,
      cols = 1:ncols,
      widths = "auto"
    )
  }
}


saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)


# ============================================================
# 18. MANUSCRIPT-SAFE METHODS TEXT
# ============================================================

methods_text <- paste0(
  
  "Benchmarking against the GSE158830 source study. ",
  
  "As an additional source-dataset benchmark, locally available ",
  "processed differential-expression outputs from GSE158830 were ",
  "audited for key features of the published vigR regulatory model. ",
  
  "The audit examined representation of vigR and the published ",
  "vigR-associated transcripts isaA and folD across the available ",
  "DESeq2 output workbooks. ",
  
  "Recovery in these processed transcriptomic files was interpreted ",
  "only as dataset-level or expression-level consistency. Direct ",
  "RNA-RNA interactions and experimentally established phenotypes ",
  "reported in the original study were classified as not tested ",
  "unless the corresponding evidence layer was directly available ",
  "and reanalysed."
)


writeLines(
  methods_text,
  file.path(
    outdir,
    "Manuscript_Safe_GSE158830_Methods.txt"
  )
)


# ============================================================
# 19. MANUSCRIPT-SAFE RESULTS TEXT
# ============================================================

results_text <- paste0(
  
  "A second analytical benchmark was performed using processed ",
  "GSE158830 outputs from the previously published vigR study. ",
  
  "The local source files contained ",
  nrow(vigR_hits),
  " row(s) matching vigR, ",
  nrow(isaA_hits),
  " row(s) matching isaA, and ",
  nrow(folD_hits),
  " row(s) matching folD. ",
  
  "These observations were treated as source-dataset consistency ",
  "checks rather than independent reproduction of the published ",
  "mechanistic model. In particular, direct vigR-mediated RNA-RNA ",
  "regulation and the experimentally established contribution of ",
  "the vigR 3′UTR to the vancomycin-intermediate phenotype were ",
  "not inferred from differential-expression tables alone."
)


writeLines(
  results_text,
  file.path(
    outdir,
    "Manuscript_Safe_GSE158830_Results.txt"
  )
)


# ============================================================
# 20. OUTPUT VERIFICATION
# ============================================================

verification <- tibble(
  
  metric = c(
    "Input workbooks",
    "Total workbook sheets",
    "vigR matching rows",
    "isaA matching rows",
    "folD matching rows",
    "Published benchmark rows",
    "Mechanistic features explicitly not tested"
  ),
  
  value = c(
    2,
    nrow(sheet_inventory),
    nrow(vigR_hits),
    nrow(isaA_hits),
    nrow(folD_hits),
    nrow(benchmark),
    sum(
      benchmark$status == "Not tested"
    )
  )
)


write_csv(
  verification,
  file.path(
    auditdir,
    "Script26_Output_Verification.csv"
  )
)


# ============================================================
# 21. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    auditdir,
    "sessionInfo.txt"
  )
)


# ============================================================
# 22. CONSOLE REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SCRIPT 26 — GSE158830 / vigR BENCHMARK COMPLETED\n")
cat("============================================================\n\n")

cat("Local processed workbooks = 2\n")
cat(
  "Workbook sheets examined =",
  nrow(sheet_inventory),
  "\n\n"
)

cat("Key published features recovered in local files:\n")
cat(
  " vigR =",
  nrow(vigR_hits),
  "matching row(s)\n"
)

cat(
  " isaA =",
  nrow(isaA_hits),
  "matching row(s)\n"
)

cat(
  " folD =",
  nrow(folD_hits),
  "matching row(s)\n\n"
)

cat("Benchmark classification:\n")

print(
  benchmark %>%
    select(
      benchmark_id,
      status,
      manuscript_role
    )
)

cat("\n")
cat("IMPORTANT:\n")
cat(
  "Direct vigR-isaA regulation and the VISA phenotype are ",
  "classified as NOT TESTED by this local DESeq2 audit.\n"
)

cat(
  "They must not be presented as independently reproduced ",
  "by the current analysis.\n\n"
)

cat("Main output:\n")
cat(
  xlsx_file,
  "\n\n"
)

cat("============================================================\n")