# ============================================================
# 29_Create_Table5_and_Figure5_IntaRNA_Structural_Evidence.R
#
# Purpose:
# 1. Create manuscript-ready Table 5 summarising the structural
#    characteristics of prioritised IntaRNA interactions.
# 2. Create Supplementary Table S10 containing the complete
#    IntaRNA structural output.
# 3. Create Figure 5:
#    A. Interaction-energy ranking.
#    B. Target-binding and seed-site architecture relative to
#       the translation start site.
#
# Primary input:
# D:/Bac-sRNA/19_IntaRNA_integration_validated_coloured/tables/
# Manuscript_IntaRNA_structural_summary.csv
#
# All figure text and exported table text are in English.
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "readr",
  "stringr",
  "tidyr",
  "ggplot2",
  "openxlsx",
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


# ------------------------------------------------------------
# 2. Load packages
# ------------------------------------------------------------

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(patchwork)
library(svglite)
library(scales)


# ------------------------------------------------------------
# 3. Project folders
# ------------------------------------------------------------

project_folder <- "D:/Bac-sRNA"

input_folder <- file.path(
  project_folder,
  "19_IntaRNA_integration_validated_coloured",
  "tables"
)

input_file <- file.path(
  input_folder,
  "Manuscript_IntaRNA_structural_summary.csv"
)

output_folder <- file.path(
  project_folder,
  "29_IntaRNA_structural_evidence"
)

if (!dir.exists(output_folder)) {
  dir.create(
    output_folder,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 4. Output files
# ------------------------------------------------------------

table5_csv <- file.path(
  output_folder,
  "Table5_IntaRNA_structural_characteristics.csv"
)

table5_excel <- file.path(
  output_folder,
  "Table5_IntaRNA_structural_characteristics.xlsx"
)

supplementary_csv <- file.path(
  output_folder,
  "Supplementary_Table_S10_complete_IntaRNA_predictions.csv"
)

supplementary_excel <- file.path(
  output_folder,
  "Supplementary_Table_S10_complete_IntaRNA_predictions.xlsx"
)

section_summary_csv <- file.path(
  output_folder,
  "Section_3_6_IntaRNA_summary_statistics.csv"
)

figure5a_png <- file.path(
  output_folder,
  "Figure5A_interaction_energy_ranking.png"
)

figure5b_png <- file.path(
  output_folder,
  "Figure5B_target_binding_map.png"
)

figure5_png <- file.path(
  output_folder,
  "Figure5_combined_IntaRNA_structural_evidence.png"
)

figure5_tiff <- file.path(
  output_folder,
  "Figure5_combined_IntaRNA_structural_evidence.tiff"
)

figure5_pdf <- file.path(
  output_folder,
  "Figure5_combined_IntaRNA_structural_evidence.pdf"
)

figure5_svg <- file.path(
  output_folder,
  "Figure5_combined_IntaRNA_structural_evidence.svg"
)

plot_data_csv <- file.path(
  output_folder,
  "Figure5_plotting_data.csv"
)


# ------------------------------------------------------------
# 5. Detect input file if exact name is unavailable
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  
  possible_files <- list.files(
    path = input_folder,
    pattern = "Manuscript.*IntaRNA.*structural.*summary.*\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(possible_files) == 0) {
    
    # Search validated folders as a fallback
    possible_files <- list.files(
      path = project_folder,
      pattern = "Manuscript.*IntaRNA.*structural.*summary.*\\.csv$",
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )
  }
  
  if (length(possible_files) == 0) {
    stop(
      paste0(
        "The manuscript IntaRNA structural summary was not found.\n",
        "Expected folder:\n",
        input_folder
      )
    )
  }
  
  input_file <- possible_files[1]
  
  message(
    "Detected IntaRNA input file:\n",
    input_file
  )
}


# ------------------------------------------------------------
# 6. Helper functions
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
  
  matches <- candidates[
    candidates %in% names(data)
  ]
  
  if (
    length(matches) > 0 &&
    !new_name %in% names(data)
  ) {
    names(data)[
      names(data) == matches[1]
    ] <- new_name
  }
  
  data
}


safe_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(x)
  )
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


first_non_missing <- function(x) {
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA)
  }
  
  x[1]
}


format_scientific_text <- function(x) {
  
  ifelse(
    is.na(x),
    "",
    format(
      x,
      scientific = TRUE,
      digits = 3
    )
  )
}


# ------------------------------------------------------------
# 7. Coordinate parser
#
# Accepts formats such as:
# "45-67", "45..67", "45:67", "45;67", "-20--5"
# ------------------------------------------------------------

parse_coordinate_start <- function(x) {
  
  x <- clean_text(x)
  
  if (is.na(x)) {
    return(NA_real_)
  }
  
  numbers <- str_extract_all(
    x,
    "-?\\d+"
  )[[1]]
  
  if (length(numbers) < 1) {
    return(NA_real_)
  }
  
  as.numeric(numbers[1])
}


parse_coordinate_end <- function(x) {
  
  x <- clean_text(x)
  
  if (is.na(x)) {
    return(NA_real_)
  }
  
  numbers <- str_extract_all(
    x,
    "-?\\d+"
  )[[1]]
  
  if (length(numbers) < 2) {
    return(NA_real_)
  }
  
  as.numeric(numbers[2])
}


# ------------------------------------------------------------
# 8. Convert target-sequence coordinates to positions relative
#    to the translation start site.
#
# Extracted target sequence:
# 1     = -150 nt
# 151   = translation start site, position 0
# 251   = +100 nt
#
# Coordinates already containing negative values are retained.
# ------------------------------------------------------------

convert_to_relative_position <- function(x) {
  
  if (is.na(x)) {
    return(NA_real_)
  }
  
  if (x < 0) {
    return(x)
  }
  
  if (x >= 1 && x <= 251) {
    return(x - 151)
  }
  
  x
}


# ------------------------------------------------------------
# 9. Read input
# ------------------------------------------------------------

raw_data <- read_csv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
) |>
  clean_column_names()


# ------------------------------------------------------------
# 10. Standardise column names
# ------------------------------------------------------------

raw_data <- raw_data |>
  
  rename_first_match(
    "rank",
    c(
      "rank",
      "structural_rank",
      "energy_rank"
    )
  ) |>
  
  rename_first_match(
    "interaction_id",
    c(
      "interaction_id",
      "interaction",
      "pair_id"
    )
  ) |>
  
  rename_first_match(
    "regulatory_rna",
    c(
      "regulatory_rna",
      "srna",
      "query",
      "query_id"
    )
  ) |>
  
  rename_first_match(
    "target_gene",
    c(
      "target_gene",
      "target",
      "gene",
      "target_id"
    )
  ) |>
  
  rename_first_match(
    "interaction_energy",
    c(
      "interaction_energy_kcal_mol",
      "interaction_energy",
      "energy",
      "e"
    )
  ) |>
  
  rename_first_match(
    "hybridisation_energy",
    c(
      "hybridization_energy_kcal_mol",
      "hybridisation_energy_kcal_mol",
      "hybridization_energy",
      "hybridisation_energy",
      "hybrid_energy"
    )
  ) |>
  
  rename_first_match(
    "target_unfolding_energy",
    c(
      "target_unfolding_energy_kcal_mol",
      "target_unfolding_energy",
      "target_ed"
    )
  ) |>
  
  rename_first_match(
    "query_unfolding_energy",
    c(
      "query_unfolding_energy_kcal_mol",
      "query_unfolding_energy",
      "srna_unfolding_energy",
      "query_ed"
    )
  ) |>
  
  rename_first_match(
    "total_unfolding_energy",
    c(
      "total_unfolding_energy_kcal_mol",
      "total_unfolding_energy",
      "unfolding_energy"
    )
  ) |>
  
  rename_first_match(
    "target_binding_site",
    c(
      "target_binding_site",
      "target_binding_region",
      "target_site",
      "target_coordinates"
    )
  ) |>
  
  rename_first_match(
    "query_binding_site",
    c(
      "query_binding_site",
      "query_binding_region",
      "srna_binding_site",
      "query_site",
      "query_coordinates"
    )
  ) |>
  
  rename_first_match(
    "target_seed_site",
    c(
      "target_seed_site",
      "target_seed_region",
      "seed_target"
    )
  ) |>
  
  rename_first_match(
    "query_seed_site",
    c(
      "query_seed_site",
      "query_seed_region",
      "seed_query"
    )
  ) |>
  
  rename_first_match(
    "seed_length",
    c(
      "seed_length_nt",
      "seed_length",
      "seed"
    )
  ) |>
  
  rename_first_match(
    "target_binding_length",
    c(
      "target_binding_length_nt",
      "target_binding_length"
    )
  ) |>
  
  rename_first_match(
    "query_binding_length",
    c(
      "query_binding_length_nt",
      "query_binding_length",
      "srna_binding_length"
    )
  ) |>
  
  rename_first_match(
    "target_interaction_sequence",
    c(
      "target_interaction_sequence",
      "target_sequence",
      "target_seq"
    )
  ) |>
  
  rename_first_match(
    "query_interaction_sequence",
    c(
      "query_interaction_sequence",
      "query_sequence",
      "srna_sequence",
      "query_seq"
    )
  ) |>
  
  rename_first_match(
    "hybrid_dot_bracket",
    c(
      "hybrid_dot_bracket",
      "dot_bracket",
      "hybrid_structure"
    )
  )


# ------------------------------------------------------------
# 11. Validate required fields
# ------------------------------------------------------------

required_columns <- c(
  "regulatory_rna",
  "target_gene",
  "interaction_energy"
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
# 12. Create optional columns when absent
# ------------------------------------------------------------

optional_columns <- c(
  "rank",
  "interaction_id",
  "hybridisation_energy",
  "target_unfolding_energy",
  "query_unfolding_energy",
  "total_unfolding_energy",
  "target_binding_site",
  "query_binding_site",
  "target_seed_site",
  "query_seed_site",
  "seed_length",
  "target_binding_length",
  "query_binding_length",
  "target_interaction_sequence",
  "query_interaction_sequence",
  "hybrid_dot_bracket"
)

for (column_name in optional_columns) {
  
  if (!column_name %in% names(raw_data)) {
    raw_data[[column_name]] <- NA
  }
}


# ------------------------------------------------------------
# 13. Clean and convert values
# ------------------------------------------------------------

intarna_data <- raw_data |>
  mutate(
    
    regulatory_rna =
      clean_text(regulatory_rna),
    
    target_gene =
      clean_text(target_gene),
    
    interaction_id =
      clean_text(interaction_id),
    
    interaction_energy =
      safe_numeric(interaction_energy),
    
    hybridisation_energy =
      safe_numeric(hybridisation_energy),
    
    target_unfolding_energy =
      safe_numeric(target_unfolding_energy),
    
    query_unfolding_energy =
      safe_numeric(query_unfolding_energy),
    
    total_unfolding_energy =
      safe_numeric(total_unfolding_energy),
    
    seed_length =
      safe_numeric(seed_length),
    
    target_binding_length =
      safe_numeric(target_binding_length),
    
    query_binding_length =
      safe_numeric(query_binding_length),
    
    rank =
      safe_numeric(rank)
  ) |>
  
  filter(
    !is.na(regulatory_rna),
    !is.na(target_gene)
  )


# ------------------------------------------------------------
# 14. Derive total unfolding energy where absent
# ------------------------------------------------------------

intarna_data <- intarna_data |>
  mutate(
    
    total_unfolding_energy = case_when(
      
      !is.na(total_unfolding_energy) ~
        total_unfolding_energy,
      
      !is.na(target_unfolding_energy) &
        !is.na(query_unfolding_energy) ~
        target_unfolding_energy +
        query_unfolding_energy,
      
      TRUE ~
        NA_real_
    )
  )


# ------------------------------------------------------------
# 15. Retain one minimum-energy configuration per pair
# ------------------------------------------------------------

intarna_best <- intarna_data |>
  mutate(
    pair_label = paste0(
      regulatory_rna,
      " \u2013 ",
      target_gene
    )
  ) |>
  
  group_by(
    regulatory_rna,
    target_gene
  ) |>
  
  arrange(
    is.na(interaction_energy),
    interaction_energy,
    .by_group = TRUE
  ) |>
  
  slice_head(
    n = 1
  ) |>
  
  ungroup()


# ------------------------------------------------------------
# 16. Recalculate structural rank from interaction energy
#
# Rank 1 = most negative and most energetically favourable.
# ------------------------------------------------------------

intarna_best <- intarna_best |>
  arrange(
    interaction_energy
  ) |>
  
  mutate(
    structural_rank =
      row_number()
  )


# ------------------------------------------------------------
# 17. Parse target binding and seed coordinates
# ------------------------------------------------------------

intarna_best <- intarna_best |>
  rowwise() |>
  
  mutate(
    
    target_binding_start_raw =
      parse_coordinate_start(
        target_binding_site
      ),
    
    target_binding_end_raw =
      parse_coordinate_end(
        target_binding_site
      ),
    
    target_seed_start_raw =
      parse_coordinate_start(
        target_seed_site
      ),
    
    target_seed_end_raw =
      parse_coordinate_end(
        target_seed_site
      ),
    
    target_binding_start =
      convert_to_relative_position(
        target_binding_start_raw
      ),
    
    target_binding_end =
      convert_to_relative_position(
        target_binding_end_raw
      ),
    
    target_seed_start =
      convert_to_relative_position(
        target_seed_start_raw
      ),
    
    target_seed_end =
      convert_to_relative_position(
        target_seed_end_raw
      )
  ) |>
  
  ungroup()


# ------------------------------------------------------------
# 18. Correct reversed coordinates where necessary
# ------------------------------------------------------------

intarna_best <- intarna_best |>
  mutate(
    
    binding_min =
      pmin(
        target_binding_start,
        target_binding_end,
        na.rm = TRUE
      ),
    
    binding_max =
      pmax(
        target_binding_start,
        target_binding_end,
        na.rm = TRUE
      ),
    
    seed_min =
      pmin(
        target_seed_start,
        target_seed_end,
        na.rm = TRUE
      ),
    
    seed_max =
      pmax(
        target_seed_start,
        target_seed_end,
        na.rm = TRUE
      ),
    
    binding_min = ifelse(
      is.infinite(binding_min),
      NA_real_,
      binding_min
    ),
    
    binding_max = ifelse(
      is.infinite(binding_max),
      NA_real_,
      binding_max
    ),
    
    seed_min = ifelse(
      is.infinite(seed_min),
      NA_real_,
      seed_min
    ),
    
    seed_max = ifelse(
      is.infinite(seed_max),
      NA_real_,
      seed_max
    )
  )


# ------------------------------------------------------------
# 19. Determine structurally feasible predictions
# ------------------------------------------------------------

intarna_best <- intarna_best |>
  mutate(
    
    structurally_feasible = case_when(
      
      !is.na(interaction_energy) &
        interaction_energy < 0 ~
        "Yes",
      
      TRUE ~
        "No"
    ),
    
    near_translation_start = case_when(
      
      !is.na(binding_min) &
        !is.na(binding_max) &
        binding_max >= -50 &
        binding_min <= 50 ~
        "Yes",
      
      !is.na(binding_min) &
        !is.na(binding_max) ~
        "No",
      
      TRUE ~
        "Unavailable"
    )
  )


# ------------------------------------------------------------
# 20. Construct concise manuscript Table 5
# ------------------------------------------------------------

table5 <- intarna_best |>
  arrange(
    structural_rank
  ) |>
  
  transmute(
    
    `Structural rank` =
      structural_rank,
    
    `Regulatory RNA` =
      regulatory_rna,
    
    `Target gene` =
      target_gene,
    
    `Interaction energy (kcal/mol)` =
      round(
        interaction_energy,
        2
      ),
    
    `Hybridisation energy (kcal/mol)` =
      round(
        hybridisation_energy,
        2
      ),
    
    `Total unfolding energy (kcal/mol)` =
      round(
        total_unfolding_energy,
        2
      ),
    
    `sRNA-binding region` =
      query_binding_site,
    
    `Target-binding region` =
      target_binding_site,
    
    `Seed length (nt)` =
      seed_length
  )


# ------------------------------------------------------------
# 21. Construct complete Supplementary Table S10
# ------------------------------------------------------------

supplementary_table <- intarna_best |>
  arrange(
    structural_rank
  ) |>
  
  transmute(
    
    `Structural rank` =
      structural_rank,
    
    `Interaction ID` =
      interaction_id,
    
    `Regulatory RNA` =
      regulatory_rna,
    
    `Target gene` =
      target_gene,
    
    `Structurally feasible` =
      structurally_feasible,
    
    `Interaction energy (kcal/mol)` =
      interaction_energy,
    
    `Hybridisation energy (kcal/mol)` =
      hybridisation_energy,
    
    `Target unfolding energy (kcal/mol)` =
      target_unfolding_energy,
    
    `Query unfolding energy (kcal/mol)` =
      query_unfolding_energy,
    
    `Total unfolding energy (kcal/mol)` =
      total_unfolding_energy,
    
    `Target-binding region` =
      target_binding_site,
    
    `sRNA-binding region` =
      query_binding_site,
    
    `Target seed region` =
      target_seed_site,
    
    `sRNA seed region` =
      query_seed_site,
    
    `Seed length (nt)` =
      seed_length,
    
    `Target-binding length (nt)` =
      target_binding_length,
    
    `sRNA-binding length (nt)` =
      query_binding_length,
    
    `Target interaction sequence` =
      target_interaction_sequence,
    
    `sRNA interaction sequence` =
      query_interaction_sequence,
    
    `Hybrid dot-bracket structure` =
      hybrid_dot_bracket,
    
    `Target-binding start relative to translation start` =
      binding_min,
    
    `Target-binding end relative to translation start` =
      binding_max,
    
    `Seed start relative to translation start` =
      seed_min,
    
    `Seed end relative to translation start` =
      seed_max,
    
    `Binding region near translation start` =
      near_translation_start
  )


# ------------------------------------------------------------
# 22. Calculate Section 3.6 summary statistics
# ------------------------------------------------------------

feasible_data <- intarna_best |>
  filter(
    structurally_feasible == "Yes"
  )

n_feasible <- nrow(feasible_data)

least_negative_energy <- if (
  n_feasible > 0
) {
  max(
    feasible_data$interaction_energy,
    na.rm = TRUE
  )
} else {
  NA_real_
}

most_negative_energy <- if (
  n_feasible > 0
) {
  min(
    feasible_data$interaction_energy,
    na.rm = TRUE
  )
} else {
  NA_real_
}

median_energy <- if (
  n_feasible > 0
) {
  median(
    feasible_data$interaction_energy,
    na.rm = TRUE
  )
} else {
  NA_real_
}

strongest_pair <- if (
  n_feasible > 0
) {
  feasible_data |>
    arrange(
      interaction_energy
    ) |>
    slice_head(
      n = 1
    ) |>
    pull(
      pair_label
    )
} else {
  NA_character_
}

strongest_energy <- if (
  n_feasible > 0
) {
  feasible_data |>
    arrange(
      interaction_energy
    ) |>
    slice_head(
      n = 1
    ) |>
    pull(
      interaction_energy
    )
} else {
  NA_real_
}

minimum_seed_length <- suppressWarnings(
  min(
    intarna_best$seed_length,
    na.rm = TRUE
  )
)

maximum_seed_length <- suppressWarnings(
  max(
    intarna_best$seed_length,
    na.rm = TRUE
  )
)

if (is.infinite(minimum_seed_length)) {
  minimum_seed_length <- NA_real_
}

if (is.infinite(maximum_seed_length)) {
  maximum_seed_length <- NA_real_
}

n_near_translation_start <- sum(
  intarna_best$near_translation_start == "Yes",
  na.rm = TRUE
)


section_summary <- data.frame(
  
  Metric = c(
    "Prioritised pairs evaluated",
    "Structurally feasible predictions",
    "Least negative interaction energy (kcal/mol)",
    "Most negative interaction energy (kcal/mol)",
    "Median interaction energy (kcal/mol)",
    "Strongest predicted interaction",
    "Strongest interaction energy (kcal/mol)",
    "Minimum seed length (nt)",
    "Maximum seed length (nt)",
    "Binding regions near translation start"
  ),
  
  Value = c(
    nrow(intarna_best),
    n_feasible,
    least_negative_energy,
    most_negative_energy,
    median_energy,
    strongest_pair,
    strongest_energy,
    minimum_seed_length,
    maximum_seed_length,
    n_near_translation_start
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 23. Generate ready-to-use Section 3.6 text
# ------------------------------------------------------------

section_text <- paste0(
  "A structurally feasible prediction was obtained for ",
  n_feasible,
  " of the ",
  nrow(intarna_best),
  " prioritised pairs. For each interaction, the minimum-energy ",
  "configuration was retained. Interaction energies ranged from ",
  sprintf("%.2f", least_negative_energy),
  " to ",
  sprintf("%.2f", most_negative_energy),
  " kcal/mol, with a median of ",
  sprintf("%.2f", median_energy),
  " kcal/mol. The strongest predicted interaction was ",
  strongest_pair,
  ", with an interaction energy of ",
  sprintf("%.2f", strongest_energy),
  " kcal/mol. Seed lengths ranged from ",
  minimum_seed_length,
  " to ",
  maximum_seed_length,
  " nucleotides. Predicted binding regions overlapped or occurred ",
  "within 50 nucleotides of the translation start site for ",
  n_near_translation_start,
  " interactions."
)

section_summary <- bind_rows(
  section_summary,
  data.frame(
    Metric = "Ready-to-use Results sentence",
    Value = section_text,
    stringsAsFactors = FALSE
  )
)


# ------------------------------------------------------------
# 24. Save CSV outputs
# ------------------------------------------------------------

write_csv(
  table5,
  table5_csv,
  na = ""
)

write_csv(
  supplementary_table,
  supplementary_csv,
  na = ""
)

write_csv(
  section_summary,
  section_summary_csv,
  na = ""
)

write_csv(
  intarna_best,
  plot_data_csv,
  na = ""
)


# ============================================================
# TABLE 5 EXCEL FORMATTING
# ============================================================


# ------------------------------------------------------------
# 25. Create Table 5 workbook
# ------------------------------------------------------------

wb_table5 <- createWorkbook()

addWorksheet(
  wb_table5,
  "Table 5",
  gridLines = FALSE
)


# ------------------------------------------------------------
# 26. Table 5 styles
# ------------------------------------------------------------

title_style <- createStyle(
  fontName = "Arial",
  fontSize = 12,
  fontColour = "#FFFFFF",
  fgFill = "#1F4E78",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE
)

header_style <- createStyle(
  fontName = "Arial",
  fontSize = 10,
  fontColour = "#000000",
  fgFill = "#D9EAF7",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "TopBottomLeftRight",
  borderColour = "#7F8C8D"
)

body_style <- createStyle(
  fontName = "Arial",
  fontSize = 10,
  fontColour = "#000000",
  valign = "center",
  wrapText = TRUE,
  border = "TopBottomLeftRight",
  borderColour = "#D0D7DE"
)

centre_style <- createStyle(
  fontName = "Arial",
  fontSize = 10,
  fontColour = "#000000",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "TopBottomLeftRight",
  borderColour = "#D0D7DE"
)

rna_style <- createStyle(
  fontName = "Arial",
  fontSize = 10,
  fontColour = "#1F4E78",
  fgFill = "#EAF2F8",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "TopBottomLeftRight",
  borderColour = "#D0D7DE"
)

strongest_style <- createStyle(
  fontColour = "#006100",
  fgFill = "#C6EFCE",
  textDecoration = "bold"
)

energy_style <- createStyle(
  numFmt = "0.00"
)

integer_style <- createStyle(
  numFmt = "0"
)

note_style <- createStyle(
  fontName = "Arial",
  fontSize = 9,
  fontColour = "#404040",
  fgFill = "#F3F4F6",
  textDecoration = "italic",
  valign = "top",
  wrapText = TRUE
)


# ------------------------------------------------------------
# 27. Write Table 5 title and data
# ------------------------------------------------------------

table5_title <- paste0(
  "Table 5. Structural characteristics of the ten prioritised ",
  "regulatory RNA–mRNA interactions predicted by IntaRNA."
)

writeData(
  wb_table5,
  "Table 5",
  table5_title,
  startRow = 1,
  startCol = 1
)

mergeCells(
  wb_table5,
  "Table 5",
  cols = 1:ncol(table5),
  rows = 1
)

addStyle(
  wb_table5,
  "Table 5",
  title_style,
  rows = 1,
  cols = 1:ncol(table5),
  gridExpand = TRUE
)

writeData(
  wb_table5,
  "Table 5",
  table5,
  startRow = 3,
  startCol = 1,
  withFilter = FALSE
)

data_rows_table5 <- 4:(nrow(table5) + 3)

addStyle(
  wb_table5,
  "Table 5",
  header_style,
  rows = 3,
  cols = 1:ncol(table5),
  gridExpand = TRUE
)

addStyle(
  wb_table5,
  "Table 5",
  body_style,
  rows = data_rows_table5,
  cols = 1:ncol(table5),
  gridExpand = TRUE
)

addStyle(
  wb_table5,
  "Table 5",
  centre_style,
  rows = data_rows_table5,
  cols = 1:ncol(table5),
  gridExpand = TRUE,
  stack = TRUE
)

addStyle(
  wb_table5,
  "Table 5",
  rna_style,
  rows = data_rows_table5,
  cols = 2,
  gridExpand = TRUE,
  stack = TRUE
)


# ------------------------------------------------------------
# 28. Format numeric columns
# ------------------------------------------------------------

table5_energy_columns <- c(
  which(
    names(table5) ==
      "Interaction energy (kcal/mol)"
  ),
  which(
    names(table5) ==
      "Hybridisation energy (kcal/mol)"
  ),
  which(
    names(table5) ==
      "Total unfolding energy (kcal/mol)"
  )
)

for (column_index in table5_energy_columns) {
  
  addStyle(
    wb_table5,
    "Table 5",
    energy_style,
    rows = data_rows_table5,
    cols = column_index,
    gridExpand = TRUE,
    stack = TRUE
  )
}

rank_column <- which(
  names(table5) ==
    "Structural rank"
)

seed_column <- which(
  names(table5) ==
    "Seed length (nt)"
)

for (column_index in c(
  rank_column,
  seed_column
)) {
  
  addStyle(
    wb_table5,
    "Table 5",
    integer_style,
    rows = data_rows_table5,
    cols = column_index,
    gridExpand = TRUE,
    stack = TRUE
  )
}


# ------------------------------------------------------------
# 29. Highlight strongest interaction
# ------------------------------------------------------------

strongest_row <- which(
  table5$`Structural rank` == 1
) + 3

if (length(strongest_row) > 0) {
  
  addStyle(
    wb_table5,
    "Table 5",
    strongest_style,
    rows = strongest_row,
    cols = 1:ncol(table5),
    gridExpand = TRUE,
    stack = TRUE
  )
}


# ------------------------------------------------------------
# 30. Table widths and note
# ------------------------------------------------------------

setColWidths(
  wb_table5,
  "Table 5",
  cols = 1:ncol(table5),
  widths = c(
    14,
    29,
    14,
    22,
    24,
    25,
    22,
    22,
    16
  )
)

setRowHeights(
  wb_table5,
  "Table 5",
  rows = 1,
  heights = 36
)

setRowHeights(
  wb_table5,
  "Table 5",
  rows = 3,
  heights = 42
)

setRowHeights(
  wb_table5,
  "Table 5",
  rows = data_rows_table5,
  heights = 42
)

table5_note_row <- nrow(table5) + 5

table5_note <- paste0(
  "Note: For each regulatory RNA–mRNA pair, the prediction with ",
  "the minimum interaction energy was retained. More negative ",
  "interaction energies indicate more energetically favourable ",
  "predictions. Total unfolding energy represents the combined ",
  "energetic cost of making the interacting regions accessible. ",
  "IntaRNA predictions provide structural support and do not ",
  "constitute direct functional validation."
)

writeData(
  wb_table5,
  "Table 5",
  table5_note,
  startRow = table5_note_row,
  startCol = 1
)

mergeCells(
  wb_table5,
  "Table 5",
  cols = 1:ncol(table5),
  rows = table5_note_row
)

addStyle(
  wb_table5,
  "Table 5",
  note_style,
  rows = table5_note_row,
  cols = 1:ncol(table5),
  gridExpand = TRUE
)

setRowHeights(
  wb_table5,
  "Table 5",
  rows = table5_note_row,
  heights = 58
)

freezePane(
  wb_table5,
  "Table 5",
  firstActiveRow = 4
)

pageSetup(
  wb_table5,
  "Table 5",
  orientation = "landscape",
  paperSize = 9,
  fitToWidth = 1,
  fitToHeight = 0
)

saveWorkbook(
  wb_table5,
  table5_excel,
  overwrite = TRUE
)


# ============================================================
# SUPPLEMENTARY TABLE S10 EXCEL FORMATTING
# ============================================================


# ------------------------------------------------------------
# 31. Create Supplementary workbook
# ------------------------------------------------------------

wb_s10 <- createWorkbook()

addWorksheet(
  wb_s10,
  "Supplementary Table S10",
  gridLines = FALSE
)

writeDataTable(
  wb_s10,
  "Supplementary Table S10",
  supplementary_table,
  tableStyle = "TableStyleMedium2",
  withFilter = TRUE
)

setColWidths(
  wb_s10,
  "Supplementary Table S10",
  cols = 1:ncol(supplementary_table),
  widths = "auto"
)

# Prevent excessively wide sequence columns
sequence_columns <- which(
  names(supplementary_table) %in% c(
    "Target interaction sequence",
    "sRNA interaction sequence",
    "Hybrid dot-bracket structure"
  )
)

if (length(sequence_columns) > 0) {
  
  setColWidths(
    wb_s10,
    "Supplementary Table S10",
    cols = sequence_columns,
    widths = 50
  )
}

wrap_all_style <- createStyle(
  fontName = "Arial",
  fontSize = 9,
  valign = "top",
  wrapText = TRUE
)

addStyle(
  wb_s10,
  "Supplementary Table S10",
  wrap_all_style,
  rows = 2:(nrow(supplementary_table) + 1),
  cols = 1:ncol(supplementary_table),
  gridExpand = TRUE
)

freezePane(
  wb_s10,
  "Supplementary Table S10",
  firstActiveRow = 2,
  firstActiveCol = 4
)

saveWorkbook(
  wb_s10,
  supplementary_excel,
  overwrite = TRUE
)


# ============================================================
# FIGURE 5A: INTERACTION ENERGY RANKING
# ============================================================


# ------------------------------------------------------------
# 32. Prepare Figure 5A data
# ------------------------------------------------------------

energy_plot_data <- intarna_best |>
  arrange(
    interaction_energy
  ) |>
  
  mutate(
    
    pair_label = paste0(
      regulatory_rna,
      " \u2013 ",
      target_gene
    ),
    
    pair_label = factor(
      pair_label,
      levels = rev(pair_label)
    ),
    
    energy_label = sprintf(
      "%.2f",
      interaction_energy
    )
  )


# ------------------------------------------------------------
# 33. Create Figure 5A
# ------------------------------------------------------------

figure5a <- ggplot(
  energy_plot_data,
  aes(
    x = interaction_energy,
    y = pair_label,
    fill = interaction_energy
  )
) +
  
  geom_col(
    width = 0.68,
    colour = "white",
    linewidth = 0.45
  ) +
  
  geom_text(
    aes(
      label = energy_label
    ),
    hjust = 1.12,
    colour = "white",
    fontface = "bold",
    size = 3.7
  ) +
  
  geom_vline(
    xintercept = 0,
    colour = "#333333",
    linewidth = 0.6
  ) +
  
  scale_fill_gradient(
    low = "#9ECAE1",
    high = "#084594",
    guide = "none"
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.08,
        0.04
      )
    )
  ) +
  
  labs(
    x = "Minimum interaction energy (kcal/mol)",
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
      size = 12.5
    ),
    
    axis.text.x = element_text(
      size = 10.5,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 10.5,
      colour = "black",
      face = "italic"
    ),
    
    axis.line.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    plot.margin = margin(
      12,
      18,
      12,
      12
    )
  )


# ============================================================
# FIGURE 5B: TARGET-BINDING AND SEED ARCHITECTURE
# ============================================================


# ------------------------------------------------------------
# 34. Prepare Figure 5B data
# ------------------------------------------------------------

binding_plot_data <- intarna_best |>
  arrange(
    structural_rank
  ) |>
  
  mutate(
    
    pair_label = paste0(
      regulatory_rna,
      " \u2013 ",
      target_gene
    ),
    
    pair_label = factor(
      pair_label,
      levels = rev(pair_label)
    )
  )


# Full target region background
target_region_data <- binding_plot_data |>
  transmute(
    pair_label,
    region_start = -150,
    region_end = 100
  )


# Binding regions
binding_regions <- binding_plot_data |>
  filter(
    !is.na(binding_min),
    !is.na(binding_max)
  ) |>
  
  transmute(
    pair_label,
    region_start = binding_min,
    region_end = binding_max
  )


# Seed regions
seed_regions <- binding_plot_data |>
  filter(
    !is.na(seed_min),
    !is.na(seed_max)
  ) |>
  
  transmute(
    pair_label,
    region_start = seed_min,
    region_end = seed_max
  )


# ------------------------------------------------------------
# 35. Create Figure 5B
# ------------------------------------------------------------

figure5b <- ggplot() +
  
  geom_segment(
    data = target_region_data,
    aes(
      x = region_start,
      xend = region_end,
      y = pair_label,
      yend = pair_label
    ),
    linewidth = 2.1,
    colour = "#D9D9D9",
    lineend = "round"
  ) +
  
  geom_segment(
    data = binding_regions,
    aes(
      x = region_start,
      xend = region_end,
      y = pair_label,
      yend = pair_label
    ),
    linewidth = 6.2,
    colour = "#3182BD",
    lineend = "butt"
  ) +
  
  geom_segment(
    data = seed_regions,
    aes(
      x = region_start,
      xend = region_end,
      y = pair_label,
      yend = pair_label
    ),
    linewidth = 3.2,
    colour = "#E6550D",
    lineend = "butt"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "#252525"
  ) +
  
  annotate(
    "text",
    x = 0,
    y = length(
      levels(
        binding_plot_data$pair_label
      )
    ) + 0.75,
    label = "Translation start",
    size = 3.4,
    fontface = "bold",
    hjust = -0.05
  ) +
  
  annotate(
    "segment",
    x = -135,
    xend = -110,
    y = 0.25,
    yend = 0.25,
    linewidth = 5.5,
    colour = "#3182BD"
  ) +
  
  annotate(
    "text",
    x = -106,
    y = 0.25,
    label = "Predicted binding region",
    hjust = 0,
    size = 3.2
  ) +
  
  annotate(
    "segment",
    x = -20,
    xend = 5,
    y = 0.25,
    yend = 0.25,
    linewidth = 3,
    colour = "#E6550D"
  ) +
  
  annotate(
    "text",
    x = 9,
    y = 0.25,
    label = "Seed region",
    hjust = 0,
    size = 3.2
  ) +
  
  scale_x_continuous(
    limits = c(
      -150,
      100
    ),
    breaks = c(
      -150,
      -100,
      -50,
      0,
      50,
      100
    ),
    expand = expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +
  
  scale_y_discrete(
    drop = FALSE
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  labs(
    x = "Position relative to the translation start site (nt)",
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
      size = 12.5
    ),
    
    axis.text.x = element_text(
      size = 10.5,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 10.2,
      colour = "black",
      face = "italic"
    ),
    
    axis.line.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    plot.margin = margin(
      12,
      18,
      25,
      12
    )
  )


# ============================================================
# COMBINED FIGURE 5
# ============================================================


# ------------------------------------------------------------
# 36. Combine panels
# ------------------------------------------------------------

figure5_combined <- (
  figure5a +
    figure5b
) +
  
  plot_layout(
    widths = c(
      0.88,
      1.12
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
# 37. Display plots
# ------------------------------------------------------------

print(figure5a)
print(figure5b)
print(figure5_combined)


# ------------------------------------------------------------
# 38. Save individual panels
# ------------------------------------------------------------

ggsave(
  filename = figure5a_png,
  plot = figure5a,
  width = 8,
  height = 6.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = figure5b_png,
  plot = figure5b,
  width = 9,
  height = 6.8,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 39. Save combined PNG
# ------------------------------------------------------------

ggsave(
  filename = figure5_png,
  plot = figure5_combined,
  width = 16,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 40. Save combined TIFF
# ------------------------------------------------------------

ggsave(
  filename = figure5_tiff,
  plot = figure5_combined,
  width = 16,
  height = 7.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 41. Save combined PDF
# ------------------------------------------------------------

ggsave(
  filename = figure5_pdf,
  plot = figure5_combined,
  width = 16,
  height = 7.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# ------------------------------------------------------------
# 42. Save editable SVG
# ------------------------------------------------------------

ggsave(
  filename = figure5_svg,
  plot = figure5_combined,
  width = 16,
  height = 7.5,
  units = "in",
  device = svglite,
  bg = "white"
)


# ------------------------------------------------------------
# 43. Final report
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "Table 5 and Figure 5 created successfully\n",
  "============================================\n\n",
  
  "Input file:\n",
  normalizePath(
    input_file,
    winslash = "/"
  ),
  "\n\n",
  
  "Pairs evaluated: ",
  nrow(intarna_best),
  "\n",
  
  "Structurally feasible pairs: ",
  n_feasible,
  "\n",
  
  "Interaction-energy range: ",
  sprintf("%.2f", least_negative_energy),
  " to ",
  sprintf("%.2f", most_negative_energy),
  " kcal/mol\n",
  
  "Median interaction energy: ",
  sprintf("%.2f", median_energy),
  " kcal/mol\n",
  
  "Strongest interaction: ",
  strongest_pair,
  "\n",
  
  "Strongest energy: ",
  sprintf("%.2f", strongest_energy),
  " kcal/mol\n",
  
  "Seed-length range: ",
  minimum_seed_length,
  " to ",
  maximum_seed_length,
  " nt\n",
  
  "Binding regions near the translation start: ",
  n_near_translation_start,
  "\n\n",
  
  "Table 5 Excel:\n",
  normalizePath(
    table5_excel,
    winslash = "/"
  ),
  "\n\n",
  
  "Supplementary Table S10 Excel:\n",
  normalizePath(
    supplementary_excel,
    winslash = "/"
  ),
  "\n\n",
  
  "Combined Figure 5 PNG:\n",
  normalizePath(
    figure5_png,
    winslash = "/"
  ),
  "\n\n",
  
  "Combined Figure 5 TIFF:\n",
  normalizePath(
    figure5_tiff,
    winslash = "/"
  ),
  "\n\n",
  
  "Editable Figure 5 SVG:\n",
  normalizePath(
    figure5_svg,
    winslash = "/"
  ),
  "\n\n",
  
  "Ready-to-use Results text:\n",
  section_text,
  "\n\n"
)

print(table5)
"D:/Bac-sRNA/29_IntaRNA_structural_evidence/Figure5_combined_IntaRNA_structural_evidence.png"
# ============================================================
# 30_Refine_Figure5_IntaRNA_Panels.R
#
# Purpose:
# Refine Figure 5 generated by the previous IntaRNA script.
#
# This script must be run AFTER:
# 29_Create_Table5_and_Figure5_IntaRNA_Structural_Evidence.R
#
# Required existing object:
# intarna_best
#
# Improvements:
# - Consistent pair ordering across panels
# - Clear interaction-energy labels
# - Improved binding-region rectangles
# - Improved seed visualization
# - Non-overlapping legend
# - Clear translation-start annotation
# - Journal-ready panel spacing
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "ggplot2",
  "patchwork",
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

library(dplyr)
library(ggplot2)
library(patchwork)
library(svglite)


# ------------------------------------------------------------
# 2. Check that the previous script was executed
# ------------------------------------------------------------

if (!exists("intarna_best")) {
  stop(
    paste0(
      "Object 'intarna_best' was not found.\n",
      "Run the previous IntaRNA Table 5/Figure 5 script first, ",
      "then run this refinement script."
    )
  )
}


# ------------------------------------------------------------
# 3. Output folder
# ------------------------------------------------------------

refined_output_folder <- file.path(
  "D:/Bac-sRNA",
  "29_IntaRNA_structural_evidence",
  "Figure5_refined"
)

if (!dir.exists(refined_output_folder)) {
  dir.create(
    refined_output_folder,
    recursive = TRUE
  )
}


refined_png <- file.path(
  refined_output_folder,
  "Figure5_refined_IntaRNA_structural_evidence.png"
)

refined_tiff <- file.path(
  refined_output_folder,
  "Figure5_refined_IntaRNA_structural_evidence.tiff"
)

refined_pdf <- file.path(
  refined_output_folder,
  "Figure5_refined_IntaRNA_structural_evidence.pdf"
)

refined_svg <- file.path(
  refined_output_folder,
  "Figure5_refined_IntaRNA_structural_evidence.svg"
)

refined_panel_a_png <- file.path(
  refined_output_folder,
  "Figure5A_refined_energy_ranking.png"
)

refined_panel_b_png <- file.path(
  refined_output_folder,
  "Figure5B_refined_binding_architecture.png"
)

refined_plot_data_csv <- file.path(
  refined_output_folder,
  "Figure5_refined_plotting_data.csv"
)


# ------------------------------------------------------------
# 4. Prepare common pair labels and ordering
#
# The most energetically favourable interaction is placed
# at the top in both panels.
# ------------------------------------------------------------

refined_data <- intarna_best |>
  filter(
    !is.na(regulatory_rna),
    !is.na(target_gene),
    !is.na(interaction_energy)
  ) |>
  mutate(
    pair_label = paste0(
      regulatory_rna,
      " \u2013 ",
      target_gene
    )
  ) |>
  arrange(
    interaction_energy
  )


pair_levels <- rev(
  refined_data$pair_label
)


refined_data <- refined_data |>
  mutate(
    pair_label = factor(
      pair_label,
      levels = pair_levels
    ),
    row_number_plot = as.numeric(pair_label),
    energy_label = sprintf(
      "%.2f",
      interaction_energy
    )
  )


# ------------------------------------------------------------
# 5. Check and repair binding coordinates
# ------------------------------------------------------------

refined_data <- refined_data |>
  mutate(
    
    binding_start_plot = pmin(
      binding_min,
      binding_max,
      na.rm = TRUE
    ),
    
    binding_end_plot = pmax(
      binding_min,
      binding_max,
      na.rm = TRUE
    ),
    
    seed_start_plot = pmin(
      seed_min,
      seed_max,
      na.rm = TRUE
    ),
    
    seed_end_plot = pmax(
      seed_min,
      seed_max,
      na.rm = TRUE
    ),
    
    binding_start_plot = ifelse(
      is.infinite(binding_start_plot),
      NA_real_,
      binding_start_plot
    ),
    
    binding_end_plot = ifelse(
      is.infinite(binding_end_plot),
      NA_real_,
      binding_end_plot
    ),
    
    seed_start_plot = ifelse(
      is.infinite(seed_start_plot),
      NA_real_,
      seed_start_plot
    ),
    
    seed_end_plot = ifelse(
      is.infinite(seed_end_plot),
      NA_real_,
      seed_end_plot
    )
  )


# ------------------------------------------------------------
# 6. Prevent extremely narrow regions from disappearing
#
# This changes only the visual width, not the exported
# biological coordinates.
# ------------------------------------------------------------

minimum_binding_display_width <- 5
minimum_seed_display_width <- 2


refined_data <- refined_data |>
  mutate(
    
    binding_midpoint = (
      binding_start_plot +
        binding_end_plot
    ) / 2,
    
    binding_width = abs(
      binding_end_plot -
        binding_start_plot
    ),
    
    binding_display_start = case_when(
      is.na(binding_midpoint) ~ NA_real_,
      binding_width >= minimum_binding_display_width ~
        binding_start_plot,
      TRUE ~
        binding_midpoint -
        minimum_binding_display_width / 2
    ),
    
    binding_display_end = case_when(
      is.na(binding_midpoint) ~ NA_real_,
      binding_width >= minimum_binding_display_width ~
        binding_end_plot,
      TRUE ~
        binding_midpoint +
        minimum_binding_display_width / 2
    ),
    
    seed_midpoint = (
      seed_start_plot +
        seed_end_plot
    ) / 2,
    
    seed_width = abs(
      seed_end_plot -
        seed_start_plot
    ),
    
    seed_display_start = case_when(
      is.na(seed_midpoint) ~ NA_real_,
      seed_width >= minimum_seed_display_width ~
        seed_start_plot,
      TRUE ~
        seed_midpoint -
        minimum_seed_display_width / 2
    ),
    
    seed_display_end = case_when(
      is.na(seed_midpoint) ~ NA_real_,
      seed_width >= minimum_seed_display_width ~
        seed_end_plot,
      TRUE ~
        seed_midpoint +
        minimum_seed_display_width / 2
    )
  )


# ------------------------------------------------------------
# 7. Save refined plotting data
# ------------------------------------------------------------

readr::write_csv(
  refined_data,
  refined_plot_data_csv,
  na = ""
)


# ============================================================
# PANEL A
# Interaction-energy ranking
# ============================================================


# ------------------------------------------------------------
# 8. Determine Panel A limits
# ------------------------------------------------------------

minimum_energy <- min(
  refined_data$interaction_energy,
  na.rm = TRUE
)

energy_left_limit <- floor(
  minimum_energy - 1.5
)


# ------------------------------------------------------------
# 9. Create refined Panel A
# ------------------------------------------------------------

figure5a_refined <- ggplot(
  refined_data,
  aes(
    x = interaction_energy,
    y = pair_label,
    fill = interaction_energy
  )
) +
  
  geom_col(
    width = 0.66,
    colour = "#FFFFFF",
    linewidth = 0.45
  ) +
  
  geom_text(
    aes(
      x = interaction_energy - 0.22,
      label = energy_label
    ),
    hjust = 1,
    colour = "#FFFFFF",
    fontface = "bold",
    size = 3.7
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.65,
    colour = "#303030"
  ) +
  
  scale_fill_gradient(
    low = "#A6CEE3",
    high = "#174A97",
    guide = "none"
  ) +
  
  scale_x_continuous(
    limits = c(
      energy_left_limit,
      0.6
    ),
    breaks = pretty(
      c(
        energy_left_limit,
        0
      ),
      n = 5
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  labs(
    x = "Minimum interaction energy (kcal/mol)",
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
      size = 13,
      margin = margin(
        t = 10
      )
    ),
    
    axis.text.x = element_text(
      size = 10.5,
      colour = "#000000"
    ),
    
    axis.text.y = element_text(
      size = 10.5,
      colour = "#000000",
      face = "italic",
      margin = margin(
        r = 8
      )
    ),
    
    axis.line.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    axis.line.x = element_line(
      linewidth = 0.7
    ),
    
    panel.grid = element_blank(),
    
    plot.margin = margin(
      t = 15,
      r = 15,
      b = 12,
      l = 12
    )
  )


# ============================================================
# PANEL B
# Binding-site architecture
# ============================================================


# ------------------------------------------------------------
# 10. Prepare background target-region data
# ------------------------------------------------------------

background_regions <- refined_data |>
  transmute(
    pair_label,
    row_number_plot,
    xmin = -150,
    xmax = 100,
    ymin = row_number_plot - 0.055,
    ymax = row_number_plot + 0.055
  )


binding_rectangles <- refined_data |>
  filter(
    !is.na(binding_display_start),
    !is.na(binding_display_end)
  ) |>
  transmute(
    pair_label,
    row_number_plot,
    xmin = binding_display_start,
    xmax = binding_display_end,
    ymin = row_number_plot - 0.145,
    ymax = row_number_plot + 0.145,
    Region = "Predicted binding region"
  )


seed_rectangles <- refined_data |>
  filter(
    !is.na(seed_display_start),
    !is.na(seed_display_end)
  ) |>
  transmute(
    pair_label,
    row_number_plot,
    xmin = seed_display_start,
    xmax = seed_display_end,
    ymin = row_number_plot - 0.070,
    ymax = row_number_plot + 0.070,
    Region = "Seed region"
  )


# ------------------------------------------------------------
# 11. Create refined Panel B
# ------------------------------------------------------------

figure5b_refined <- ggplot() +
  
  geom_rect(
    data = background_regions,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "#D9D9D9",
    colour = NA
  ) +
  
  geom_rect(
    data = binding_rectangles,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = Region
    ),
    colour = "#1F5F8B",
    linewidth = 0.30
  ) +
  
  geom_rect(
    data = seed_rectangles,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = Region
    ),
    colour = "#B23B00",
    linewidth = 0.20
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.70,
    colour = "#252525"
  ) +
  
  annotate(
    "label",
    x = 0,
    y = nrow(refined_data) + 0.68,
    label = "Translation start",
    fontface = "bold",
    size = 3.5,
    label.size = 0,
    fill = "#FFFFFF",
    hjust = 0.5
  ) +
  
  scale_fill_manual(
    name = NULL,
    values = c(
      "Predicted binding region" = "#3182BD",
      "Seed region" = "#E6550D"
    ),
    breaks = c(
      "Predicted binding region",
      "Seed region"
    )
  ) +
  
  scale_x_continuous(
    limits = c(
      -150,
      100
    ),
    breaks = c(
      -150,
      -100,
      -50,
      0,
      50,
      100
    ),
    expand = expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +
  
  scale_y_continuous(
    limits = c(
      0.5,
      nrow(refined_data) + 1.0
    ),
    breaks = refined_data$row_number_plot,
    labels = as.character(
      refined_data$pair_label
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  labs(
    x = "Position relative to the translation start site (nt)",
    y = NULL
  ) +
  
  guides(
    fill = guide_legend(
      direction = "horizontal",
      title.position = "top",
      label.position = "right",
      keywidth = grid::unit(
        1.25,
        "cm"
      ),
      keyheight = grid::unit(
        0.35,
        "cm"
      ),
      byrow = TRUE
    )
  ) +
  
  coord_cartesian(
    clip = "off"
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
      size = 13,
      margin = margin(
        t = 10
      )
    ),
    
    axis.text.x = element_text(
      size = 10.5,
      colour = "#000000"
    ),
    
    axis.text.y = element_text(
      size = 10.2,
      colour = "#000000",
      face = "italic",
      margin = margin(
        r = 8
      )
    ),
    
    axis.line.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    axis.line.x = element_line(
      linewidth = 0.7
    ),
    
    panel.grid = element_blank(),
    
    legend.position = "top",
    
    legend.justification = "center",
    
    legend.box.just = "center",
    
    legend.margin = margin(
      t = 0,
      r = 0,
      b = 8,
      l = 0
    ),
    
    legend.text = element_text(
      size = 10.5
    ),
    
    plot.margin = margin(
      t = 5,
      r = 18,
      b = 12,
      l = 12
    )
  )


# ============================================================
# COMBINED FIGURE
# ============================================================


# ------------------------------------------------------------
# 12. Combine panels
# ------------------------------------------------------------

figure5_refined <- (
  figure5a_refined |
    figure5b_refined
) +
  
  plot_layout(
    widths = c(
      0.88,
      1.22
    )
  ) +
  
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        family = "Arial",
        face = "bold",
        size = 18,
        colour = "#000000"
      )
    )
  )


# ------------------------------------------------------------
# 13. Display refined figures
# ------------------------------------------------------------

print(figure5a_refined)
print(figure5b_refined)
print(figure5_refined)


# ------------------------------------------------------------
# 14. Save separate Panel A
# ------------------------------------------------------------

ggsave(
  filename = refined_panel_a_png,
  plot = figure5a_refined,
  width = 7.5,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 15. Save separate Panel B
# ------------------------------------------------------------

ggsave(
  filename = refined_panel_b_png,
  plot = figure5b_refined,
  width = 9.5,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 16. Save combined PNG
# ------------------------------------------------------------

ggsave(
  filename = refined_png,
  plot = figure5_refined,
  width = 16.5,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 17. Save combined TIFF
# ------------------------------------------------------------

ggsave(
  filename = refined_tiff,
  plot = figure5_refined,
  width = 16.5,
  height = 7.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 18. Save vector PDF
# ------------------------------------------------------------

ggsave(
  filename = refined_pdf,
  plot = figure5_refined,
  width = 16.5,
  height = 7.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# ------------------------------------------------------------
# 19. Save editable SVG
# ------------------------------------------------------------

ggsave(
  filename = refined_svg,
  plot = figure5_refined,
  width = 16.5,
  height = 7.5,
  units = "in",
  device = svglite,
  bg = "white"
)


# ------------------------------------------------------------
# 20. Final report
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "Refined Figure 5 created successfully\n",
  "============================================\n\n",
  
  "Pairs displayed: ",
  nrow(refined_data),
  "\n\n",
  
  "Combined PNG:\n",
  normalizePath(
    refined_png,
    winslash = "/"
  ),
  "\n\n",
  
  "Combined TIFF:\n",
  normalizePath(
    refined_tiff,
    winslash = "/"
  ),
  "\n\n",
  
  "Vector PDF:\n",
  normalizePath(
    refined_pdf,
    winslash = "/"
  ),
  "\n\n",
  
  "Editable SVG:\n",
  normalizePath(
    refined_svg,
    winslash = "/"
  ),
  "\n\n",
  
  "Refined plotting data:\n",
  normalizePath(
    refined_plot_data_csv,
    winslash = "/"
  ),
  "\n\n"
)