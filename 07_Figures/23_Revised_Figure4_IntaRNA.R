# ============================================================
# REVISED FIGURE 4 — CLEAN PUBLICATION VERSION
# IntaRNA-based structural characterization of
# CLASH-supported RNA–mRNA interactions
#
# Scientific logic:
#   CLASH  = physical interaction evidence
#   IntaRNA = predicted duplex energetics / localization
#
# IntaRNA is NOT treated as:
#   - experimental validation
#   - proof of direct regulation
#   - proof of causality
#   - evidence of vancomycin-specific regulation
#
# Main Figure:
#   A. Interaction-energy ranking for all 10 CLASH-supported pairs
#   B. Target-site localization for all 10 pairs
#   C. Focused structural summary of RS04120–nrdF and
#      RS05585–nrdF
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "patchwork",
  "grid"
)

missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(grid)


# ============================================================
# 2. PATHS
# ============================================================

project <- "D:/Bac-sRNA"

input_file <- file.path(
  project,
  "29_IntaRNA_structural_evidence",
  "Table5_IntaRNA_structural_characteristics.csv"
)

outdir <- file.path(
  project,
  "24_FINAL_Figure4_IntaRNA"
)

figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")

dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}


# ============================================================
# 3. READ DATA
# ============================================================

dat <- read_csv(
  input_file,
  show_col_types = FALSE
)

required_cols <- c(
  "Structural rank",
  "Regulatory RNA",
  "Target gene",
  "Interaction energy (kcal/mol)",
  "Hybridisation energy (kcal/mol)",
  "Total unfolding energy (kcal/mol)",
  "sRNA-binding region",
  "Target-binding region",
  "Seed length (nt)"
)

missing_cols <- setdiff(required_cols, names(dat))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}


# ============================================================
# 4. STANDARDIZE
# ============================================================

structural <- dat %>%
  transmute(
    structural_rank =
      as.numeric(`Structural rank`),
    
    regulatory_RNA =
      as.character(`Regulatory RNA`),
    
    target_gene =
      as.character(`Target gene`),
    
    interaction_energy =
      as.numeric(`Interaction energy (kcal/mol)`),
    
    hybridisation_energy =
      as.numeric(`Hybridisation energy (kcal/mol)`),
    
    unfolding_energy =
      as.numeric(`Total unfolding energy (kcal/mol)`),
    
    srna_region =
      as.character(`sRNA-binding region`),
    
    target_region =
      as.character(`Target-binding region`),
    
    seed_length =
      as.numeric(`Seed length (nt)`)
  ) %>%
  mutate(
    display_RNA = case_when(
      str_detect(
        regulatory_RNA,
        regex("RS04120", ignore_case = TRUE)
      ) ~ "RS04120",
      
      str_detect(
        regulatory_RNA,
        regex("RS05585", ignore_case = TRUE)
      ) ~ "RS05585",
      
      TRUE ~ regulatory_RNA
    ),
    
    pair_label =
      paste0(
        display_RNA,
        " – ",
        target_gene
      ),
    
    is_nrdF_candidate =
      display_RNA %in% c(
        "RS04120",
        "RS05585"
      ) &
      str_to_lower(target_gene) == "nrdf"
  )


# ============================================================
# 5. PARSE COORDINATE INTERVALS
# ============================================================

parse_region <- function(x) {
  
  x <- as.character(x)
  
  x <- str_replace_all(
    x,
    "[–—]",
    "-"
  )
  
  nums <- str_extract_all(
    x,
    "\\d+"
  )
  
  start <- sapply(
    nums,
    function(z) {
      if (length(z) >= 1) {
        as.numeric(z[1])
      } else {
        NA_real_
      }
    }
  )
  
  end <- sapply(
    nums,
    function(z) {
      if (length(z) >= 2) {
        as.numeric(z[2])
      } else {
        NA_real_
      }
    }
  )
  
  tibble(
    start = start,
    end = end
  )
}


srna_pos <- parse_region(
  structural$srna_region
)

target_pos <- parse_region(
  structural$target_region
)

structural <- structural %>%
  mutate(
    srna_start = srna_pos$start,
    srna_end   = srna_pos$end,
    
    target_start = target_pos$start,
    target_end   = target_pos$end,
    
    target_mid =
      (target_start + target_end) / 2
  )


# ============================================================
# 6. DATA AUDIT
# ============================================================

cat("\n=============================================\n")
cat("FIGURE 4 DATA AUDIT\n")
cat("=============================================\n")

cat("\nNumber of interactions:", nrow(structural), "\n")

cat(
  "Number of nrdF candidate interactions:",
  sum(structural$is_nrdF_candidate),
  "\n"
)

cat("\nCandidate rows:\n")

print(
  structural %>%
    filter(is_nrdF_candidate) %>%
    select(
      pair_label,
      interaction_energy,
      hybridisation_energy,
      unfolding_energy,
      srna_region,
      target_region,
      seed_length
    )
)

write_csv(
  structural,
  file.path(
    tabdir,
    "Figure4_complete_IntaRNA_data.csv"
  )
)


# ============================================================
# 7. COMMON THEME
# ============================================================

publication_theme <- theme_classic(
  base_size = 10
) +
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 11.5
      ),
    
    plot.subtitle =
      element_text(
        size = 8.5,
        margin = margin(b = 8)
      ),
    
    axis.title =
      element_text(
        size = 9.5
      ),
    
    axis.text =
      element_text(
        size = 8.5
      ),
    
    legend.title =
      element_text(
        size = 8.5,
        face = "bold"
      ),
    
    legend.text =
      element_text(
        size = 8
      ),
    
    plot.margin =
      margin(
        10, 14, 10, 10
      )
  )


# ============================================================
# 8. PANEL A
# ALL 10 INTERACTIONS — ENERGY RANKING
# ============================================================

panel_A_data <- structural %>%
  arrange(interaction_energy) %>%
  mutate(
    pair_label =
      factor(
        pair_label,
        levels = rev(pair_label)
      )
  )


panel_A <- ggplot(
  panel_A_data,
  aes(
    x = interaction_energy,
    y = pair_label
  )
) +
  
  geom_segment(
    aes(
      x = 0,
      xend = interaction_energy,
      yend = pair_label
    ),
    linewidth = 0.75
  ) +
  
  geom_point(
    aes(
      shape = is_nrdF_candidate
    ),
    size = 3.2,
    stroke = 0.9
  ) +
  
  geom_text(
    aes(
      label =
        sprintf(
          "%.2f",
          interaction_energy
        )
    ),
    hjust = 1.25,
    size = 3.0
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.3
  ) +
  
  scale_shape_manual(
    values = c(
      `FALSE` = 1,
      `TRUE`  = 16
    ),
    labels = c(
      `FALSE` = "Other CLASH-supported pair",
      `TRUE`  = "nrdF candidate interaction"
    ),
    name = NULL
  ) +
  
  labs(
    title =
      "A  Predicted interaction energies",
    
    subtitle =
      paste0(
        "IntaRNA predictions for the 10 CLASH-supported pairs; ",
        "candidate nrdF interactions are highlighted"
      ),
    
    x =
      "Predicted interaction energy (kcal/mol)",
    
    y =
      NULL
  ) +
  
  publication_theme +
  
  theme(
    axis.text.y =
      element_text(
        size = 8.3,
        face = ifelse(
          levels(panel_A_data$pair_label) %in%
            c(
              "RS04120 – nrdF",
              "RS05585 – nrdF"
            ),
          "bold",
          "plain"
        )
      ),
    
    legend.position =
      "bottom",
    
    legend.justification =
      "left"
  )


# ============================================================
# 9. PANEL B
# ALL 10 INTERACTIONS — TARGET-SITE LOCALIZATION
# ============================================================

panel_B_data <- structural %>%
  arrange(target_start) %>%
  mutate(
    pair_label =
      factor(
        pair_label,
        levels = rev(pair_label)
      )
  )


max_target <- max(
  panel_B_data$target_end,
  na.rm = TRUE
)

panel_B <- ggplot(
  panel_B_data,
  aes(
    y = pair_label
  )
) +
  
  # full target-input reference line
  geom_segment(
    aes(
      x = 0,
      xend = max_target + 12,
      yend = pair_label
    ),
    linewidth = 0.35
  ) +
  
  # predicted target-binding region
  geom_segment(
    aes(
      x = target_start,
      xend = target_end,
      yend = pair_label,
      linewidth = is_nrdF_candidate
    ),
    lineend = "butt"
  ) +
  
  geom_point(
    aes(
      x = target_mid,
      shape = is_nrdF_candidate
    ),
    size = 2.6
  ) +
  
  geom_text(
    aes(
      x = target_mid,
      label =
        paste0(
          target_start,
          "–",
          target_end
        )
    ),
    vjust = -1.15,
    size = 2.8
  ) +
  
  scale_linewidth_manual(
    values = c(
      `FALSE` = 2.6,
      `TRUE`  = 5.0
    ),
    guide = "none"
  ) +
  
  scale_shape_manual(
    values = c(
      `FALSE` = 1,
      `TRUE`  = 16
    ),
    guide = "none"
  ) +
  
  scale_x_continuous(
    limits = c(
      0,
      max_target + 15
    ),
    expand = expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  
  labs(
    title =
      "B  Predicted target-site localization",
    
    subtitle =
      paste0(
        "Binding intervals are reported relative to each ",
        "IntaRNA target input sequence"
      ),
    
    x =
      "Position in target input sequence (nt)",
    
    y =
      NULL
  ) +
  
  publication_theme +
  
  theme(
    axis.text.y =
      element_text(
        size = 8.3
      ),
    
    plot.margin =
      margin(
        10, 15, 10, 15
      )
  )


# ============================================================
# 10. PANEL C
# ONLY THE TWO nrdF CANDIDATES
# ============================================================

candidate_data <- structural %>%
  filter(is_nrdF_candidate) %>%
  arrange(structural_rank)


if (nrow(candidate_data) != 2) {
  warning(
    paste(
      "Expected exactly 2 nrdF candidate interactions, found:",
      nrow(candidate_data)
    )
  )
}


# ============================================================
# 11. PANEL C — CLEAN SUMMARY TABLE
# ============================================================

candidate_summary <- candidate_data %>%
  transmute(
    Candidate =
      display_RNA,
    
    `Interaction energy` =
      sprintf(
        "%.2f kcal/mol",
        interaction_energy
      ),
    
    `Hybridisation energy` =
      sprintf(
        "%.2f kcal/mol",
        hybridisation_energy
      ),
    
    `Unfolding energy` =
      sprintf(
        "%.2f kcal/mol",
        unfolding_energy
      ),
    
    `RNA-binding region` =
      paste0(
        srna_region,
        " nt"
      ),
    
    `Target-binding region` =
      paste0(
        target_region,
        " nt"
      ),
    
    `Predicted seed` =
      paste0(
        seed_length,
        " nt"
      )
  )


write_csv(
  candidate_summary,
  file.path(
    tabdir,
    "Figure4_nrdF_candidate_structural_summary.csv"
  )
)


candidate_long <- candidate_summary %>%
  pivot_longer(
    cols = -Candidate,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric =
      factor(
        metric,
        levels = rev(
          c(
            "Interaction energy",
            "Hybridisation energy",
            "Unfolding energy",
            "RNA-binding region",
            "Target-binding region",
            "Predicted seed"
          )
        )
      ),
    
    Candidate =
      factor(
        Candidate,
        levels = c(
          "RS05585",
          "RS04120"
        )
      )
  )


panel_C <- ggplot(
  candidate_long,
  aes(
    x = Candidate,
    y = metric
  )
) +
  
  geom_tile(
    fill = "grey96",
    linewidth = 0.45
  ) +
  
  geom_text(
    aes(
      label = value
    ),
    size = 3.6
  ) +
  
  scale_x_discrete(
    position = "top"
  ) +
  
  labs(
    title =
      "C  Structural descriptors of the two nrdF candidate interactions",
    
    subtitle =
      paste0(
        "Predicted energetic and positional properties are reported ",
        "descriptively and are not treated as functional validation"
      ),
    
    x =
      NULL,
    
    y =
      NULL
  ) +
  
  theme_minimal(
    base_size = 10
  ) +
  
  theme(
    panel.grid =
      element_blank(),
    
    axis.ticks =
      element_blank(),
    
    axis.text.x =
      element_text(
        face = "bold",
        size = 10.5
      ),
    
    axis.text.y =
      element_text(
        size = 9.2
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 11.5
      ),
    
    plot.subtitle =
      element_text(
        size = 8.5,
        margin = margin(
          b = 9
        )
      ),
    
    plot.margin =
      margin(
        10, 80, 10, 80
      )
  )


# ============================================================
# 12. SCIENTIFIC INTERPRETATION STRIP
# ============================================================

interpretation_plot <- ggplot() +
  
  annotate(
    "text",
    x = 0.5,
    y = 0.68,
    
    label =
      paste0(
        "RNase III-CLASH provides the interaction evidence; ",
        "IntaRNA provides predicted duplex energetics and ",
        "binding-site localization."
      ),
    
    fontface = "bold",
    size = 3.6
  ) +
  
  annotate(
    "text",
    x = 0.5,
    y = 0.27,
    
    label =
      paste0(
        "IntaRNA predictions are not interpreted as independent ",
        "functional validation or evidence of causality."
      ),
    
    size = 3.35
  ) +
  
  xlim(0, 1) +
  ylim(0, 1) +
  
  theme_void()


# ============================================================
# 13. COMBINE FINAL FIGURE
# ============================================================

top_row <-
  panel_A +
  panel_B +
  plot_layout(
    widths = c(
      1.0,
      1.05
    )
  )


figure4 <-
  top_row /
  panel_C /
  interpretation_plot +
  
  plot_layout(
    heights = c(
      1.65,
      1.0,
      0.28
    )
  ) +
  
  plot_annotation(
    title =
      paste0(
        "IntaRNA-based structural characterization of ",
        "CLASH-supported RNA–mRNA interactions"
      ),
    
    subtitle =
      paste0(
        "Predicted interaction energetics and binding-site ",
        "localization are evaluated as a structural evidence layer ",
        "separate from RNase III-CLASH interaction evidence."
      ),
    
    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 16
          ),
        
        plot.subtitle =
          element_text(
            size = 9.5,
            margin = margin(
              b = 10
            )
          ),
        
        plot.margin =
          margin(
            15, 20, 12, 20
          )
      )
  )


# ============================================================
# 14. DISPLAY
# ============================================================

print(figure4)


# ============================================================
# 15. EXPORT FINAL FIGURE
# ============================================================

pdf_file <- file.path(
  figdir,
  "Figure4_FINAL_IntaRNA_structural_characterization.pdf"
)

png_file <- file.path(
  figdir,
  "Figure4_FINAL_IntaRNA_structural_characterization_600dpi.png"
)

tiff_file <- file.path(
  figdir,
  "Figure4_FINAL_IntaRNA_structural_characterization_600dpi.tiff"
)


ggsave(
  filename = pdf_file,
  plot = figure4,
  width = 13.5,
  height = 10.5,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)


ggsave(
  filename = png_file,
  plot = figure4,
  width = 13.5,
  height = 10.5,
  units = "in",
  dpi = 600,
  limitsize = FALSE
)


ggsave(
  filename = tiff_file,
  plot = figure4,
  width = 13.5,
  height = 10.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  limitsize = FALSE
)


# ============================================================
# 16. SUPPLEMENTARY TABLE
# Full 10-pair structural results
# ============================================================

supplementary_table <- structural %>%
  arrange(structural_rank) %>%
  transmute(
    `Structural rank` =
      structural_rank,
    
    `Regulatory RNA` =
      display_RNA,
    
    `Target gene` =
      target_gene,
    
    `Interaction energy (kcal/mol)` =
      interaction_energy,
    
    `Hybridisation energy (kcal/mol)` =
      hybridisation_energy,
    
    `Total unfolding energy (kcal/mol)` =
      unfolding_energy,
    
    `RNA-binding region` =
      srna_region,
    
    `Target-binding region` =
      target_region,
    
    `Seed length (nt)` =
      seed_length
  )


write_csv(
  supplementary_table,
  file.path(
    tabdir,
    "Supplementary_Table_IntaRNA_all_10_pairs.csv"
  )
)


# ============================================================
# 17. FINAL FIGURE LEGEND
# ============================================================

legend_text <- paste0(
  
  "Figure 4. IntaRNA-based structural characterization of ",
  "RNase III-CLASH-supported RNA–mRNA interactions. ",
  
  "(A) Predicted interaction energies for the 10 CLASH-supported ",
  "RNA–mRNA pairs subjected to IntaRNA analysis. The RS05585–nrdF ",
  "and RS04120–nrdF interactions showed predicted interaction ",
  "energies of −13.25 and −10.82 kcal/mol, respectively. ",
  
  "(B) Predicted target-binding intervals for the same interactions, ",
  "reported relative to the corresponding target input sequences ",
  "used for IntaRNA analysis. The predicted nrdF-binding intervals ",
  "were positions 219–241 for RS05585 and 190–207 for RS04120. ",
  
  "(C) Focused structural descriptors for the two putative ",
  "3′UTR-associated RNA candidate interactions with nrdF. ",
  "Both predictions contained a 7-nt seed. ",
  
  "IntaRNA was used to characterize predicted duplex energetics ",
  "and binding-site localization of interactions already supported ",
  "by RNase III-CLASH and was not considered independent ",
  "experimental validation of direct regulation or functional ",
  "causality."
)


writeLines(
  legend_text,
  file.path(
    outdir,
    "Figure4_FINAL_Legend.txt"
  )
)


# ============================================================
# 18. MANUSCRIPT-SAFE INTERPRETATION
# ============================================================

interpretation <- c(
  
  "FIGURE 4 — MANUSCRIPT-SAFE INTERPRETATION",
  "==========================================",
  "",
  
  "RNase III-CLASH:",
  "Provides the experimental RNA-RNA interaction evidence.",
  "",
  
  "IntaRNA:",
  "Provides predicted duplex energetics and binding-site localization.",
  "",
  
  "RS05585-nrdF:",
  "Interaction energy: -13.25 kcal/mol",
  "Hybridisation energy: -19.04 kcal/mol",
  "Unfolding energy: 5.79 kcal/mol",
  "RNA-binding region: 7-27",
  "Target-binding region: 219-241",
  "Seed length: 7 nt",
  "",
  
  "RS04120-nrdF:",
  "Interaction energy: -10.82 kcal/mol",
  "Hybridisation energy: -17.91 kcal/mol",
  "Unfolding energy: 7.09 kcal/mol",
  "RNA-binding region: 2-17",
  "Target-binding region: 190-207",
  "Seed length: 7 nt",
  "",
  
  "IMPORTANT:",
  "The more negative predicted energy for RS05585-nrdF does not",
  "establish stronger biological evidence than RS04120-nrdF.",
  "",
  
  "The evidence layers remain separate:",
  "- CLASH interaction support",
  "- genomic context",
  "- target transcriptomic response",
  "- Term-seq boundary evidence",
  "- IntaRNA structural prediction",
  "",
  
  "Not established by IntaRNA:",
  "- direct functional regulation",
  "- causality",
  "- independent RNA processing",
  "- vancomycin-specific regulatory activity"
)


writeLines(
  interpretation,
  file.path(
    outdir,
    "Figure4_Interpretation_Notes.txt"
  )
)


# ============================================================
# 19. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "sessionInfo.txt"
  )
)


# ============================================================
# 20. FINAL CONSOLE REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FINAL REVISED FIGURE 4 COMPLETED\n")
cat("============================================================\n\n")

cat("Main figure structure:\n")
cat("A = Interaction-energy ranking, all 10 pairs\n")
cat("B = Target-site localization, all 10 pairs\n")
cat("C = Focused RS04120/RS05585 structural summary\n\n")

cat("Scientific interpretation:\n")
cat("CLASH = interaction evidence\n")
cat("IntaRNA = structural prediction only\n")
cat("No functional-validation claim\n")
cat("No biological-priority score\n\n")

cat("Output PNG:\n")
cat(png_file, "\n\n")

cat("Output PDF:\n")
cat(pdf_file, "\n\n")

cat("Output TIFF:\n")
cat(tiff_file, "\n\n")

cat("============================================================\n")
list.files(
  "D:/Bac-sRNA",
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.(csv|xlsx)$"
)