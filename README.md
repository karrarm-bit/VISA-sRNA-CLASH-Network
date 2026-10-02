# VISA sRNA–mRNA Regulatory Network

Reproducible analysis repository for the study of vancomycin-responsive small RNA (sRNA) regulatory networks in *Staphylococcus aureus*.

This repository contains the analysis scripts used for data preparation, transcriptomic analysis, RNA–RNA interaction analysis, evidence integration, 3′UTR/Term-seq assessment, IntaRNA structural characterization, and generation of manuscript figures and tables.

---

## Overview

The analysis integrates multiple complementary evidence layers to characterize candidate sRNA–mRNA regulatory relationships associated with vancomycin response.

The analytical framework was designed to keep transcriptomic response, physical RNA–RNA interaction evidence, transcript-boundary information, functional annotation, and RNA interaction-energy predictions as distinct evidence layers rather than combining them into a single biological-confidence score.

The repository is organized according to the analytical workflow used in the manuscript.

---

## Repository Structure

```text
VISA-sRNA-CLASH-Network/
│
├── 01_Data_Preparation/
│   └── Dataset discovery, metadata processing, experimental-design
│       curation, and preparation of processed GEO data.
│
├── 02_RNA_Seq_Analysis/
│   └── RNA-seq differential-expression and gene-annotation analyses.
│
├── 03_CLASH_Analysis/
│   └── Processing, quality assessment, benchmarking, and recovery
│       analysis of RNase III-CLASH interaction data.
│
├── 04_Evidence_Integration/
│   └── Integration and prioritization of independent evidence layers,
│       functional annotation, target resolution, and final evidence
│       integration.
│
├── 05_TermSeq_3UTR_Analysis/
│   └── Assessment of transcript boundaries and 3′UTR-associated
│       evidence using Term-seq information.
│
├── 06_IntaRNA_Analysis/
│   └── Sequence extraction and IntaRNA-based structural characterization
│       of selected RNA–RNA interaction candidates.
│
├── 07_Figures/
│   └── Scripts used to generate the manuscript figures.
│
└── 08_Tables/
    └── Scripts used to generate manuscript tables and supplementary
        table content.
