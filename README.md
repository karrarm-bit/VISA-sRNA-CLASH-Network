# VISA sRNA–mRNA Regulatory Network

Reproducible computational workflow for the analysis of vancomycin-responsive small RNA (sRNA) regulatory networks in *Staphylococcus aureus*.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23109159.svg)](https://doi.org/10.5281/zenodo.23109159)

## Overview

This repository contains the R scripts used for the computational analysis of the VISA sRNA–mRNA regulatory network study.

The workflow integrates complementary evidence layers to characterize candidate sRNA–mRNA regulatory relationships, including transcriptomic response, RNase III-CLASH interaction evidence, transcript-boundary/3′UTR evidence, functional annotation, and RNA–RNA structural analysis.

## Workflow

The repository is organized according to the computational workflow used in the study:

1. **Data Preparation**
   - Dataset discovery
   - Metadata processing
   - Experimental-design curation
   - Processed GEO data preparation

2. **RNA-seq Analysis**
   - Differential-expression analysis
   - Gene annotation
   - Transcriptomic response characterization

3. **CLASH Analysis**
   - RNase III-CLASH data processing
   - Interaction extraction and quality assessment
   - Benchmarking and recovery analysis

4. **Evidence Integration**
   - Integration of complementary evidence layers
   - Target identification and prioritization
   - Functional annotation

5. **Term-seq / 3′UTR Analysis**
   - Transcript-boundary assessment
   - 3′UTR evidence integration

6. **IntaRNA Analysis**
   - Sequence extraction
   - RNA–RNA interaction prediction
   - Structural characterization

7. **Figures**
   - Manuscript figure-generation scripts

8. **Tables**
   - Manuscript and supplementary table-generation scripts

## Reproducibility

The scripts are provided to document and facilitate reproduction of the computational analyses associated with the study.

The repository is versioned through GitHub releases, and the corresponding release has been archived on Zenodo.

## Data Availability

The analyses use publicly available datasets and resources referenced within the corresponding study and analysis scripts.

Raw or third-party datasets are not redistributed where their original repositories provide the authoritative source.

## Software

The workflow is implemented primarily in **R** and uses packages appropriate to the individual analysis steps.

Package requirements and computational dependencies are specified within the relevant scripts.

## Citation

If you use this repository or its computational workflow, please cite the associated study and the archived software release:

> karrarm-bit (2026). VISA sRNA–mRNA Regulatory Network — Reproducibility Release v1.0.1. Zenodo. https://doi.org/10.5281/zenodo.23109159

## Repository

GitHub:
https://github.com/karrarm-bit/VISA-sRNA-CLASH-Network

Zenodo:
https://doi.org/10.5281/zenodo.23109159

## License

This repository is distributed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.
