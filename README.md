# VISA sRNA–mRNA Regulatory Network

Analysis scripts, processed results, and supporting files for an evidence-layered reanalysis of regulatory RNA–mRNA interactions in vancomycin-intermediate *Staphylococcus aureus*.

## Authors

1. **Karrar Mahmood Shaker Al-Khfaji** — corresponding author
2. **Sarah Mohammed Mohsin**
3. **Zeina Haider Abbas**

Department of Clinical and Laboratory Sciences, College of Pharmacy, University of Kerbala, Karbala, Iraq.

Correspondence: Karrar.m@uokerbala.edu.iq

## Important: corrected Term-seq analysis

Use **`scripts/correct_termseq.py`** for the corrected candidate-boundary analysis and **`scripts/plot_corrected_figure.py`** to regenerate Figure 3. The archived R boundary and figure scripts document the earlier analysis; they have not been reconciled with this correction and can regenerate superseded results.

The correction maps candidate loci defined on JKD6008 to the homologous JKD6009 reference sequence before extracting Term-seq signal. The Term-seq samples GSM4811618 and GSM4811619 are untreated parental-strain JKD6009 cultures.

| Candidate | JKD6008 coordinates | Mapped JKD6009 coordinates | Corrected Term-seq finding |
|---|---|---|---|
| RS04120 | 811819–811847 (+) | 812356–812384 (+) | Strong local peak at 812386, 2 nt from the mapped candidate 3′ boundary, in both replicates |
| RS05585 | 1084557–1084613 (−) | 1085084–1085140 (−) | Weak boundary-proximal signal; dominant local peak at 1084998, 86 nt from the mapped candidate 3′ boundary |

Coordinates are 1-based inclusive. Local parental-strain 3′-end support does not establish independent RNA processing, direct regulation, or a vancomycin-dependent effect in JKD6008. Original track signals are not normalized comparisons between replicates.

## Study and datasets

The workflow keeps transcriptomic response, physically captured RNA–RNA interactions, genomic context, transcript-boundary evidence, and predicted interaction energies as separate evidence layers.

| GEO series | Role |
|---|---|
| GSE254530 | Primary RNA-seq differential-expression analysis |
| GSE254532 | RNase III-CLASH interaction-network reconstruction |
| GSE254531 | Untreated baseline ribosome-profiling context |
| GSE158830 | Supporting annotation and parental-strain Term-seq evidence |

The original deposited datasets remain the authoritative source of third-party data. This repository includes selected processed results; it is not a complete redistribution of all source datasets.

## Repository organization

| Path | Contents |
|---|---|
| `01_Data_Preparation/` | Dataset discovery, metadata curation, and download scripts |
| `02_RNA_Seq_Analysis/` | Differential-expression and annotation scripts |
| `03_CLASH_Analysis/` | CLASH processing and benchmark scripts |
| `04_Evidence_Integration/` | Evidence integration and annotation scripts |
| `05_TermSeq_3UTR_Analysis/` | Earlier R boundary-analysis scripts; see correction notice above |
| `06_IntaRNA_Analysis/` | Sequence preparation and structural-analysis scripts |
| `07_Figures/`, `08_Tables/` | Original figure and table scripts |
| `scripts/` | Corrected Term-seq extraction and Figure 3 plotting scripts |
| `references/`, `inputs/` | Reference sequences, provenance, and genomic-neighborhood input |
| `tables/` | Corrected coordinate mapping, boundary windows, metrics, and integration table |
| `figures/` | Corrected Figure 3 as PDF and PNG |
| `17A_TermSeq_3UTR_boundary_audit/tables/` | Corrected boundary windows, per-replicate metrics, and candidate summary; JKD6008 source-coordinate records remain separately identified |
| `17D_final_3UTR_evidence/` | Corrected final evidence-integration table and earlier supporting audit |
| `project_replacements/` | Copies of the four corrected tables for replacement in the original project layout |
| `README_CORRECTION.md` | Detailed correction methods, interpretation, and reproduction instructions |

Other numbered data directories contain selected processed analysis outputs and supporting records. Earlier notes and audit files document their original analysis stage and should be read with the correction notice.

## Reproduce the corrected Term-seq step

Run from the repository root using Python 3. The extraction script uses the Python standard library. Supply the original four strand-specific Term-seq tracks in `07_processed_data/downloaded_files/GSE158830/`, or provide the original `07_processed_data.zip` in the project root.

```bash
python scripts/correct_termseq.py --project-root /path/to/Bac-sRNA --output tables
```

Replace `/path/to/Bac-sRNA` with the directory containing those source data. The script uses the included reference sequences and performs no network requests. If the original 17A audit tables are available at their expected project location, it also compares the old-coordinate calculations against those records.

To regenerate Figure 3 from the supplied corrected tables, install Matplotlib and run:

```bash
python -m pip install matplotlib
python scripts/plot_corrected_figure.py
```

The figure is written to `figures/`. The final evidence-integration table is supplied with the correction; these two commands do not rerun the complete R integration workflow.

## Reproducibility scope

The original analysis was implemented primarily in R; individual scripts specify package requirements and may require changing local project paths such as `D:/Bac-sRNA`. Saved session information is included with selected outputs.

The correction verified unique exact sequence matches across each candidate and 50 nt on either side. At the old coordinates, the correction implementation reproduced the supplied original numerical boundary metrics. Corrected windows contain 804 rows: 201 positions × two replicates × two candidates. Input checksums and an audit record are supplied in `tables/`.

Only the candidate Term-seq coordinate/window step was recomputed for this correction. RNA-seq, CLASH, and IntaRNA were not rerun as part of it. This verification does not establish that the entire study workflow has been independently reproduced.

## Citation

Please include all three authors when citing this repository:

> Al-Khfaji, Karrar Mahmood Shaker; Mohsin, Sarah Mohammed; Abbas, Zeina Haider. (2026). *VISA sRNA–mRNA Regulatory Network* [Analysis code and processed results]. GitHub. https://github.com/karrarm-bit/VISA-sRNA-CLASH-Network

For a version-specific citation, add the release version and its corresponding Zenodo DOI after that release has been archived. The earlier archive is available at https://doi.org/10.5281/zenodo.23109159; do not assume that this earlier DOI identifies the corrected release.

## Author contributions

**Karrar Mahmood Shaker Al-Khfaji:** Conceptualization, Methodology, Software, Formal analysis, Data curation, Investigation, Visualization, Project administration, Writing – original draft, and Writing – review & editing.

**Sarah Mohammed Mohsin:** Writing – original draft and Writing – review & editing.

**Zeina Haider Abbas:** Writing – original draft and Writing – review & editing.

## License

Repository content is distributed under the Creative Commons Attribution 4.0 International license (CC BY 4.0), as declared for this project. Third-party datasets and reference resources remain subject to their original terms.
