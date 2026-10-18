# Corrected candidate Term-seq boundary analysis

This correction retains Term-seq in the study. It changes candidate-window coordinates, saved boundary metrics, candidate interpretation and Figure 3. It does not rerun RNA-seq, CLASH or IntaRNA and does not change the JKD6008 genomic-context positions or the original 53/177-nt distances.

## Findings

- RS04120: JKD6008 811819-811847 (+) maps to JKD6009 812356-812384 (+). Both parental JKD6009 replicates have their strongest local peak at 812386 (+2 nt from the mapped candidate boundary); maxima are 849 and 1034. Both meet the original descriptive strong-local-support criteria.
- RS05585: JKD6008 1084557-1084613 (-) maps to JKD6009 1085084-1085140 (-). Boundary-window maxima are 3 in both replicates. Because the median background is zero, the original pseudocount formula gives a ratio of 4; this is weak descriptive local enrichment, not a strong sharp boundary. Dominant peaks are at 1084998 (-86 nt), with signals of 411 and 304.
- These are untreated parental JKD6009 profiles. A replicated 3-prime end is not proof of independent processing, direct regulation, expression in JKD6008, or vancomycin-dependent processing. Counts are original unnormalized track signals; peak magnitudes are not normalized comparisons between replicates.

## Verification

Each candidate with 50 nt on either side matches one unique JKD6009 locus exactly (129 and 157 nt). Mapping was fixed by reference-sequence identity before evaluating signal. All coordinates are 1-based inclusive. The same formulas and thresholds as archived Script17A were used. At the original coordinates, the new implementation reproduced all numeric metrics in the four supplied original candidate-replicate records. Corrected windows contain 804 rows (201 positions x 2 replicates x 2 candidates). SHA-256 checksums of inputs are recorded. See tables/Correction_Audit.json.

## Use the files

1. Use the separately delivered Manuscript_TermSeq_corrected.docx, which is based on the last manuscript supplied in this conversation. It replaces Figure 3 and revises affected abstract, methods, results, Table 2 note, caption, alt text, discussion, conclusion and data-availability passages. Other requested journal-formatting changes were not performed. Manuscript_TermSeq_replacements.txt provides the replacement paragraphs if you have since edited your own copy.
2. Replace the old Figure 3 with figures/Fig3_TermSeq_corrected.pdf or the 600-dpi PNG. All five panels are retained: a/b use JKD6008 genomic coordinates; c/d use mapped JKD6009 boundaries; e reports these evidence layers separately.
3. Copy project_replacements/ contents into the project root to replace the four corresponding CSV files. Back up old versions first. The original JKD6008 Candidate_Coordinates.csv must NOT be globally relabeled as JKD6009: genomic context uses JKD6008. Use the two separately named source/target coordinate tables in tables/ and their mapping file.
4. The original archived Script17A still hard-codes the old coordinates. Do not rerun it over the corrected outputs. Use scripts/correct_termseq.py as the replacement for this candidate-window step. Other old figure/integration scripts may hard-code the previous No/No interpretation; use the supplied corrected figure and integration table until those scripts are reconciled.
5. Do not continue to cite old derived copies such as the Term-seq profile table in folder 18 or any old negative boundary summaries as current results. The tables and Figure 3 in this correction supersede them. Upload the corrected tables, mapping and scripts together to GitHub and archive a new corresponding release before changing the manuscript DOI to a new version. The existing Zenodo record has not been modified.

## Reproduce this correction

With Python 3 and the existing project data in 07_processed_data/downloaded_files/GSE158830 (or the original 07_processed_data.zip):

    python scripts/correct_termseq.py --project-root /path/to/Bac-sRNA --output tables

The original 17A audit, if supplied at its original path or ZIP, is used for an explicit implementation-equivalence check. References needed for the exact-match mapping are included. There is no network request in the correction script. To redraw Figure 3 after running the correction, install matplotlib if necessary and run:

    python scripts/plot_corrected_figure.py

The corrected final evidence-integration table is supplied; its CLASH and source genomic fields are retained from your original table, while its Term-seq fields and interpretations are updated. The full study workflow has not been independently rerun.
