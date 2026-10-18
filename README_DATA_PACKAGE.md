# Essential data supplement for the VISA sRNA CLASH study

## Contents and placement

This package contains exactly 99 files: 98 selected ORIGINAL project data/result/provenance files and this README. It contains no executable analysis scripts. It accompanies the scripts already uploaded to GitHub. Original file bytes and relative paths are preserved.

Extract the CONTENTS into the GitHub repository root, beside the existing script directories. Upload the folders and their contents, not only this ZIP. Do not add a new wrapping directory around them: the scripts refer to the project paths below. Folder numbering intentionally follows the existing project data paths, not a newly invented order.

## Folder guide

| Folder | Files | Purpose |
|---|---:|---|
| 03_selected_datasets | 1 | Dataset selection provenance |
| 04_metadata | 5 | Original GEO XML metadata and dataset list |
| 05_parsed_metadata | 2 | Sample metadata and SRA accession crosswalk |
| 06_curated_design | 3 | Final experimental design, correction log and design checks |
| 07_processed_data | 17 | Original processed input files: RNA counts, CLASH, Term-seq, baseline Ribo-seq and benchmark source workbooks |
| 09_RNAseq_analysis | 12 | Merged raw counts, sample mapping, DESeq2 results, transformations and original software session |
| 10_annotation | 7 | Exact reference records, parsed annotations and gene mapping |
| 11_publication_figures_improved | 6 | Figure 1 PDF and plotting source data/session |
| 12_master_knowledgebase | 2 | Integrated gene knowledgebase and field definitions |
| 13_CLASH_network | 11 | CLASH reconstruction, endpoint mapping audits, nodes and collapsed edges |
| 14_revised_evidence_layered_network | 9 | Complete revised network, responsive subset and benchmark/candidate evidence |
| 16_revised_functional_annotation | 3 | Functional annotations and unresolved annotation audit |
| 17A_TermSeq_3UTR_boundary_audit | 6 | Candidate coordinates, strand-specific windows and replicate metrics |
| 17C_true_genomic_neighborhood | 2 | Candidate neighboring genes and nearest features |
| 17D_final_3UTR_evidence | 2 | Final saved candidate evidence and claim audit |
| 23_Revised_Figure4_IntaRNA | 4 | Saved IntaRNA figure, numeric plotting table and methodology/session |
| 25_published_result_benchmark | 1 | Published-feature recovery audit |
| 26_GSE158830_vigR_benchmark | 1 | vigR-associated published-result benchmark |
| 28_Revised_Table2 | 2 | Final Table 2 CSV and target-annotation audit |
| 29_Revised_Supplementary_Table_S1 | 2 | CLASH scoring rules and support categories |

## Scope of verification

The larger source collection was checked for consistency before this 98-file selection. The six merged raw-count columns match the six supplied original sample count files. Applying the count filter retains 2,529 of 2,783 genes. Reclassifying the saved expression results gives 221 upregulated and 235 downregulated genes. The CLASH source contains 1,356 records; the revised and collapsed networks agree on 211 source-target pairs, 207 source RNA nodes and 182 target nodes. The responsive subset contains 10 pairs and seven target nodes. Supplied S1 scoring rules reproduce all 211 saved support scores/categories. Final Table 2 matches the network. The 804 saved Term-seq window rows match the source tracks at the SAME numeric positions. These are static saved-file checks; the full R workflow and IntaRNA were NOT rerun.

## Before describing the repository as fully reproducible

- Make archived scripts portable: replace D:/Bac-sRNA with a common configurable project root; document execution order and package/IntaRNA versions. No scripts were modified or added here.
- Supply sequence FASTAs, sequence metadata, exact IntaRNA commands/settings/version and raw/combined predictions. Folder 23 holds saved figure data, not a complete executable IntaRNA analysis. The archived extraction script expects 17_network_validation, 18_reference_genome and 18_structural_interaction; Script23 expects 29_IntaRNA_structural_evidence/Table5_IntaRNA_structural_characteristics.csv. These dependencies remain absent from the supplied archives. Their exact final replacements should be reconciled with the actual code.
- Supply the original generating R code for 17C and benchmark steps 25/26, which was not identifiable in the inspected archived script release. The absence of this code is not resolved by their saved output tables.
- Identify/provide the final manuscript Figures 2 and 3 and their exact generating steps. Older alternative figure images were omitted. Figure 1 and the revised IntaRNA Figure 4 PDFs are retained as inspection outputs.
- Resolve reference-coordinate correspondence before interpreting Term-seq candidate boundaries: the input metadata and window chromosome identify JKD6009 / Saa_6009_C1, whereas candidate coordinates use JKD6008. Numerical extraction consistency does NOT demonstrate homologous positions across strains. Verify mapping and rerun affected windows if needed.
- Document apeglm shrinkage/fallback behavior used by the archived RNA-seq script and the saved log2_fold_change column. Preserve source-reference annotation provenance.
- Legacy combined-score fields in some network tables are retained for provenance; they are not the revised biological-confidence framework and should be explicitly marked as unused for revised interpretation.

This selection prioritizes inspectable inputs and analysis tables. Extra TIFF/PNG copies, Word/Excel renderings of available CSV tables, temporary files, preview-only outputs, and older intermediate figure versions are omitted to keep the file count at 99. It does not create missing experimental or computational evidence.

## Exact inventory and SHA-256

The following 98 paths cover every original project file in this package; the 99th file is this README. The SHA-256 values let a reviewer verify preserved file bytes.

| File | SHA-256 |
|---|---|
| `03_selected_datasets/Staphylococcus_sRNA_unique_GEO_shortlist.csv` | `b00650d0d3417606a434691495ae9454e0bda477d3ffa3f9757d1e28231f7f57` |
| `04_metadata/BacRegRNA_dataset_list.csv` | `c6260adf93ad7f38accabd06f737a6b66b2c4e1663f9b3d93535be48f3306ce2` |
| `04_metadata/extracted/GSE158830_family/GSE158830_family.xml` | `7b9322a627ea1d445fbb4001f87397a34e955a36f7cb594da017eae397906b72` |
| `04_metadata/extracted/GSE254530_family/GSE254530_family.xml` | `f92e6bcec429d879b4ef9978f7851c71736f53b8d12bc39f4fa5f2b42e1fa42a` |
| `04_metadata/extracted/GSE254531_family/GSE254531_family.xml` | `4c2aa49df4e8f2afca05139f439a36602b10871c888e347513bea57adf1da391` |
| `04_metadata/extracted/GSE254532_family/GSE254532_family.xml` | `85e13ab64a1e0e15c0b0828264f4a7302d82441b39ec04b01ec28077a96ec2ec` |
| `05_parsed_metadata/BacRegRNA_all_sample_metadata.csv` | `56a3123e32486fb5ea340f359e51fa5121fd3475728fa67041b45c57fe83b5c2` |
| `05_parsed_metadata/BacRegRNA_sample_SRA_table.csv` | `eea5ad4f42e69dc46ebe97491a3e0b12b4a710068b50204523ea6e082bc43dd9` |
| `06_curated_design/BacRegRNA_design_validation_checks.csv` | `fd358d2d04907052100b793c2705dbbc12ec5d2c58744575bae131449faed7d3` |
| `06_curated_design/BacRegRNA_final_curated_sample_design.csv` | `1f0929b1adffab18529a6ba57d03bc2bb7d2d44551ff5f6da9776f5ea583f088` |
| `06_curated_design/BacRegRNA_manual_corrections_log.csv` | `07c29402b86067e6338754e094290b8944f93f86b9ee5fca49776c5951763023` |
| `07_processed_data/downloaded_files/GSE158830/GSE158830_2021.03.16_MUTvsWT.DEseq2_combined_output.xlsx` | `2a7b116119c75187394e09e7ae1cab6dda849e571557f88dea103296ece2747d` |
| `07_processed_data/downloaded_files/GSE158830/GSE158830_DEseq_output_vigR_MutvsREPAIR.xlsx` | `1426b5991915afa534651cd6911ae2a2ac6a6b46495e13c72c6c4a42c6d2b887` |
| `07_processed_data/downloaded_files/GSE158830/GSE158830_MasterTable_GEO.tsv.gz` | `93ed4eb54307c2d8f2081c04ee68c64f8083b77fb2bf954af54dd0f839fc3ea9` |
| `07_processed_data/downloaded_files/GSE158830/GSM4811618_Staph_termseq_cont1_minus_strand_single_3end_base_only.gr.gz` | `cd9d1df40ad396937999926d18e080227c6fb3d53b9bddbca7f2745b0951de36` |
| `07_processed_data/downloaded_files/GSE158830/GSM4811618_Staph_termseq_cont1_plus_strand_single_3end_base_only.gr.gz` | `fc73a0f6d37897b88173f8b0699223ca026713447a4d44a80549cc0208a4bc88` |
| `07_processed_data/downloaded_files/GSE158830/GSM4811619_Staph_termseq_cont2_minus_strand_single_3end_base_only.gr.gz` | `25dcbec808866b8d3fa953d05d0b612e2575aa5776f852e2e83cd528e4db0b5c` |
| `07_processed_data/downloaded_files/GSE158830/GSM4811619_Staph_termseq_cont2_plus_strand_single_3end_base_only.gr.gz` | `f55574ecc1c1a8042afd97742ccab944650a8f37786a07818d817f22213ce884` |
| `07_processed_data/downloaded_files/GSE254530/GSE254530_RNAseq_VAN30_DEseq2.txt.gz` | `b49e4272a350145cf68758addb42abcd5cf6425f2ec51179835a8dcf0e29ce27` |
| `07_processed_data/downloaded_files/GSE254530/GSM8045612_Cont1_30mins_RNA.txt.gz` | `aa73a621373dd775a6c91db8b23090ba7febeea44b4640b3175ad6bc0d53be80` |
| `07_processed_data/downloaded_files/GSE254530/GSM8045613_Cont2_30mins_RNA.txt.gz` | `9abbf2485f106159d9b43b2fbb660b643976499ccc01222da1ffcfcf281423d3` |
| `07_processed_data/downloaded_files/GSE254530/GSM8045614_Cont3_30mins_RNA.txt.gz` | `474630b16a4aab9492de1e870f02ca11e806966afd3cd532ff1a27500bdfa099` |
| `07_processed_data/downloaded_files/GSE254530/GSM8045615_Van1_30mins_RNA.txt.gz` | `a1fa0cb17532880ccb078342515f20b45e51d03e8c99e69671266c87dfb1bc1b` |
| `07_processed_data/downloaded_files/GSE254530/GSM8045616_Van2_30mins_RNA.txt.gz` | `e40dfb4a984156f25cb07b926459c55f37f09b68ead3ddce70c0e2ada9f2ce1d` |
| `07_processed_data/downloaded_files/GSE254530/GSM8045617_Van3_30mins_RNA.txt.gz` | `0affd12bedbeff67dd6cae2ceba367aa18198b31dd9420b185dc7b704265cbb0` |
| `07_processed_data/downloaded_files/GSE254531/GSM8045618_WW_Riboseq_CONT_1.txt.gz` | `960687b209b1501b2800d39a0fb356b26d5e3e0f7ec3b3a2e612400594bc35e7` |
| `07_processed_data/downloaded_files/GSE254531/GSM8045619_WW_Riboseq_CONT_2.txt.gz` | `95f0cc562fce4467732b6d21b3500e927a81d21907fac5cb310cebbe387c5bda` |
| `07_processed_data/downloaded_files/GSE254532/GSE254532_Final_RATT_VISACLASH.txt.gz` | `a5855f355fb28650432f1662446cc6d55a4eb85bc0855fb9bda295f80a12510d` |
| `09_RNAseq_analysis/GSE254530_DESeq2_object.rds` | `33912f3520a4d530ff6e6999625318cbed9bd681ea135df6e5dea28eafec43e9` |
| `09_RNAseq_analysis/GSE254530_VST_object.rds` | `a708cd43cd5fb56eb9fed2709b7d6262b1f3af51337f4fc619ff38d7dc0b1bea` |
| `09_RNAseq_analysis/sessionInfo.txt` | `ba772b5b50dc89795bf1eae234881201a261560cc2d15dc63803a87394cbd222` |
| `09_RNAseq_analysis/tables/GSE254530_DESeq2_results_all_genes.csv` | `9c597f6d5ce44590be74458ec648f75f943b1be603c2b433009f63812af41c98` |
| `09_RNAseq_analysis/tables/GSE254530_DESeq2_significant_genes.csv` | `1dbd958d604316ae54de56da9378e72490cdcc7ab5f2eb222bfa3dbc5d6ccbe0` |
| `09_RNAseq_analysis/tables/GSE254530_DE_summary.csv` | `f25b8e7a66e197559b7974225e1e1b5493388af9eb8f8c85a03e5f2417aad2f7` |
| `09_RNAseq_analysis/tables/GSE254530_PCA_coordinates.csv` | `0fc172f6f1696a587752ec4009d2f37622aaa615b3d485b668f07e05527cd743` |
| `09_RNAseq_analysis/tables/GSE254530_VST_expression_matrix.csv` | `0a4196a413c4a208e1e54b6bf6205f127daef904ace532e32664fab5a00def69` |
| `09_RNAseq_analysis/tables/GSE254530_filter_summary.csv` | `4efbe1454119d6a465638d61ac87616c6c5e31bc8aba14d23e0eaf2348ebd302` |
| `09_RNAseq_analysis/tables/GSE254530_normalized_counts.csv` | `2765bc4987c63da5ace7cd664c9a3f152c21e9d70b4ec68d50fc5a24508ab284` |
| `09_RNAseq_analysis/tables/GSE254530_raw_count_matrix.csv` | `99db225c8b62e22dde34fee6982f04cbb8860c4a4e4a26c3f9763812bd9399f1` |
| `09_RNAseq_analysis/tables/GSE254530_sample_design.csv` | `061f470f71e0beffcb5e372062419b3e971acf0ddc24f1071188c9e2985c6134` |
| `10_annotation/reference_files/CP002120.1_JKD6008.gb` | `078c57bb690f59bbd2cde50f58fb3211f0fa03f512fdd445be755f808fef1ddf` |
| `10_annotation/reference_files/NC_017341.1_JKD6008.gb` | `627b5e49a2d4b42b82ec647af23fda34801b6e25999a1438372da0fa1a0f62b3` |
| `10_annotation/tables/GSE254530_DESeq2_annotated_all_genes.csv` | `84cbe5462b91522d8ee60e099ba8149afee43c168a6ad1392e434db25178c88d` |
| `10_annotation/tables/GSE254530_DESeq2_annotated_significant_genes.csv` | `2ab9cbd2bb97021506a6a86049ff762ef8435a550813ad15da8f5eaaad3cc4c4` |
| `10_annotation/tables/GSE254530_unmapped_gene_ids.csv` | `2361f9d556ca99b7e178a33404410409e280031fa06d114571af2d102b2c564d` |
| `10_annotation/tables/JKD6008_all_parsed_GenBank_features.csv` | `66a243bf2d8495cc1876c9a4063fce0112276816e7b9c5cecbb03ec3a3acd1ab` |
| `10_annotation/tables/JKD6008_genome_annotation.csv` | `be884c3004bf85a655f78523d9a3e94f7376b2f8ea716083ea8f8295d9babc14` |
| `11_publication_figures_improved/Figure1A_PCA_coordinates.csv` | `eb76c11e6d644d44d23f0eeb3c15345e076601b846480b7ae619e47bbba260d7` |
| `11_publication_figures_improved/Figure1C_final_displayed_labels.csv` | `4d39b0c4a478e3a7d26e311400c125598fbb717cabda6d7ca1b4314a759d7b22` |
| `11_publication_figures_improved/Figure1D_final_heatmap_genes.csv` | `208367ef07cfdb4db792a1f8af882e57fbf6c31cbcce79487c24a25735d70bbe` |
| `11_publication_figures_improved/Figure1_sample_key.csv` | `9ba737f402e86e0feb07e6e0ea893dfb389a498374d1b804ff583501249c58c4` |
| `11_publication_figures_improved/Figure_1_ABCD_improved.pdf` | `260de44b0884c663d5e8b9703a67a845c626fd2517f45a8f2ebae32bd88a5e1e` |
| `11_publication_figures_improved/sessionInfo.txt` | `d6dd821759dd78ab0718ba384389020de413253404a87c73f9062fc339f3f77e` |
| `12_master_knowledgebase/audit/master_knowledgebase_field_dictionary.csv` | `bb474056bbf6d58fd0baa8ae61045a4d67a03d6f616ebabaacc6e6c3bf58c9fe` |
| `12_master_knowledgebase/tables/JKD6008_master_gene_knowledgebase.csv` | `b2bee8adca4ad15b1ee6f977d200630fea2c50cec5ac5a51c4ed5f041113ff61` |
| `13_CLASH_network/audit/GSE254532_CLASH_summary.csv` | `22281d9fccf6edd4639f1dcde84a2de15c4a346e9b9fcdbc36ff64a2bf6abcaa` |
| `13_CLASH_network/audit/GSE254532_column_dictionary.csv` | `161c57b4c523b9f28ea02392caf1a500e3d66ae674e1c7752bd361e4444c1b26` |
| `13_CLASH_network/audit/GSE254532_endpoint_mapping_audit.csv` | `8920ddeadcd9d56bf3805d1860485bb25337886c335de2782c7ff1b17f98d256` |
| `13_CLASH_network/audit/GSE254532_unmapped_endpoints.csv` | `c1bbaff9f46d2a3ac10886535f7e017575ec4a5e7fb98e965d285db9f7f5f288` |
| `13_CLASH_network/audit/ambiguous_master_identifiers.csv` | `7961b44c0f57ab0cf96ee492c6b139f47c8f32639a084e8a7f1b54cbe857a637` |
| `13_CLASH_network/tables/GSE254532_CLASH_network_nodes.csv` | `ca52a43a6bc13412aac8b609d100d0ed86a39286369596e5b37f53ee1c2a7aeb` |
| `13_CLASH_network/tables/GSE254532_CLASH_standardized_all_interactions.csv` | `6da2ea7da0da6df77e23e75ba23d285bc47974fd6420693c25bfbca837faed20` |
| `13_CLASH_network/tables/GSE254532_regulatory_RNA_summary.csv` | `f858b287bab524b37e2375893e11f765accad09b0c1e41d557cbd7246710767f` |
| `13_CLASH_network/tables/GSE254532_sRNA_mRNA_edges_all_rows.csv` | `fd91b6ca83eb0f805e85313f52a0a2c5c15604c0c0844766dde909aa78e645e1` |
| `13_CLASH_network/tables/GSE254532_sRNA_mRNA_edges_collapsed.csv` | `1e4c13b0cac22db29b94dd546c1a616a1a87eac01c025a20a47a8ca31bb8652c` |
| `13_CLASH_network/tables/JKD6008_master_gene_knowledgebase_CLASH_updated.csv` | `2c4759ee232492d51460b8de417070f16a0c264e184fc64c87443197f659f949` |
| `14_revised_evidence_layered_network/tables/Benchmark_vs_nrdF_Descriptive_Comparison.csv` | `8281976a77261c2d21d2d7bc3170b48886980d47a2ff3fb527f9e0d9c465f427` |
| `14_revised_evidence_layered_network/tables/Evidence_Layer_Summary.csv` | `b7d3bb416d7ad3e1d616146fcb14ae2079b5b414516da7df152f150b397f4a49` |
| `14_revised_evidence_layered_network/tables/Evidence_Layered_CLASH_Network_ALL.csv` | `c5bd0b4226071d2b015ce901d6ff5a94193282cf4e8b5f9ccd1f0c24563e6f81` |
| `14_revised_evidence_layered_network/tables/Known_sRNA_Benchmark_Subset.csv` | `c70ae6f655ae814cb74b6738929a864aac92b9c48fb1136f82c789b7e57138fa` |
| `14_revised_evidence_layered_network/tables/RsaOI_CLASH_Interactions.csv` | `a92e47828bcf64a11f30c3a49111071596c5323e79d55ca95e7928de3ee5c323` |
| `14_revised_evidence_layered_network/tables/SprA2_RsaJ_CLASH_Interactions.csv` | `22e0d13980c39fefa9a190bafdcbf65a31b9923a375fb759d0176dc4c519730b` |
| `14_revised_evidence_layered_network/tables/Target_Transcriptomic_Layer_Distribution.csv` | `17f06e1ce379d6834b4721d777fbab005ca07d06770100a9a24714b687da751e` |
| `14_revised_evidence_layered_network/tables/Vancomycin_Responsive_Target_Subset.csv` | `92ef7c8ead246da27cba116d8aa359b5076fda0ee4e66ac7993f6af64dadb4b1` |
| `14_revised_evidence_layered_network/tables/nrdF_Candidate_Interactions.csv` | `83251ebe48630884c37e2ed55f67f37d5b63a78ff9f253c45e0302b1bb49a006` |
| `16_revised_functional_annotation/audit/Missing_or_Unclassified_Functional_Annotations.csv` | `78c7461fda7670cce0e86373d0a52cfd8ac62e2796e71877fe7b1972f6f30314` |
| `16_revised_functional_annotation/tables/DEG_Associated_Functional_Subset.csv` | `37ee2457cc49c23159c3b8f03d6195390fd4d059d36609f5c83d54bba636d537` |
| `16_revised_functional_annotation/tables/Evidence_Layered_Network_Functionally_Annotated.csv` | `257af55c55002854a72210f757110f727108ccba799eaab87f7ac01bf65b65b4` |
| `17A_TermSeq_3UTR_boundary_audit/Script17A_Methodological_Notes.txt` | `505d5caa73526cbc1d23dc941c70651bcc6847a27e29c1dcae0cf14d0625c2bc` |
| `17A_TermSeq_3UTR_boundary_audit/tables/Candidate_Coordinates.csv` | `5ea6302f04d3eeea0f066c15c9c6a4af7117ca49fe8129487737becf9a7d6feb` |
| `17A_TermSeq_3UTR_boundary_audit/tables/TermSeq_Candidate_Boundary_Summary.csv` | `27f80214d8f1d42b6ba694a8572a02f2667cb5628359c74ce02b6c020c896a69` |
| `17A_TermSeq_3UTR_boundary_audit/tables/TermSeq_Data_Audit.csv` | `589362f97363626e5303865bfd67f498b0c780f40a9b03f3c314396b98a360e1` |
| `17A_TermSeq_3UTR_boundary_audit/tables/TermSeq_Local_Windows.csv` | `acd053b3f387a3e3fc85c0192a63f70599d5665b70b0cd96cf648e1c07415ab9` |
| `17A_TermSeq_3UTR_boundary_audit/tables/TermSeq_Per_Replicate_Boundary_Metrics.csv` | `5ad43b39e23c793813b5b6e603e18feb3cafd78b8439140a0203d1bbab42fd47` |
| `17C_true_genomic_neighborhood/Candidate_5kb_Genomic_Neighborhood.csv` | `a5fde6836920ef1bd99a129bb61e0f124536b3a4a74830891413f14883e28d60` |
| `17C_true_genomic_neighborhood/Candidate_Closest_Features.csv` | `a9b1918e0523a0c76e7c7b04c05f4142c6b8b233f4cb691468386450c86fc98e` |
| `17D_final_3UTR_evidence/3UTR_Claim_Audit.csv` | `8e0d8a89c4b50394f27cead45298fe436654d0e08c676e6d03d3c7a2878f7312` |
| `17D_final_3UTR_evidence/FINAL_3UTR_Candidate_Evidence_Integration.csv` | `f056a6910ac33468ea28ba88f645f30c368d49617f500d5b9653992734eba65d` |
| `23_Revised_Figure4_IntaRNA/Figure4_Methodological_Notes.txt` | `c8d3efa129b58174f9ec9ab15ee5a4a7676394ea45a7943605a8102e5ea94655` |
| `23_Revised_Figure4_IntaRNA/figures/Figure4_REVISED_IntaRNA_structural_characterization.pdf` | `6ff11d389ed0f7a4e2c2802a16640db8cefd88084ba1aae63c9404cd3da3d0e5` |
| `23_Revised_Figure4_IntaRNA/sessionInfo.txt` | `2afc23dd44f9b99c78a9cc2b62534a14f1458365dc284c37d2921dd9250b67fc` |
| `23_Revised_Figure4_IntaRNA/tables/Figure4_IntaRNA_structural_data.csv` | `93f7b8ce6e51866af39f06a9ade2afae82ca60097ed928e92d4de2214b8c5b6d` |
| `25_published_result_benchmark/tables/Published_Result_Benchmark_Audit_FULL.csv` | `95327adea1befdf0bc7f589d574addf09edfe818cdfeed5ee7f7b0fbea06b4e3` |
| `26_GSE158830_vigR_benchmark/tables/GSE158830_Published_Result_Benchmark_FULL.csv` | `1ab72c6f144e6bc8eb8f1278fbebd2abf7d658b5f0006cee754f2a254b1ed244` |
| `28_Revised_Table2/Table2_Revised_CLASH_Responsive_Interactions.csv` | `b815c465b0e2a0a859a31de6aef1c19f83d505f9cfae68770634c88dd96ddf71` |
| `28_Revised_Table2/Table2_Target_Annotation_Audit.csv` | `f1946e48bda741c3b8c4667c75fb7044a406415fc71e506b5fd8c562255eaaad` |
| `29_Revised_Supplementary_Table_S1/Supplementary_Table_S1_CLASH_Scoring_Framework.csv` | `2628c2dede52e44cf6bcac89ee61a8945e0ea0dfb5f1c8d99a6dfa018daa7d04` |
| `29_Revised_Supplementary_Table_S1/Supplementary_Table_S1_CLASH_Support_Categories.csv` | `83f9ae4c4305550eaacf1f0fb155b64a2d6f46687a95fb6310d5943c55f31d7b` |
