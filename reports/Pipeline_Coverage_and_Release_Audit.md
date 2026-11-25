# PIE-SAE pipeline coverage and release audit

25 November 2026. The comparison uses Ihor Kendiukhov’s March 2026 preprint, version 1, “Sparse autoencoders reveal organized biological knowledge but minimal regulatory logic in single-cell foundation models.” Section numbers refer to that version. [Paper](https://arxiv.org/abs/2603.02952v1) · [PIE scientific snapshot](https://github.com/vs-dx/pie-sae/tree/c70dccff03d617187812db72ddf48a9e5062876f)

## Method-by-method comparison

| Paper method | PIE evidence | Extent |
|---|---|---|
| 5.1 Models and data | One pretrained PIE checkpoint, three contexts in the selected validation panel, 128 training and 32 validation conditions. | Adapted study on one model. The Geneformer–scGPT head-to-head comparison remains specific to the paper. |
| 5.2 Activation extraction | 01–02 restore the checkpoint, verify hook identity and gene ordering, and extract 14 block outputs to HDF5. | Implemented for PIE. Encoder and processor positions correspond to latent slots, decoder positions to genes. |
| 5.3 TopK SAE | 03 fits 14 dictionaries of width 3,072 at input dimension 768 and k=32, using training-set means. | Shared design with different data, optimizer schedule and batch size. The pinned upstream implementation normalizes decoder-weight rows, whereas the paper specifies unit-norm feature columns. |
| 5.4 Feature analysis, SVD and annotation | 04 measures usage and alignment to individual top-50 SVD axes, ranks decoder genes and computes enrichment and network summaries. 06b adds a global ontology-family analysis. | Adapted. GO BP, KEGG and Reactome are ontology families. Historical 04 includes one eligible curated TRRUST term. STRING interaction summaries and directed TRRUST edges have distinct roles. |
| 5.5 Within-layer modules | 05 and 06b examine co-activation of selected cross-layer pairs. | The paper’s all-pairs **within-layer** PMI network and Leiden communities remain open. |
| 5.6 Feature patching and specificity | 07 subtracts a decoded contribution, compares random and learned directions, and transfers a coefficient between selected conditions. | Adapted intervention. The paper zeroes a latent feature in a reconstructed state and measures changes for target versus other genes. PIE’s ratio to random directions estimates a different quantity. |
| 5.7 Perturbation response mapping | PIE predicts perturbation outcomes. 07 measures sensitivity of those predictions to internal edits. | External CRISPRi target specificity using responding features and known TF targets remains open. The paper’s 3/48 result describes Geneformer. |
| 5.8 Cross-layer computational graph | 05 matches adjacent dictionaries by signed cosine and tests activation agreement. 07 examines a selected pair across decoder blocks. | Partial overlap. The multi-distance graph of activation PMI and the fraction of source features with an outgoing PMI > 3 edge remain open. |
| 5.9 Multi-tissue SAE control | Source experiments reuse fixed block-specific SAEs. | No matched K562-only versus K562 + Tabula Sapiens retraining experiment is reported. |
| 5.10 Unannotated features | Atlas cards expose historical annotation coverage and individual gene rankings. | Jaccard/Leiden grouping and functional inference from community membership remain open. |
| 5.11 Cell-type enrichment | PIE analyses stratify some comparisons by observed context. | Feature enrichment across 56 Tabula Sapiens cell types remains open. |
| 5.12 Interactive web atlas | Nine sections cover blocks, features, genes, pathways, transitions, intervention summaries, comparisons and provenance. | Implemented as a static JavaScript application with local JSON. |
| 5.13 Visualizations | PIE reports and Atlas display layer statistics, SVD comparisons and cross-layer measurements. | The paper’s module force-directed layouts and annotation TF-IDF UMAP/t-SNE visualizations remain open. |

## Release checks

- **Inputs.** The source was exporter 0.1.1 and the supplied 33-file archive, SHA-256 `7d4eb9a24c0b49d018c82a97fb77051bd47a8fc4c2b0e3cd5c112c790a705fcd`. The compact web dataset contains 36 JSON files. Provenance records the hashes of 40 input artifacts.
- **Fresh numerical audit.** A separate calculation from the saved training Top-20 catalogs and reference sets reproduced all historical feature–term keys and overlap genes. The maximum absolute q-value deviation was below 7 × 10⁻¹⁸. Per-layer historical BH gave 73,764 associations and 6,525 features. Corrected global 06b gave 73,657 and 6,492 across 54,325,248 hypotheses. Both maxT controls identified the same 19 features in the 98-feature panel. [Correction record](Case_Normalization_Correction.md)
- **CSV and card checks.** The compact layer reconstruction scores, historical decoder counts by source, and mutual-pair counts agree with the three [published CSVs](tables/). For decoder.0 feature 189, the browser card shows ASH1L at rank 1 on both splits and the historical Reactome histone-lysine methylation q-value of 1.09 × 10⁻⁶. The global 06b output reports approximately 1.063 × 10⁻⁶ for the same term. [Notebook 06b](../notebooks/06b_full_statistical_validation.ipynb)
- **Model intervention.** Re-reading the saved 07 rows reproduces the mean-change ratios to five random directions, 1.206 for DE, 1.300 for LFC and 1.358 for delta-p. Donor-directed patching improved distance in 3/9, 2/9 and 3/9 directed pairs respectively. The 189 → 413 response was measured on source-active positions in six chosen conditions. A geometry-matched control and a direct mediation experiment remain open. [Notebook 07](../notebooks/07_causal_validation.ipynb)
- **Scope.** This audit used the saved artifacts. A clean-room rerun of PIE forward passes, fitting of the 14 SAEs, and empirical biological validation remain separate work. The official validation split informed checkpoint or protocol selection, and the maxT panel was partly selected using earlier BH results.

## Reading order

The [project README](../README.md) outlines the complete PIE workflow. The [preliminary report](Preliminary_Findings.md) covers 01–05, the [full report](Research_Report.md) adds 06a–07, and the [interactive Atlas](https://vs-dx.github.io/pie-sae/atlas/) exposes the saved records at feature level.
