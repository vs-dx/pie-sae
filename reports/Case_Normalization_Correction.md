# PIE-SAE: Case Normalization Correction and Reproducibility Record

**Execution status.** The corrected 06a and 06b notebooks ran locally on 21 November 2026. Commit [`eea7171`](https://github.com/vs-dx/pie-sae/commit/eea71712a62cf2e5b9061fcf132e5b730cc6e52c) records the executed notebooks, with 33 and 41 sequential code cells respectively. The original results are available at [`0812a7c`](https://github.com/vs-dx/pie-sae/commit/0812a7c0a7efb6926784138b5cf39580b3008dee).

**Basis:** Independent critical audit dated 21 November 2026, supplied as `PIE_SAE_Audit_RU.pdf` and `pie_sae_audit.py`. The independent audit used saved Top-20 catalogs, external reference sets and selected official split files. It did **not** independently rerun PIE forward passes or train SAE checkpoints.

## Root cause

The original annotation function in `src/sae_biological_annotation.py` converts reference-set and PIE gene symbols to uppercase before comparing gene identities. In notebook 06a/06b the reference members were already uppercase, while `set(vocab["genes"])` retained mixed-case PIE gene names. Exact string matching then lost valid memberships. This was a gene-identity mismatch, not evidence that the historical 04 annotations used additional unknown biology.

The audit identified 70 mixed-case symbols among the 6,642 output genes, with no collisions after uppercase normalization. Twelve affected symbols restored relevant memberships; 49 ontology memberships and a previously filtered-out term were recovered. The GO BP boundary example is `GO:0010648`: `C6orf89` must match `C6ORF89` before the 5-gene minimum-size test.

The original `04` ontology family is **4,421 GO BP/KEGG/Reactome terms plus one curated TRRUST TF term = 4,422 terms**. The historical explicit overlap-test residual was incorrectly described as 337; its component totals sum to **340**. Corrected normalization and the TRRUST contribution explain this completely (**zero** residual).

## Distinct statistical families

| Quantity | Original saved 06b (case-sensitive, global BH) | Corrected local 06b and independent audit (uppercase, global BH) | Historical 04 (uppercase, includes TRRUST, per-layer BH) |
|---|---:|---:|---:|
| Reference terms | 4,420 | 4,421 | 4,422 |
| Tested feature–term hypotheses | 54,312,960 | 54,325,248 | 4 × 3,072 × 4,422 |
| Explicit tests with overlap ≥ 2 | 609,033 | 609,373 | 612,812 |
| Significant feature–term associations | 73,646 | 73,657 | 73,764 |
| Features with ≥1 reported significant annotation | 6,491 | 6,492 | 6,525 |

The historical 04 CSVs report a per-layer BH adjustment across 4,422 terms including TRRUST. The corrected 06b tests form one global family across the four decoder layers and 4,421 ontology terms. Case normalization added 340 explicit overlap tests in 06b. Adding the 3,439 TRRUST overlaps recovers the 612,812 tests in 04, with zero residual in each decoder layer.

In the local runs, **7/10** features in 06a and **19/98** in 06b passed both the global and annotation-degree-stratified fixed-panel maxT controls. These counts agree with the independent recalculation. The panel was partly selected on earlier BH results. In the independent recalculation, rebuilding the 98-feature panel after normalization replaced decoder.3 feature 1419 with feature 1421. The same 19 features passed both controls. The p-values apply to the selected panel.

## What is changed on this branch

- `notebooks/06a_statistical_significance.ipynb`: canonical uppercase gene universe and Top-20 matching, corrected ontology-family size, saved local results.
- `notebooks/06b_full_statistical_validation.ipynb`: canonical gene universe and reference matching, corrected ontology, hypothesis and association counts, saved local results.
- `README.md`, `reports/Research_Report.md`, `reports/Preliminary_Findings.md`, and `docs/index.html`: document the historical annotation family and the normalization correction. The preliminary report covers notebooks 01–05.
- `scripts/verify_gene_case.py`: lightweight read-only check of vocabulary normalization and local reference-term counts.

## Local verification

The original outputs are preserved at [`0812a7c`](https://github.com/vs-dx/pie-sae/commit/0812a7c0a7efb6926784138b5cf39580b3008dee). After the source correction, the affected outputs were cleared. The new run used the local project vocabulary, reference files, training catalogs and existing model artifacts. The project author saved the executed notebooks at [`eea7171`](https://github.com/vs-dx/pie-sae/commit/eea71712a62cf2e5b9061fcf132e5b730cc6e52c).

The vocabulary precheck can be repeated from the project root with `python .\scripts\verify_gene_case.py`. On the local run it reported 6,642 vocabulary genes, 70 symbols requiring normalization, and 4,421 eligible ontology terms (2,839 GO BP, 286 KEGG, 1,296 Reactome). The `C6ORF89` boundary check for `GO:0010648` passed.

The saved 06b output reports 73,657 rejections by both the BH step-up threshold and adjusted q-values. Its 19 internal consistency checks passed. The independent audit and the local rerun used the same underlying catalogs and references.

## Reading the executed notebooks

At execution commit [`eea7171`](https://github.com/vs-dx/pie-sae/commit/eea71712a62cf2e5b9061fcf132e5b730cc6e52c), the code cells contain the corrected run, while some Markdown cells and the notebook-level `pending_rerun` field still refer to the earlier revision. Cell numbers below are one-based and include Markdown cells.

| Saved output | Result to use |
|---|---|
| 06a §12.5, cell 55 | For feature 0, the global and degree-stratified nulls have 95th percentiles of 4.0671 and 3.8918, maxima of 5.0797 and 5.2967, and empirical p = 1/501 under both. |
| 06a §12.6 and §13, cells 58 and 62 | Seven of ten features passed both maxT controls. For decoder.0 feature 2120, global and degree-stratified maxT p-values are 0.03497 and 0.04296. |
| 06b §C.2–C.3, cells 53, 56 and 59 | All four historical layer counts reconcile. The corrected global BH run gives 73,657 associations across 6,492 features. |
| 06b §C.4, cell 63 | In the 98-feature panel, 22 passed global maxT, 19 passed degree-stratified maxT, and 19 passed both. The respective adjusted p-values are 0.00599 and 0.03497 for feature 0, and 0.35165 and 0.55445 for feature 2120. |
| 06b §E.3, cell 100 | All 19 numerical checks passed. The two fixed `NOT YET PERFORMED` messages describe the scope of the original 06b analysis. The subsequent independent audit verified annotations and enrichment. Notebook [07](../notebooks/07_causal_validation.ipynb) reports interventions in the fixed model. |

## Remaining independent concerns

**Validation dependence.** The pretrained PIE checkpoint was selected using the official validation set, and SAE configuration comparisons inspected validation performance. The sampled 128 train and 32 validation context–perturbation pairs do not overlap, but performance on the latter is not a fully untouched independent test estimate.

**Post-selection.** The ten-feature pilot and 98-feature panel do not establish selection-adjusted significance for all 12,288 decoder SAE features.

**Mechanistic inference.** Notebook 07 reports a strong downstream response for decoder.0 feature 189 → decoder.1 feature 413 in six chosen conditions. The target pair has cosine about 0.3023 versus control absolute cosine at most 0.0702. The controls are not geometry-matched; causal specificity and mediation are hypotheses, not findings. Donor patching improved the output-distance metric in only 3/9, 2/9, and 3/9 pairs across DE, LFC, and Δp.

**Independent computation.** Full independent reruns of PIE and SAE weights, HDF5 reconstruction, and all original 07 interventions remain outstanding.

## Provenance

- Original repository snapshot: `0812a7c`, 17 November 2026.
- Independent audit: `PIE_SAE_Audit_RU.pdf`, 21 November 2026; supplementary `pie_sae_audit.py`.
- Corrected local executions on branch `audit/case-normalization-2026-11-21`, recorded at [`eea7171`](https://github.com/vs-dx/pie-sae/commit/eea71712a62cf2e5b9061fcf132e5b730cc6e52c).
