# PIE-SAE: Case Normalization Correction and Reproducibility Record

**Status:** Corrected source code prepared on branch `audit/case-normalization-2026-11-21`; affected notebook outputs cleared **pending local re-execution**. Historical published outputs are preserved on `main` at commit `0812a7c`.

**Basis:** Independent critical audit dated 21 November 2026, supplied as `PIE_SAE_Audit_RU.pdf` and `pie_sae_audit.py`. The independent audit used saved Top-20 catalogs, external reference sets and selected official split files. It did **not** independently rerun PIE forward passes or train SAE checkpoints.

## Root cause

The original annotation function in `src/sae_biological_annotation.py` converts reference-set and PIE gene symbols to uppercase before comparing gene identities. In notebook 06a/06b the reference members were already uppercase, while `set(vocab["genes"])` retained mixed-case PIE gene names. Exact string matching then lost valid memberships. This was a gene-identity mismatch, not evidence that the historical 04 annotations used additional unknown biology.

The audit identified 70 mixed-case symbols among the 6,642 output genes, with no collisions after uppercase normalization. Twelve affected symbols restored relevant memberships; 49 ontology memberships and a previously filtered-out term were recovered. The GO BP boundary example is `GO:0010648`: `C6orf89` must match `C6ORF89` before the 5-gene minimum-size test.

The original `04` ontology family is **4,421 GO BP/KEGG/Reactome terms plus one curated TRRUST TF term = 4,422 terms**. The historical explicit overlap-test residual was incorrectly described as 337; its component totals sum to **340**. Corrected normalization and the TRRUST contribution explain this completely (**zero** residual).

## Distinct statistical families

| Quantity | Original saved 06b (case-sensitive, global BH) | Corrected external audit (uppercase, global BH) | Historical 04 (uppercase, includes TRRUST, per-layer BH) |
|---|---:|---:|---:|
| Reference terms | 4,420 | 4,421 | 4,422 |
| Tested feature–term hypotheses | 54,312,960 | 54,325,248 | 4 × 3,072 × 4,422 |
| Explicit tests with overlap ≥ 2 | 609,033 | 609,373 | 612,812 |
| Significant feature–term associations | 73,646 | 73,657 | 73,764 |
| Features with ≥1 reported significant annotation | 6,491 | 6,492 | 6,525 |

Do not overwrite the historical 04 CSV counts or silently mix the historical per-layer BH tests with the corrected global-BH 06b procedure. Per-layer and global BH adjust different families, so totals do not need to agree.

For the selected 06a 10-feature panel and 06b 98-feature panel, the independent audit reports that the same **7/10** and **19/98** features passed both tested fixed-panel maxT procedures after case normalization. The 98-feature panel was partly selected on earlier BH results, and its p-values are not selection-adjusted for the entire SAE dictionary.

## What is changed on this branch

- `notebooks/06a_statistical_significance.ipynb`: canonical uppercase gene universe and Top-20 matching, updated ontology-family dimension and explanatory notes.
- `notebooks/06b_full_statistical_validation.ipynb`: canonical gene universe and reference matching, recalculated expected ontology/hypothesis/association counts, updated explanatory notes.
- `README.md`, `reports/Research_Report.md`, `reports/Preliminary_Findings.md`, and `docs/index.html`: separate historical saved results from independent audit corrections and eliminate the claim of unexplained provenance.
- `scripts/verify_gene_case.py`: lightweight read-only check of vocabulary normalization and local reference-term counts.

**Output policy:** The affected code-cell outputs were cleared from the first modified code cell onward in both notebooks. The original outputs remain retrievable at `0812a7c`. Corrected values in this branch are attributed to the independent audit and are *not yet* newly executed notebook results. After an actual rerun, review every regenerated output and make sure the narrative matches the observed results.

## Local verification sequence

From the project root:

```powershell
python .\scripts\verify_gene_case.py
```

Then open the corrected 06a and 06b notebooks from this branch in their normal Python environment, execute the analysis sequentially using the original local gene vocabulary, reference files and training catalogs, and compare the results against the corrected column of the table above. Do not reset original local model artifacts, HDF5 files or checkpoints. Keep the old main-branch versions available for a direct comparison.

Only after the new notebook outputs and all diagnostic assertions have been verified should the branch be merged into `main` and the production atlas data regenerated/republished.

## Remaining independent concerns

**Validation dependence.** The pretrained PIE checkpoint was selected using the official validation set, and SAE configuration comparisons inspected validation performance. The sampled 128 train and 32 validation context–perturbation pairs do not overlap, but performance on the latter is not a fully untouched independent test estimate.

**Post-selection.** The ten-feature pilot and 98-feature panel do not establish selection-adjusted significance for all 12,288 decoder SAE features.

**Mechanistic inference.** Notebook 07 reports a strong downstream response for decoder.0 feature 189 → decoder.1 feature 413 in six chosen conditions. The target pair has cosine about 0.3023 versus control absolute cosine at most 0.0702. The controls are not geometry-matched; causal specificity and mediation are hypotheses, not findings. Donor patching improved the output-distance metric in only 3/9, 2/9, and 3/9 pairs across DE, LFC, and Δp.

**Independent computation.** Full independent reruns of PIE and SAE weights, HDF5 reconstruction, and all original 07 interventions remain outstanding.

## Provenance

- Original repository snapshot: `0812a7c`, 17 November 2026.
- Independent audit: `PIE_SAE_Audit_RU.pdf`, 21 November 2026; supplementary `pie_sae_audit.py`.
- Updated source: branch `audit/case-normalization-2026-11-21`.
- Do not interpret GitHub-hosted notebook code, markdown text, and audit-reported totals as a substitute for executing affected cells on the original local data.
