# PIE-SAE: Sparse Autoencoder Analysis of PIE Representations

This repository investigates the internal representations of the pretrained **PIE (Perturbation Is Everything)** model using layer-wise **sparse autoencoders (SAEs)**, statistical controls, and targeted computational interventions.

Fourteen TopK SAEs are fitted to PIE's encoder, processor, and decoder blocks. The study measures reconstruction, feature activation and geometry, gene-associated functional annotations, cross-layer correspondence, statistical robustness, and the effects of modifying selected SAE directions during PIE inference.

**Research status:** The planned exploratory analyses in notebooks 01–07 are complete. The findings describe properties of a fixed pretrained model and fixed SAE dictionaries. They **do not establish biological causality, unique feature-level mechanisms, or independent experimental replication**.

The project uses pinned upstream implementations of [PIE (Arc Institute)](https://github.com/ArcInstitute/pie) and [bio-sae (Biodyn-AI)](https://github.com/Biodyn-AI/bio-sae), included as Git submodules. PIE-specific extraction, analysis, statistical checks, and interventions are implemented in this repository.

## Research questions

1. How faithfully do layer-wise SAEs reconstruct PIE's internal representations?
2. How do feature usage, dictionary geometry, and gene-associated annotations vary across blocks?
3. Do geometrically matched features in adjacent blocks have related activations and gene associations?
4. Which patterns survive bootstrap, permutation, gene-label, and random-direction controls?
5. What changes in PIE's predictions and downstream representations when selected SAE components are ablated or patched?

## Experimental setup

| Component | Configuration |
|---|---|
| Pretrained model | PIE, K562 held-out-context checkpoint |
| Representation blocks | 4 encoder, 6 processor, 4 decoder (14 total) |
| SAE implementation | Pinned upstream `TopKSAE` / `SAETrainer` from `bio-sae` |
| Representation dimension | 768 |
| Dictionary size | 3,072 features per block |
| Sparsity | TopK = 32 |
| Selected activation dataset | 128 training and 32 validation perturbation conditions |
| SAE optimization | Adam; 20 epochs; batch size 256; learning rate 3e-4 |
| SAE initialization | Layer-specific seeds 1000–1013; batch-order seed 42 |
| Encoder / processor positions | 512 latent slots per condition |
| Decoder positions | 6,642 output-gene positions per condition |
| Feature analysis | Fixed pretrained PIE model and 14 fixed trained SAE checkpoints |

Means used to center SAE inputs are estimated from the relevant training activations. Observations at different gene or latent positions within a condition are **not independent biological replicates**. The 32-condition validation panel is a previously selected subset of a larger validation collection; follow-up analyses on that panel should be interpreted as exploratory rather than independent confirmatory tests. The pretrained PIE checkpoint was itself selected using the official validation split, and SAE protocol comparisons also inspected this validation data. Thus, reported hold-out SAE reconstruction metrics are not a pristine independent test assessment.

## Notebook pipeline

The recommended reading order is 01 → 02 → 03 → 04 → 05 → 06a → 06b → 07.

| Notebook | Main contents |
|---|---|
| [01 — PIE architecture](notebooks/01_pie_architecture.ipynb) | Architecture, checkpoint restoration, forward hooks, activation shapes, identity checks |
| [02 — Activation extraction](notebooks/02_extract_activations.ipynb) | Condition selection, HDF5 extraction, exact gene-axis checks, manifests and data validation |
| [03 — Layer-wise SAE training](notebooks/03_train_layerwise_sae.ipynb) | Independent initialization, 14 TopK SAE fits, held-out reconstruction and checkpoint verification |
| [04 — SAE feature analysis](notebooks/04_sae_feature_analysis.ipynb) | Feature usage, SVD comparisons, decoder gene catalogs, annotation and gene-label controls |
| [05 — Cross-layer correspondence](notebooks/05_cross_layer_analysis.ipynb) | Mutual-nearest-neighbor dictionary matching, co-activation, gene/term overlap and null comparisons |
| [06a — Statistical significance pilot](notebooks/06a_statistical_significance.ipynb) | Shared condition variation, random projections, selected-feature enrichment and two permutation nulls |
| [06b — Statistical validation](notebooks/06b_full_statistical_validation.ipynb) | Provenance audit, context-stratified bootstrap, cross-split stability, full-dictionary enrichment, selected-panel maxT and cross-layer controls |
| [07 — Computational interventions](notebooks/07_causal_validation.ipynb) | Identity checks, SAE-direction ablations, matched controls, activation patching, downstream responses and synthesis |

Reusable Python modules are in [`src/`](src/). Published tabular summaries from the earlier pipeline are in [`reports/tables/`](reports/tables/). [`Preliminary_Findings.md`](reports/Preliminary_Findings.md) is an **earlier-stage summary of notebooks 01–05**, not a synthesis of the later statistical and intervention experiments. For the later results, consult notebooks 06a, 06b and 07 directly.

## Results: reconstruction, annotation and correspondence (01–05)

**SAE reconstruction.** Held-out reconstruction R² ranges from **0.8119 to 0.9821** across the 14 SAE dictionaries ([training summary](reports/tables/sae_training_summary.csv)). TopK fixes approximately 32 active features per representation; this sparsity setting is not itself an independent quality metric.

**Dictionary geometry.** Across 14 dictionaries, **200 of 43,008** feature directions have absolute cosine similarity above 0.5 with an individual top-50 SVD direction; **20** exceed 0.7. This is an axis-wise comparison, not a test of alignment with the complete leading singular subspace.

**Decoder annotations.** The historical analysis reports **1,694 / 1,616 / 1,603 / 1,612** annotated features (out of 3,072) in decoder blocks 0–3. The original enrichment workflow uses the 6,642 PIE output genes as its background and one-sided over-representation tests with Benjamini–Hochberg adjustment. These historical per-layer counts were independently reproduced from the saved feature-gene lists and pinned annotation references in the 21 November 2026 audit. The historical procedure includes TRRUST and applies BH separately within each decoder layer.

**Adjacent-block correspondence.** For the 11 within-component adjacent transitions, geometrically matched SAE features show sampled conditional co-activation of **0.339–0.518**, compared with **0.009–0.015** under fixed-pair feature-reassignment controls ([cross-layer summary](reports/tables/cross_layer_summary.csv)). Each of the 11 comparisons had an unadjusted empirical p-value of 1/201 from 200 permutations. That permutation resolution does not permit a Bonferroni-adjusted rejection at 0.05 across all eleven tests (minimum adjusted value approximately 0.0547). The null does not repeat the original geometric matching procedure.

Across the three decoder transitions, matched features share on average **4.437 / 4.577 / 4.894** validation Top-20 genes and **3.873 / 4.103 / 4.435** significant functional terms. Gene-list and annotation overlap are dependent summaries of related information, not independent confirmations.

## Results: statistical analyses (06a–06b)

### Selected-feature pilot (06a)

For decoder.0 feature 0, the observed common-condition variation fraction was **55.476%** on a ten-gene panel; a within-context measure based on a different variance denominator was **68.854%**. Random linear projections of the underlying PIE representations also showed substantial shared variation, so these observations do not establish SAE-specific structure.

The original 06a run examined **4,420** case-sensitive GO Biological Process, KEGG, and Reactome terms. The case-corrected family contains **4,421** terms. In both versions, ten selected decoder.0 features were evaluated. Under global gene-label and annotation-degree-stratified permutation controls with 1,000 permutations each, **seven of the ten** passed both fixed-panel maxT procedures. The panel and its controls are exploratory, not a full-dictionary selection-adjusted test.

### Expanded statistical validation (06b)

The follow-up analysis compared all **12,288 decoder SAE features** across training and validation catalogs; **12,285** had complete Top-20 gene-association catalogs in both splits. Individual features differed markedly in stability and firing support, notably feature 2120.

The **original saved 06b run** evaluated **54,312,960 feature–term hypotheses** and reported **73,646 associations across 6,491 features** at q ≤ 0.05. An **independent 21 November 2026 audit** identified case-sensitive gene matching and recomputed the ontology-only global-BH family as **4,421 terms / 54,325,248 hypotheses / 73,657 associations across 6,492 features**. Corrected notebook code is available on the `audit/case-normalization-2026-11-21` branch, with new outputs **pending a local rerun**. These are results under the stated model; gene selection, dependent ontology terms, and assumptions required for formal FDR interpretation remain substantive limitations.

A subsequently examined **98-feature exploratory panel** yielded **19 features** passing both fixed-panel global and degree-aware maxT controls. The panel was partly selected using earlier results: these counts do **not** establish selection-adjusted significance over the complete 12,288-feature decoder population.

The expanded analysis also revisited feature-matching, co-activation, and gene-overlap results under distinct permutation controls. The controls test different null hypotheses and do not demonstrate causal coupling between features.

**Resolved historical annotation provenance (independent audit):** The original 04 analysis was reproduced as **4,421 case-normalized ontology terms plus one curated TRRUST TF term = 4,422 terms**, with **73,764 significant feature–term associations across 6,525 features**, using per-layer BH. The previous unexplained explicit-test residual was **340**, not 337, and falls to **zero** after correcting gene case and accounting for TRRUST. The corrected 06b global-BH results are a **different hypothesis family and correction protocol** from historical 04; see [the correction record](reports/Case_Normalization_Correction.md) and [06b](notebooks/06b_full_statistical_validation.ipynb).

## Results: computational interventions (07)

The intervention experiments kept PIE and the SAE dictionaries fixed, modifying hidden representations at decoder.0 during forward computation. The operator subtracts a decoded feature contribution from the original hidden state; it is **not an exact `do(z_j = 0)` intervention**, since a zero activation after SAE re-encoding was not enforced.

**Identity and baseline controls.** The restored PIE checkpoint reproduced the saved decoder.0 baseline in the initial condition, and identity interventions reproduced the original model outputs.

**Feature ablation.** Ablating decoder.0 feature **189** changed PIE outputs in all **32 selected validation conditions**. Mean absolute changes were about **1.21–1.36 times** the means for five fixed random, norm-matched directions across the three output heads. The feature exceeded the mean random-direction control in only **17–19 of 32 conditions**, depending on the output head; its advantage against other learned SAE directions was weaker. Effects are measured within each output head, whose numerical scales differ.

Features **0, 189 and 2120** were compared. Feature 2120 had effects concentrated in a subset of conditions: **19/32** had nonzero effects, and its top three conditions contributed approximately **72%** of its total DE-logit effect. This heterogeneity prevents a uniform per-condition interpretation.

**Activation patching.** Transferring feature-189 coefficients between selected donor and recipient conditions changed predictions but did not reliably move them toward donor predictions. Across nine directed pairs, donor distance improved in only **2–3 pairs** per output head. Alternative-direction, donor-specificity, and reverse-displacement controls did not support reliable donor-state transfer.

**Cross-layer response.** Feature **189** in decoder.0 was previously matched geometrically to feature **413** in decoder.1 (signed cosine approximately **0.3023**). In **six selected conditions**, ablating feature 189 made feature 413 the largest decoder.1 SAE activation responder among **3,072 features**, measured on source-feature-active gene positions. The response exceeded sixteen fixed learned-direction controls; however, those controls were **not matched for alignment with feature 413**. Geometry, TopK competition, and post-selection remain plausible explanations. **Unique mediation through 189 → 413 has not been demonstrated.**

The notebook reports **15/15 internal intervention consistency checks**; notebook 06b reports **19/19 internal statistical consistency checks**. These checks establish internal numerical consistency of the saved analyses, **not independent methodological or biological validation**.

## Interpretation and limitations

- **Association, computational intervention, and biological causation are different evidentiary levels.** The interventions demonstrate sensitivity of a pretrained model to modifications of selected hidden-state directions, not causal regulation in cells.
- Features, gene lists, and several validation conditions were selected using earlier analyses. Follow-up p-values, ranks and intervention comparisons are not independent confirmatory evidence.
- Gene positions within a perturbation condition are dependent; context-stratified bootstrap intervals rely on the observed contexts, selected conditions and stated resampling assumptions.
- Direction controls do not always match geometry, sparsity, coefficient magnitude, or the full distribution of interventions. SAE features are neither guaranteed independent causal variables nor necessarily unique across retrainings.
- Annotation terms overlap and hypothesis families are large. The historical term/test-count discrepancy has been **explained by case normalization**, but the revised 06b outputs still need to be regenerated locally and independently checked after rerun.
- SAE dictionaries were not independently refitted across random seeds for a stability analysis. Neither independent clean-environment replication nor wet-lab biological validation has been completed.

Negative findings — including inconsistent donor transfer during activation patching — are retained as results rather than omitted.

## Installation and local artifacts

Clone the repository with its pinned upstream source trees:

```bash
git clone --recurse-submodules https://github.com/vs-dx/pie-sae.git
cd pie-sae
```

If submodules are missing:

```bash
git submodule update --init --recursive
```

The original experiments ran with **Python 3.12 on Windows**. A clone of this source repository is **not an independently runnable reproduction package**: it does not bundle all Python dependency versions, pretrained weights, biological reference data, or the extracted activations. Follow the pinned upstream projects' environment instructions and the individual notebook imports, and provision the research artifacts before running the sequence.

| Local artifact | Role |
|---|---|
| `models/replogle_wdataset/k562/best_auprc.ckpt` | Original pretrained PIE checkpoint (acquired separately) |
| `data/activations/expanded_128_32/` | Train/validation HDF5 activations produced by notebook 02 |
| `runs/sae_protocol_tests/` | Fourteen trained SAE checkpoints and associated run outputs |
| `results/` | Local manifests, centering means, catalogs, annotation references, analysis outputs and figures |
| `external/pie/`, `external/bio-sae/` | Pinned upstream repositories |

Large model weights, HDF5 activation arrays (approximately 14.5 GiB for the saved expanded dataset), reference databases, and intermediate analysis artifacts are intentionally excluded from Git. Portions of the local manifests contain Windows-specific paths and require adaptation in other environments. Reproducing notebook 07 also requires the previously saved validation inputs and SAE checkpoints.

The original SAE implementation normalizes **rows** of its decoder weight matrix (`dim=1`) rather than independently normalizing the feature columns; the saved models retain this behavior. The weight layout and transposition used to obtain individual decoder directions are explicitly checked in notebook 07.

## Code and provenance

- [`src/`](src/) — reusable extraction, SAE training, feature statistics, catalog, reference-data and annotation helpers.
- [`reports/tables/`](reports/tables/) — published CSV summaries for SAE training, annotation coverage and cross-layer correspondence.
- [`reports/Preliminary_Findings.md`](reports/Preliminary_Findings.md) — historical preliminary report (01–05 only).
- [PIE upstream](https://github.com/ArcInstitute/pie), pinned submodule revision `21a526add9ff3a81b4f962d9b16bbae570052c63`.
- [bio-sae upstream](https://github.com/Biodyn-AI/bio-sae), pinned submodule revision `3ac41b1e7725c224dc1225a3311e8802eb96da85`.

The original local project checks reported matching SHA-256 digests for the pretrained PIE checkpoint and all fourteen trained SAE checkpoints. Those binary artifacts are not hosted in this repository and were not independently re-executed for this documentation update.

**Next research steps:** rerun the case-corrected notebooks and verify their outputs; independently reconstruct the environment and audit statistical assumptions; run selection-aware and geometry-matched controls; repeat SAE fits and validation in independent conditions; pursue biological experiments only where supported by a well-defined hypothesis.
