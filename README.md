# PIE-SAE: Sparse Autoencoder Analysis of PIE Representations

This repository investigates the internal representations of the pretrained **PIE (Perturbation Is Everything)** model using layer-wise **sparse autoencoders (SAEs)**, statistical controls, and targeted computational interventions.

Fourteen TopK SAEs are fitted to PIE's encoder, processor, and decoder blocks. The study measures reconstruction, feature activation and geometry, gene-associated functional annotations, cross-layer correspondence, statistical robustness, and the effects of modifying selected SAE directions during PIE inference.

**Research status:** The planned exploratory analyses in notebooks 01–07 are complete. The findings describe properties of a fixed pretrained model and fixed SAE dictionaries. They **do not establish biological causality, unique feature-level mechanisms, or independent experimental replication**.

The project uses pinned upstream implementations of [PIE (Arc Institute)](https://github.com/ArcInstitute/pie) and [bio-sae (Biodyn-AI)](https://github.com/Biodyn-AI/bio-sae), included as Git submodules. PIE-specific extraction, analysis, statistical checks, and interventions are implemented in this repository.

[Open the interactive Atlas](https://vs-dx.github.io/pie-sae/atlas/) · [Research overview](https://vs-dx.github.io/pie-sae/) · [Full research report](reports/Research_Report.md) · [Pipeline and release audit](reports/Pipeline_Coverage_and_Release_Audit.md)

**Start with decoder.0 feature 189.** In the Atlas, choose **Explore feature 189**. ASH1L leads its training and validation gene lists. The annotation table shows the original notebook-04 per-layer q-values. The **Intervention synthesis** link follows the same feature into notebook 07, where subtracting its decoded contribution changes PIE outputs in 32 selected conditions. Feature 413 in decoder.1 shows the strongest response in six selected conditions. Alignment with the selected direction limits the specificity of that result.

## Relationship to the Geneformer and scGPT study

This project adapts Ihor Kendiukhov’s layer-wise interpretability workflow to PIE. PIE uses 512 latent slots per condition in the encoder and processor and a gene-indexed decoder. Biological feature catalogs therefore cover the four decoder blocks. The published study also measures regulatory specificity in external CRISPRi data and evaluates a multi-tissue training control. Those questions remain open for PIE. The audit maps each method in §§5.1–5.13 to code and results, including the different ablation and co-activation protocols. [Reference paper, v1](https://arxiv.org/abs/2603.02952v1) · [Method-by-method comparison](reports/Pipeline_Coverage_and_Release_Audit.md#method-by-method-comparison)

## What the Atlas contains

| Entry point | Data and interpretation |
|---|---|
| Overview and layer explorer | Reconstruction and validation firing statistics for 14 dictionaries of 3,072 features each. |
| Feature, gene and pathway explorers | Training and validation Top-20 lists for decoder features, a training gene index, and 73,764 historical notebook-04 annotations across 6,525 features. |
| Corrected 06b summary | Global BH gives 73,657 associations across 6,492 features. Individual q-values in cards retain their notebook-04 provenance. |
| Cross-layer flow and interventions | 12,411 mutual cosine matches over 11 adjacent transitions and a summary of the saved model interventions, with links to notebook 07. |

The static research overview is in `docs/index.html`. The interactive application and compact JSON are in `docs/atlas/`. For local browsing from the repository root, run `python -m http.server 8000 --bind 127.0.0.1 --directory docs` and open `http://127.0.0.1:8000/atlas/`. A browser can explore the saved data without model weights or a running PIE process.

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

The corrected local run examined **4,421** GO Biological Process, KEGG and Reactome terms. The original case-sensitive run examined 4,420. In each run, **seven of ten** selected decoder.0 features passed both fixed-panel maxT controls at p < 0.05, using 1,000 global gene-label permutations and 1,000 annotation-degree-stratified permutations. The inference applies to this exploratory panel.

### Expanded statistical validation (06b)

The follow-up analysis compared all **12,288 decoder SAE features** across training and validation catalogs; **12,285** had complete Top-20 gene-association catalogs in both splits. Individual features differed markedly in stability and firing support, notably feature 2120.

Across **54,325,248 hypotheses** (12,288 features × 4,421 ontology terms), the corrected local 06b run found **73,657 feature–term associations across 6,492 features** at global Benjamini–Hochberg q ≤ 0.05. These counts agree with the independent recalculation of 21 November 2026. The original case-sensitive run tested 54,312,960 hypotheses and found 73,646 associations across 6,491 features. Some Markdown cells in the execution commit still describe the earlier run. The [correction record](reports/Case_Normalization_Correction.md#reading-the-executed-notebooks) identifies the affected sections. See the [saved 06a and 06b runs](https://github.com/vs-dx/pie-sae/commit/eea71712a62cf2e5b9061fcf132e5b730cc6e52c).

In a **98-feature exploratory panel**, 22 features passed global fixed-panel maxT, 19 passed degree-stratified maxT, and **19 passed both**. All 19 came from the group already significant under BH. The panel was partly selected on earlier enrichment results, which limits inference to the selected features.

The expanded analysis also revisited feature-matching, co-activation, and gene-overlap results under distinct permutation controls. The controls test different null hypotheses and do not demonstrate causal coupling between features.

**Historical annotations (04).** The independent audit reproduced 4,422 terms, 73,764 significant associations and 6,525 annotated features using per-layer BH. Compared with original 06b, normalizing gene identifiers restored 340 feature–term pairs with at least two shared genes. Adding the curated TRRUST transcription-factor set accounts for a further 3,439 pairs. The corrected 06b output reconciles the explicit-test counts in every decoder layer. The 06b global BH analysis uses only the 4,421 GO BP, KEGG and Reactome terms. See the [correction record](reports/Case_Normalization_Correction.md).

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
- Ontology terms share genes, and each feature's Top-20 list was selected from training data. The interpretation of global BH and the permutation tests depends on their respective null and dependence assumptions.
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

Arc Institute provides the prepared inputs for notebooks 01–02 in the [Replogle/Nadig dataset](https://huggingface.co/datasets/arcinstitute/PIE_replogle_nadig_essential/tree/20c9faef76fc96fdc809871340bd93499413dfe7/preprocessed) and the [K562 held-out-context splits](https://huggingface.co/datasets/arcinstitute/PIE_splits/tree/396ab9563175ee887750c9eed7ccaea6f5fdbf50/replogle_wdataset/unseen_ctx/k562).

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

The original local project checks reported matching SHA-256 digests for the pretrained PIE checkpoint and all fourteen trained SAE checkpoints. The local 06a and 06b reruns used existing project artifacts. Independent reconstruction of the complete pipeline still requires access to the model weights and activation data.

**Next research steps.** Independently reconstruct the environment and examine the statistical assumptions. Evaluate prespecified feature panels, geometry-matched controls, repeat SAE fits, and independent perturbation conditions before pursuing biological validation.
