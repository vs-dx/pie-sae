# PIE-SAE: Layer-wise Sparse Representations, Statistical Controls, and Computational Interventions in PIE

**Research report — exploratory study**  
**Repository:** [vs-dx/pie-sae](https://github.com/vs-dx/pie-sae)  
**Report date:** 17 November 2026  
**Scope:** Notebooks 01–07 (including 06a and 06b)  
**Status:** Completed exploratory computational analyses; independent third-party replication, biological validation, and unique feature-level mediation remain outstanding.

## Abstract

Sparse autoencoders (SAEs) provide a way to approximate high-dimensional neural representations using a comparatively small set of active dictionary features. Whether these learned features identify biologically meaningful computations, rather than stable gene identities or geometric regularities, requires additional statistical and intervention-based analysis. We applied layer-wise TopK SAEs to the pretrained Perturbation Is Everything (PIE) model, fitting one 3,072-feature dictionary to each of its 14 encoder, processor, and decoder blocks. The models reconstructed selected held-out representations with R² values between 0.8119 and 0.9821. Decoder SAE features exhibited extensive associations with functional gene sets and measurable correspondence across adjacent blocks.

We subsequently examined gene-identity confounds, condition-level uncertainty, functional-enrichment null models, and cross-layer matching controls. The original saved 06b analysis of 12,288 decoder features against 4,420 case-sensitive ontology terms reported 73,646 Benjamini–Hochberg-significant feature–term associations across 6,491 features. An independent subsequent audit resolved the gene-identifier mismatch and recomputed 4,421 ontology terms, 73,657 associations, and 6,492 annotated features under the same global BH protocol. The corrected notebook has not yet been rerun locally. In a separately selected 98-feature panel, 19 features passed both fixed-panel maxT controls. These results are exploratory because of dependency, null-model, and post-selection limitations.

Targeted decoder.0 SAE-direction interventions changed PIE predictions. Feature 189 ablation affected all 32 selected validation conditions, but its advantage over matched controls was inconsistent. Feature-coefficient patching did not reliably move recipient predictions toward donor predictions. Ablating feature 189 consistently produced a large response in the geometrically matched downstream decoder.1 feature 413 in six selected conditions; nevertheless, the controls were not matched for the pair's geometric alignment, so unique causal mediation was not established. The study identifies reproducible computational phenomena and important negative controls, but it does not demonstrate biological causality or an independently validated mechanistic explanation.

## 1. Introduction

PIE is a pretrained model for biological perturbation prediction. Its internal representations may encode gene identity, perturbation context, and other structure. SAE-based analysis aims to express these representations through an overcomplete dictionary with sparse activations, thereby exposing candidate directions for interpretation and computational testing.

Three questions motivate this work. First, can SAEs reconstruct the internal states of PIE with a restricted number of active features? Second, are identified features stable across conditions and adjacent network blocks, and are their gene associations unusual under appropriate null models? Third, do interventions along these directions have distinctive consequences for model predictions and downstream representations?

The work is based on the pinned upstream [PIE](https://github.com/ArcInstitute/pie) and [bio-sae](https://github.com/Biodyn-AI/bio-sae) implementations. The SAE training procedure uses the upstream TopKSAE and SAETrainer; PIE-specific extraction, feature analysis, null comparisons, and interventions were developed in this repository. The biological interpretation of a gene-associated SAE feature is not assumed from its annotation alone.

## 2. Model, data, and experimental design

### 2.1 Pretrained PIE and internal states

The audited PIE checkpoint contains 79,106,918 parameters and uses 14 principal Transformer blocks: four encoder blocks, six processor blocks, and four decoder blocks. The hidden dimension is 768. Encoder and processor activations are associated with 512 latent positions per perturbation condition. Decoder activations are associated with 6,642 output-gene positions, processed during model inference in chunks.

The original PIE checkpoint was restored with strict parameter loading. Hook-based extraction and identity interventions were checked against unchanged forward predictions. The model weights were kept fixed throughout subsequent SAE analyses and interventions. [Source: notebook 01](../notebooks/01_pie_architecture.ipynb).

### 2.2 Activation extraction and condition sampling

The primary training/validation activation collection uses 128 training and 32 validation perturbation conditions. It was preceded by a 16/4-condition extraction pilot. The expanded HDF5 collection occupies approximately 14.52 GiB and covers all 14 blocks. Gene indices were checked against the PIE data pipeline and its ordered 6,642-gene vocabulary; the extracted HDF5 arrays were scanned for expected shapes, data types, complete condition coverage, and nonfinite values.

The 32 validation conditions were previously selected from a larger validation collection and include HepG2, Jurkat, and RPE1 contexts. They are not independent biological replicates; multiple latent or output-gene vectors from the same condition are correlated. The official PIE checkpoint was selected using validation data, and the same official validation pool was later inspected while comparing SAE training protocols. The sampled 128 train and 32 validation pairs have no overlapping context–perturbation pair, but are not an untouched model-selection-independent test set. Dataset and extraction provenance are recorded in local manifests that are not all stored in this public source repository. [Source: notebook 02](../notebooks/02_extract_activations.ipynb).

### 2.3 Layer-wise SAE training

One TopK SAE was trained for each block using a 768-dimensional input and a 3,072-dimensional sparse feature dictionary, with TopK = 32. The recorded training configuration was Adam, 20 epochs, batch size 256, and learning rate 3e-4. Layers used separate fixed initialization seeds 1000–1013 and a minibatch-order seed of 42.

Training used centering means estimated from the corresponding training split. All 14 SAE checkpoints were recorded and verified in the local training manifest; each layer was fitted once, so this series does not quantify feature stability across repeated same-layer initializations. The upstream SAE implementation normalizes rows of its saved decoder weight matrix rather than independently normalizing individual decoder feature columns. [Source: notebook 03](../notebooks/03_train_layerwise_sae.ipynb) and [training table](tables/sae_training_summary.csv).

## 3. Analysis methods

### 3.1 Reconstruction, usage and dictionary geometry

Held-out SAE reconstruction was assessed through R² and other saved metrics. Feature usage, activation frequency inequality, dictionary cosine similarities, and alignment of individual SAE directions with the top 50 singular-vector axes were evaluated across blocks. The singular-vector comparison concerns separate axes, not the principal subspace as a whole.

### 3.2 Decoder gene associations and functional annotations

Decoder SAE feature–gene catalogs rank output-gene positions using feature activation statistics. The leading 20 genes per feature were compared across training and validation conditions. Gene sets from GO Biological Process, KEGG, and Reactome were tested for over-representation using the PIE output-gene vocabulary as the background. The original annotation pipeline also included curated transcription-factor information and separate descriptive network analyses involving STRING and TRRUST.

Gene identity is an important confound: a feature that consistently responds at a set of output-gene positions need not track the variation induced by biological perturbations. Functional enrichment supplies an annotation under the specified statistical model, not evidence that the feature implements or regulates the associated process. [Source: notebook 04](../notebooks/04_sae_feature_analysis.ipynb).

### 3.3 Adjacent-block correspondence

SAE directions in adjacent blocks within the same architectural component were matched using signed cosine similarity and mutual nearest neighbors. Matched pairs were compared using aligned held-out feature co-activation, overlap of their decoder Top-20 gene lists, and overlap of enriched functional terms. Random partner reassignment and, in the statistical follow-up, signed-coordinate-permutation controls provided distinct baselines. Matching was not necessarily repeated inside every null procedure, limiting selection-aware inference. [Sources: notebooks 05](../notebooks/05_cross_layer_analysis.ipynb) [and 06b](../notebooks/06b_full_statistical_validation.ipynb).

### 3.4 Statistical robustness and null models

The 06a pilot studied shared condition variation, leave-one-out robustness, cell-context effects, and random projection controls. It originally evaluated ten selected decoder.0 features against 4,420 case-sensitive ontology terms (4,421 after case normalization) using 1,000 global gene-label and 1,000 annotation-degree-stratified permutations. A single-step maxT approach controlled the selected panel under each specified permutation scheme.

The 06b study reconstructed the annotation hypothesis family, examined training-versus-validation gene ranks, used context-stratified bootstrap procedures, tested the complete set of 12,288 decoder features against the original 4,420-term family, subsequently corrected to 4,421, and evaluated a separate 98-feature panel under fixed-panel maxT controls. It also revisited cross-layer null comparisons. Neither the 10-feature nor the 98-feature panel can be interpreted as an independent confirmatory sample after earlier feature selection. [Sources: notebooks 06a](../notebooks/06a_statistical_significance.ipynb) [and 06b](../notebooks/06b_full_statistical_validation.ipynb).

### 3.5 Computational interventions

Notebook 07 operated on the fixed pretrained PIE model and fixed SAE checkpoints. At decoder.0, a feature-specific intervention subtracted the decoded feature contribution from the original hidden representation. In schematic notation, for a representation h, centering mean μ, SAE feature coefficient z_j, and decoded feature direction w_j:

**h_intervened = h − z_j(h − μ) · w_j.**

This operator modifies the PIE hidden state; it does not enforce an exact zero value for feature j after re-encoding by the SAE. The procedure is therefore more precisely described as a decoded-component ablation than as an exact intervention on a causally independent latent variable.

Outcome comparisons used the model's DE logits, log-fold-change predictions, and delta-p logits, evaluated within each head rather than treating absolute differences across heads as directly comparable. Controls included identity hooks, fixed random directions, alternative learned SAE directions, alternative donors, reversed patching, and selected downstream feature responses.

Activation patching transferred the difference between donor and recipient feature coefficients along a selected SAE decoder direction. Downstream decoder.1 responses were measured at source-feature-active gene positions. [Source: notebook 07](../notebooks/07_causal_validation.ipynb).

## 4. Results

### 4.1 Reconstruction and dictionary properties

All 14 SAE fits completed. The observed held-out reconstruction R² was **0.8119–0.9821**, with the highest saved value at encoder.0 and the lowest at processor.5. The mean L0 was 32 in each dictionary; no inactive feature was reported among the evaluated validation positions for this fixed training series.

Across all 14 dictionaries, **200 of 43,008** feature directions had absolute cosine above 0.5 with an individual leading-50 singular-vector axis, while **20** exceeded 0.7. These values do not demonstrate singular-subspace equivalence or biological disentanglement.

The historical decoder annotation summaries reported the following feature counts:

| Decoder block | Total SAE features | Features with at least one historical annotation |
|---|---:|---:|
| decoder.0 | 3,072 | 1,694 |
| decoder.1 | 3,072 | 1,616 |
| decoder.2 | 3,072 | 1,603 |
| decoder.3 | 3,072 | 1,612 |

The 21 November 2026 independent audit reproduced these historical per-layer results and their term-overlap genes and adjusted p-values. The original full family includes 4,421 case-normalized ontology terms and one curated TRRUST TF term, with per-layer BH adjustment. Compare with Section 5. [Sources: training table](tables/sae_training_summary.csv), [annotation summary](tables/decoder_annotation_summary.csv), and [notebook 04](../notebooks/04_sae_feature_analysis.ipynb).

### 4.2 Cross-layer geometric and activation correspondence

Across **11 within-component adjacent transitions**, mutual-nearest-neighbor feature pairs were found using the fixed dictionaries. The median best-match signed cosine similarity ranged from **0.410–0.554** for encoder transitions and approximately **0.620–0.640** for processor transitions; decoder transitions had lower medians, approximately **0.241–0.269**.

On the saved held-out subsample, conditional co-activation for matched features was **0.339–0.518**, compared with **0.009–0.015** under feature-partner reassignment. All eleven 200-permutation tests reached the minimum unadjusted empirical p-value, **1/201 ≈ 0.00498**. With eleven comparisons, even the minimum Bonferroni-adjusted value is approximately **0.0547**; the sampled permutation budget does not establish family-wise rejection at 0.05.

For the three decoder transitions, mean validation Top-20 gene intersections were **4.437, 4.577, and 4.894** genes. Mean shared functional terms were **3.873, 4.103, and 4.435**. These observables are related and cannot be regarded as statistically independent confirmations. [Sources: notebook 05](../notebooks/05_cross_layer_analysis.ipynb) and [cross-layer table](tables/cross_layer_summary.csv).

### 4.3 Pilot-level statistical controls (06a)

For decoder.0 feature 0 on a ten-gene panel, the measured common-condition variation fraction was **55.476%**. A second within-context statistic, calculated with a different variance denominator, was **68.854%**; the two values are not a before/after improvement estimate. Context-stratified bootstrap intervals in notebook 06b were **41.133%–64.680%** and **54.998%–76.907%**, respectively, based on 2,000 bootstrap repetitions.

The original analysis found shared variation under condition permutations, but random linear projections of PIE also showed considerable condition-level variation. These findings do not identify a unique contribution from SAE feature 0.

Under the two fixed-panel, 1,000-permutation maxT nulls, **seven of ten** exploratory decoder.0 features passed the 0.05 threshold: **0, 2678, 2586, 189, 672, 3027, and 2120**. Feature 189 had the largest observed pilot enrichment statistic (**T = 9.464**). The inference is conditional on the previously explored ten-feature panel and the two null models. [Sources: notebooks 06a](../notebooks/06a_statistical_significance.ipynb) [and 06b](../notebooks/06b_full_statistical_validation.ipynb).

### 4.4 Full-dictionary statistical analysis (06b)

Cross-split catalogs covered **12,288 decoder features**, with **12,285** complete feature-level Top-20 catalogs in both train and validation splits. Feature-level stability was heterogeneous; features 0 and 189 showed comparatively high gene-ranking support, whereas feature 2120 was more sensitive to the association statistic and had sparse firing support among top-ranked genes.

The **original saved 06b run** used **12,288 × 4,420 = 54,312,960** feature–ontology hypotheses and reported **73,646** significant associations involving **6,491** features under a global BH procedure. An **independent audit dated 21 November 2026** found a gene-identifier case mismatch and recomputed **12,288 × 4,421 = 54,325,248** hypotheses, **609,373 explicit overlap-based tests**, **73,657 significant associations**, and **6,492 annotated features**. These audit-derived corrected values are **pending local rerun of the revised notebook**. The family-wide result is not a claim that all such features mediate biological processes: correlated ontology terms, selected gene rankings, and uncertain p-value dependence affect interpretation.

In the subsequent **98-feature exploratory panel**, **19** features passed both global and annotation-degree-aware fixed-panel maxT procedures at p < 0.05. The panel had been partly informed by the same enrichment findings; no selection-aware degree-controlled family-wise conclusion over the entire 12,288-feature population follows.

The statistical reanalysis also found that observed dictionary alignment can exceed selected orientation nulls without implying one-to-one feature identity. Mutual-nearest-neighbor fractions were lower than their geometric null means in ten of the eleven adjacent transitions. [Source: notebook 06b](../notebooks/06b_full_statistical_validation.ipynb).

### 4.5 Feature ablation and alternative-direction controls (07)

Decoder.0 feature **189** was ablated over all **32 selected validation conditions**, spanning three cell contexts. Measurable changes occurred in the three PIE output heads for all 32 conditions. Compared with five fixed norm-matched random directions using the same source activation gate, mean absolute changes were approximately **1.206× (DE logits), 1.300× (LFC), and 1.358× (delta-p logits)** the corresponding mean controls. Feature 189 exceeded the mean random-direction control in only **17/32, 19/32, and 18/32** conditions, respectively.

The relative advantages were less consistent under learned SAE-direction controls, and bootstrap intervals depend on the fixed control directions and selected validation conditions. Thus, these interventions demonstrate a computational response but do not establish a universally distinctive or biologically specific role for feature 189.

Comparison among decoder.0 features **0, 189, and 2120** revealed marked heterogeneity. Feature 2120 produced nonzero effects in **19/32** conditions, with the three largest conditions accounting for approximately **72%** of its total DE-logit effect. A mean over all conditions obscures this concentration. [Source: notebook 07, Sections B.5–B.9](../notebooks/07_causal_validation.ipynb).

### 4.6 Activation patching: a negative result

Feature-189 coefficient patching was examined in nine directed donor–recipient pairs. Under the localized comparison, patching reduced distance to donor predictions in only **3/9 pairs for DE logits**, **2/9 for LFC**, and **3/9 for delta-p logits**. Comparisons with alternative learned directions, decoy donors, and reversed coefficient displacement did not show reliable donor-state transfer.

The donor-specificity controls were not fully matched for displacement magnitude. Some local metrics also used different affected-gene masks across subexperiments; the separate C.4 union-mask controls are therefore not numerically interchangeable with the earlier C.3 measurements. These caveats do not reverse the qualitative conclusion: the selected patch did **not** provide consistent evidence of transferring donor predictions. [Source: notebook 07, Section C](../notebooks/07_causal_validation.ipynb).

### 4.7 Downstream response of the matched 189–413 pair

An earlier dictionary analysis had matched **decoder.0 feature 189** with **decoder.1 feature 413** as mutual nearest neighbors; the signed decoder-direction cosine similarity was approximately **0.3023**.

Across **six selected conditions** representing HepG2, Jurkat, and RPE1, upstream feature-189 ablation made downstream feature 413 the **highest-ranked absolute activation responder out of 3,072 decoder.1 SAE features** when measured at source-feature-active gene positions. The target upstream intervention also produced a larger feature-413 response than sixteen fixed, norm-matched learned SAE-direction controls in each selected condition.

However, the chosen source direction was more geometrically aligned with downstream feature 413 than the controls. The experimental design did not match the controls for that cosine alignment, nor did it implement a direct mediation test of feature 413. TopK competition may also affect apparent SAE activation changes. The result is a repeatable **computational response in the six selected conditions**, not proof that a unique feature-mediated causal pathway has been identified. [Source: notebook 07, Section E](../notebooks/07_causal_validation.ipynb).

## 5. Interpretation and unresolved limitations

### 5.1 Three levels of evidence

This report distinguishes:

1. **Association:** SAE activation patterns, enriched gene terms, and cross-layer correspondences under particular statistical or geometric baselines.
2. **Computational intervention:** Changes in the fixed PIE network after replacing or subtracting selected components of its hidden states.
3. **Biological causation:** Claims that the corresponding genes, pathways, or experimental cells are causally affected in the biological system.

The current study contains evidence at the first two levels under substantial qualifications. It does **not** provide independent evidence at the third level.

### 5.2 Historical annotation provenance

An independent audit (21 November 2026) traced the former discrepancy to inconsistent uppercase normalization of PIE gene symbols. The original 04 hypothesis family was reproduced with **4,421 GO BP / KEGG / Reactome terms plus one curated TRRUST TF term = 4,422 terms**, **73,764 significant feature–term associations**, and **6,525 annotated decoder features**, under *per-layer* BH correction. The previously described residual was **340 explicit tests**, not 337; its corrected residual is **zero**. This resolves historical annotation provenance, but does not remove the enrichment model's inferential limitations or substitute for a fresh run of the corrected 06b notebook. Details and the audit distinction are in [Case_Normalization_Correction.md](Case_Normalization_Correction.md).

### 5.3 Selection, dependence, and null-model assumptions

Features and conditions were examined during prior analyses, so later robustness, intervention comparisons, and empirical p-values must not be promoted to independent confirmatory tests. Output-gene positions are correlated within each perturbation condition; fixed-panel and condition-stratified uncertainty estimates do not generalize automatically to new contexts or unseen perturbations.

The permutation controls preserve different structures (feature overlap, annotation degree, fixed pairings, or aspects of geometric orientation). Rejection under one null does not imply rejection under all plausible nulls. Several tests involve highly dependent ontology terms and feature-associated gene lists; formal error-rate guarantees require assumptions not independently validated here.

### 5.4 Representation geometry and SAE identifiability

A learned SAE direction is not guaranteed to be a uniquely identifiable causal variable. Features can be correlated and sensitive to initialization, dictionary geometry, TopK competition, and off-manifold changes. The intervention operator subtracts a selected decoded component but does not enforce exact latent isolation; the target direction may overlap geometrically with other SAE features.

The 189–413 results therefore motivate **geometry-matched controls and downstream mediation experiments**, rather than establishing a unique causal chain.

### 5.5 Reproducibility boundaries

The original main-branch repository at commit `0812a7c` includes source code, all eight notebooks with saved results, three compact summary tables, pinned upstream submodule revisions, and documentation. The corrected 06a/06b versions in the audit branch have the dependent outputs cleared pending rerun. It excludes the approximately 14.5-GiB extracted activation data, pretrained PIE checkpoint, trained SAE weights, and some local analysis artifacts. A complete independent clean-environment rerun has not been performed. The saved 06b and 07 final cells report **19/19** and **15/15** internal consistency checks, respectively; these checks do not substitute for independent reproduction or peer review.

## 6. Conclusions

Fourteen TopK SAE dictionaries approximately reconstructed the sampled internal PIE representations at the chosen sparsity level. The decoder dictionaries exhibited extensive gene-associated enrichment and adjacent-block correspondence under the recorded statistical and permutation procedures. Expanded statistical auditing revealed both robust descriptive associations and limits due to feature selection, gene identity, dependence, and the need to rerun the corrected statistical notebooks.

Computational interventions demonstrated that selected SAE decoder directions affect PIE predictions and downstream representations. Nevertheless, feature-189 ablation was not uniformly superior to controls, and feature-coefficient patching failed to achieve consistent donor-prediction transfer. A strong downstream response was observed for the geometrically matched feature pair **189–413**, but the available controls do not distinguish special mediation from geometric alignment or TopK effects.

The principal contribution of this work is a connected, provenance-aware exploratory investigation—from activation extraction and SAE training to statistical checks and controlled interventions—that explicitly retains both positive and negative results. A mechanistic biological explanation remains an open hypothesis.

## 7. Recommended next experiments

1. **Rerun and verify the correction.** Re-execute the corrected 06a/06b notebooks using the original references and catalogs, compare the outputs with the independent audit, and preserve the historical results at commit `0812a7c`.
2. **Improve selection-aware inference.** Define held-out feature panels, fresh perturbation conditions, and appropriate family-level controls before testing prioritized findings.
3. **Test SAE specificity.** Compare with reconstruction-matched linear and nonlinear baselines, repeated same-layer SAE fits, and controls matched for activation prevalence, magnitude, and geometric alignment.
4. **Test mediation directly.** Compare upstream ablation, downstream feature-413 interventions, and joint or rescue interventions while accounting for off-manifold representations.
5. **Reproduce externally.** Export a dependency-locked environment and artifact manifest, then rerun the analysis in a clean computational environment.
6. **Seek biological validation cautiously.** Only advance selected candidate hypotheses to experimental validation after computational specificity and statistical robustness have been demonstrated.

## 8. Reproducibility and source index

All observations in this report are drawn from the published repository and its saved notebook outputs. The report is a synthesis of the available experimental record, **not an independent rerun**. Methods and output details should be checked against the corresponding notebooks.

| Source | Evidentiary role |
|---|---|
| [01 — Architecture](../notebooks/01_pie_architecture.ipynb) | PIE checkpoint, hooks, dimensions, identity checks |
| [02 — Extraction](../notebooks/02_extract_activations.ipynb) | Selected conditions, gene axis, HDF5 validation and provenance |
| [03 — Training](../notebooks/03_train_layerwise_sae.ipynb) | SAE protocol, checkpoints, reconstruction |
| [04 — Feature analysis](../notebooks/04_sae_feature_analysis.ipynb) | Dictionaries, SVD comparison, gene catalogs and annotations |
| [05 — Cross-layer](../notebooks/05_cross_layer_analysis.ipynb) | Matching and cross-layer correspondence |
| [06a — Statistical pilot](../notebooks/06a_statistical_significance.ipynb) | Selected-feature controls and maxT pilot |
| [06b — Statistical validation](../notebooks/06b_full_statistical_validation.ipynb) | Full-family statistical analyses, provenance audit and re-evaluation |
| [07 — Computational interventions](../notebooks/07_causal_validation.ipynb) | Ablation, patching, downstream response and final checks |
| [Training table](tables/sae_training_summary.csv) | 14 model reconstruction summaries |
| [Historical annotation table](tables/decoder_annotation_summary.csv) | Historical 04 decoder annotation counts; independently reproduced in the 21 November audit |
| [Cross-layer table](tables/cross_layer_summary.csv) | Eleven transition summaries |
| [Preliminary findings](Preliminary_Findings.md) | Earlier-stage report describing notebooks 01–05 |

**Upstream revisions:** PIE `21a526add9ff3a81b4f962d9b16bbae570052c63`; bio-sae `3ac41b1e7725c224dc1225a3311e8802eb96da85`.

**Citation status:** Exploratory repository report; not peer-reviewed. A formal manuscript would require additional reproducibility checks, explicit source bibliography, and independent methodological scrutiny.
