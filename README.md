# PIE-SAE: Layer-wise Sparse Autoencoder Analysis of PIE

This project studies the internal representations of the pretrained **PIE (Perturbation Is Everything)** model with sparse autoencoders (SAEs). We train a separate TopK SAE for each of PIE's 14 encoder, processor, and decoder blocks, characterize the resulting feature dictionaries, and examine their biological annotations and correspondence across adjacent blocks.

The experiments build on the upstream [PIE](https://github.com/ArcInstitute/pie) and [bio-sae](https://github.com/Biodyn-AI/bio-sae) implementations. Their source revisions are pinned as Git submodules; all PIE-specific extraction, analysis, and reporting are provided in this repository.

## Research questions

- How accurately can sparse dictionaries reconstruct representations from different PIE blocks?
- How do SAE feature frequencies and dictionary geometry vary across model depth?
- Which output-gene associations and functional annotations occur in decoder SAE features?
- To what extent do matched features on adjacent blocks exhibit geometric, activation, gene-association, and functional correspondence?

## Experimental setup

| Component | Configuration |
|---|---|
| Pretrained model | PIE, K562 held-out context checkpoint |
| Representation blocks | 4 encoder, 6 processor, 4 decoder (14 total) |
| SAE | Upstream `TopKSAE` / `SAETrainer` from `bio-sae` |
| Input dimension / dictionary size | 768 / 3072 |
| Sparsity | TopK = 32 |
| Train / validation | 128 / 32 perturbation conditions |
| Optimization | Adam; 20 epochs; batch size 256; learning rate 3e-4 |
| Initialization | Layer-specific seeds 1000–1013; batch-order seed 42 |
| Encoder / processor positions | 512 latent slots per condition |
| Decoder positions | 6,642 output genes per condition |

The pretrained PIE model and the 14 trained SAEs remain fixed during feature analysis. The mean activation used for centering is estimated separately on each training split. Gene and latent positions within one perturbation condition are correlated observations, not independent experimental replicates.

## Notebook pipeline

| Notebook | Analysis |
|---|---|
| [`01_pie_architecture.ipynb`](notebooks/01_pie_architecture.ipynb) | Model architecture, checkpoint loading, forward hooks, activation shapes |
| [`02_extract_activations.ipynb`](notebooks/02_extract_activations.ipynb) | Train/validation activation extraction, HDF5 validation, gene map, provenance |
| [`03_train_layerwise_sae.ipynb`](notebooks/03_train_layerwise_sae.ipynb) | TopK SAE training and checkpoint manifests for all 14 blocks |
| [`04_sae_feature_analysis.ipynb`](notebooks/04_sae_feature_analysis.ipynb) | Feature statistics, SVD, decoder gene catalogs, ontology/network analyses |
| [`05_cross_layer_analysis.ipynb`](notebooks/05_cross_layer_analysis.ipynb) | Dictionary correspondence, co-activation, gene/annotation overlap, null comparisons |

Reusable analysis modules are under [`src/`](src/). A concise summary is available in [`reports/Preliminary_Findings.md`](reports/Preliminary_Findings.md) and machine-readable results in [`reports/tables/`](reports/tables/).

## Results

**SAE reconstruction.** Held-out reconstruction R² spans **0.8119–0.9821** across the 14 fitted SAE dictionaries. TopK constrains the number of active features per representation to approximately 32; it is a design choice rather than an independent quality metric.

**Feature statistics and dictionary geometry.** Decoder feature catalogs contain 3,072 features per layer. Across all 14 dictionaries, 200 of 43,008 directions have absolute cosine similarity above 0.5 with an individual top-50 SVD direction (20 above 0.7). The axis-wise test does not measure agreement with the entire leading singular subspace.

**Decoder functional annotations.** Respectively **1,694, 1,616, 1,603 and 1,612** features in decoder blocks 0–3 have at least one reported significant annotation from the GO Biological Process / KEGG / Reactome / curated TF analysis. Enrichments use the 6,642 PIE output genes as background, one-sided Fisher tests, and Benjamini–Hochberg adjustment; see notebook 04 for the exact testing protocol.

**Cross-layer correspondence.** Across 11 within-component adjacent transitions, conditional co-activation for geometrically matched SAE features is **0.339–0.518**, versus **0.009–0.015** under the feature-reassignment controls on a fixed validation subsample. The corresponding 200-permutation empirical p-values are 1/201 for each transition. This permutation count does not resolve family-wise significance at α=0.05 across all 11 comparisons (minimum possible Bonferroni-adjusted p≈0.0547).

Matched decoder features share on average **4.437, 4.577 and 4.894** genes among their validation Top-20 lists for the three decoder transitions, and **3.873, 4.103 and 4.435** significant functional terms, respectively. Gene-list overlap and annotation overlap are statistically related because the latter is derived from the former.

## Installation and local data

Clone the project with pinned upstream source trees:

```bash
git clone --recurse-submodules https://github.com/vs-dx/pie-sae.git
cd pie-sae
```

If submodules were not initialized:

```bash
git submodule update --init --recursive
```

The original experiments ran using **Python 3.12 on Windows**. Source-only checkout does not install Python dependencies or download the pretrained checkpoint and biological reference files. The complete set of installed dependency versions for this run has not been verified in the published repository; upstream installation requirements and notebook imports provide the current installation references.

The full local workflow expects the following artifacts:

| Local path | Source |
|---|---|
| `models/replogle_wdataset/k562/best_auprc.ckpt` | Pretrained PIE checkpoint, obtained separately |
| `data/activations/expanded_128_32/` | HDF5 activations generated in notebook 02 |
| `runs/sae_protocol_tests/` | SAE checkpoints generated in notebook 03 |
| `results/` | Analytic summaries, annotations, manifests and reference files generated during notebooks 02–05 |

Some references in notebook 04 are downloaded from external ontology/network sources. Checkpoints, large HDF5 files, and full reference databases are not included in this GitHub repository. The original run used Windows paths in portions of the checkpoint manifests; other operating systems may require path adaptation.

## Methods and interpretation

SAE fitting uses the pinned upstream TopK implementation. In the decoder matrix of shape `(768, 3072)`, the upstream training code normalizes **rows** (`dim=1`), not individual feature columns. The saved models retain that behavior.

For cross-layer comparisons, this project uses **signed cosine similarity** and **mutual nearest neighbors**, followed by held-out co-activation and decoder gene/annotation comparisons. The upstream `bio-sae` cross-layer study uses an absolute-cosine threshold and feature chains. These methods address related but distinct matching questions.

Functional enrichment and feature co-activation provide evidence of associations under the specified statistical controls. They do not by themselves establish perturbation-specific regulation or causal influence on PIE outputs. Correlation within experimental conditions, stable gene identity, feature prevalence, overlapping ontologies, and the finite permutation budget limit the scope of inference. SAE feature stability across repeated fits of the same block was not assessed.

## Source provenance

- [PIE — Arc Institute](https://github.com/ArcInstitute/pie), pinned revision `21a526add9ff3a81b4f962d9b16bbae570052c63`.
- [bio-sae — Biodyn-AI](https://github.com/Biodyn-AI/bio-sae), pinned revision `3ac41b1e7725c224dc1225a3311e8802eb96da85`.

During the original project audit, the pretrained PIE checkpoint and 14 saved SAE checkpoints matched the SHA-256 digests recorded in the local manifests, and both submodules matched the expected clean revisions. The large binary artifacts are retained outside Git.
