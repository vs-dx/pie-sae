# PIE-SAE: Layer-wise Sparse Autoencoder Analysis of PIE

**Research status:** completed exploratory experiments; initial source-code and metadata audit performed. This repository documents analyses of the pretrained **PIE** model using independently initialized layer-wise sparse autoencoders (SAEs). The study is **not** a causal validation of biological mechanisms and is **not** a complete replication of all methods in the upstream `bio-sae` project.

## Research question

How do sparse features of the pretrained PIE (Perturbation Is Everything) model vary across its encoder, processor and decoder blocks? We study their activation statistics, dictionary geometry, gene associations and correspondence between adjacent blocks.

The project uses the upstream implementations as Git submodules:

- [PIE — Arc Institute](https://github.com/ArcInstitute/pie), pinned at commit `21a526add9ff3a81b4f962d9b16bbae570052c63`.
- [bio-sae — Biodyn-AI](https://github.com/Biodyn-AI/bio-sae), pinned at commit `3ac41b1e7725c224dc1225a3311e8802eb96da85`.

## Experimental design

| Aspect | Configuration |
|---|---|
| Pretrained model | PIE, K562 held-out context checkpoint |
| Representation blocks | 4 encoder, 6 processor, 4 decoder (14 total) |
| SAE type | Upstream `TopKSAE` and `SAETrainer` (`bio-sae`) |
| Input dimension / SAE features | 768 / 3072 |
| TopK | 32 |
| Training / validation conditions | 128 / 32 (fixed split) |
| SAE optimization | 20 epochs, batch size 256, Adam learning rate 3e-4 |
| Initialization seeds | 1000–1013 (one per block); minibatch-order seed 42 |

Each encoder/processor condition has 512 latent positions; each decoder condition has 6,642 gene-aligned positions. **Positions are not independent perturbation experiments.** Individual SAE dictionaries are trained separately for the 14 blocks, using training-only activation centering.

## Analysis notebooks

Run notebooks in numerical order from the project root or the `notebooks/` directory:

| Notebook | Purpose |
|---|---|
| [`01_pie_architecture.ipynb`](notebooks/01_pie_architecture.ipynb) | Inspect PIE model, 14 blocks, activation hooks and checkpoint. |
| [`02_extract_activations.ipynb`](notebooks/02_extract_activations.ipynb) | Build and validate train/validation HDF5 activation datasets; record model and dataset provenance. |
| [`03_train_layerwise_sae.ipynb`](notebooks/03_train_layerwise_sae.ipynb) | Train 14 independent SAEs and save model/metric manifests. |
| [`04_sae_feature_analysis.ipynb`](notebooks/04_sae_feature_analysis.ipynb) | Feature activation statistics, cosine geometry, SVD, decoder feature–gene catalogs and functional annotations. |
| [`05_cross_layer_analysis.ipynb`](notebooks/05_cross_layer_analysis.ipynb) | Mutual-nearest-neighbor dictionary matching, held-out feature co-activation, gene/annotation overlap and permutation controls. |

Supporting modules are in [`src/`](src/). Exploratory notebook archives, generated checkpoints, raw data and caches are intentionally excluded from the public repository.

## Selected results

- Across 14 SAEs, reported held-out reconstruction **R² ranges from 0.812 to 0.982**. The setting TopK=32 enforces roughly 32 active features per sample; this is not itself evidence of learned quality.
- For the 11 adjacent within-component block transitions, geometrically matched features have conditional co-activation probabilities approximately **0.34–0.52**, compared with mean shuffled-assignment controls of approximately **0.009–0.015** (fixed held-out subsample). The experiment used **200 permutations per transition**, so its minimum empirical unadjusted p-value is **1/201 ≈ 0.00498**. These p-values do **not** establish family-wise significance for 11 simultaneous tests (the minimum Bonferroni-adjusted p is ≈0.0547).
- Decoder mutual-nearest-neighbor feature pairs share an average **4.44–4.89 genes among their top 20 validation-associated genes** and approximately **3.87–4.43 shared enriched functional terms**, depending on the adjacent block pair. These correspondences are descriptive and not causal.
- Between **1,603 and 1,694** features per decoder layer have at least one enrichment annotation in the reported GO Biological Process / KEGG / Reactome / TRRUST-based workflow. Annotation dependencies mean that overlap among terms does not constitute independent mechanistic confirmation.

Small auditable summary tables are in [`reports/tables/`](reports/tables/); a concise interpretation and limitations are in [`reports/Preliminary_Findings.md`](reports/Preliminary_Findings.md).

## Getting started (source and environment)

```powershell
git clone --recurse-submodules <REPOSITORY-URL>
cd pie-sae
```

If already cloned without submodules:

```powershell
git submodule update --init --recursive
```

The development run used **Python 3.12** on Windows with a project-local virtual environment. The third-party dependencies and model-specific prerequisites must be installed from the pinned upstream projects; **a complete environment lockfile has not yet been published**. Do not assume cloning this repo alone downloads the PIE checkpoint or underlying biological data.

The notebooks expect (locally, **not** on GitHub):

- A pretrained PIE checkpoint at `models/replogle_wdataset/k562/best_auprc.ckpt` and the required PIE data assets.
- Extracted activation HDF5 files under `data/activations/expanded_128_32/` after notebook 02.
- Saved SAE checkpoints under `runs/sae_protocol_tests/` after notebook 03.
- Derived statistics, annotations and manifests under `results/` after notebooks 02–05.

The training manifest stores checkpoint paths in Windows form, so Linux/macOS reproduction may require path adaptation. Some biological reference annotations in notebook 04 require externally downloaded databases or network access. **The repository provides documented research code and small summary outputs, not a fully self-contained one-command reproduction.**

## Methodological distinctions and limitations

1. **SAE normalization:** the pinned upstream `bio-sae` implementation calls `F.normalize(W_dec.weight, dim=1)`, which normalizes rows of the `(768, 3072)` decoder matrix rather than each feature column. This behavior was retained for all existing checkpoints. Do not silently change it without retraining and reporting the protocol change.
2. **Cross-layer matching:** this project uses **signed cosine similarity** and **mutual nearest neighbors** between dictionary directions. That differs from the upstream `bio-sae` cross-layer procedure using absolute cosine with a 0.7 threshold and feature chains; they must not be described as identical methods.
3. **Statistical controls:** the permutation tests compare feature correspondences under a specified reassignment null. Correlated conditions, gene identity, feature prevalence and annotation dependencies can influence apparent agreement. The recorded p-value resolution and multiple comparisons restrict inferential claims.
4. **Biological interpretation:** feature–gene association and pathway enrichment do not demonstrate gene regulation or causal intervention mechanisms. No SAE-based causal patching was performed in this study.
5. **Reproducibility:** fixed manifests and SHA-256 checks were verified in the original environment, but external users need independently sourced checkpoints/data and compatible dependencies. Only summary tables are included in this release.

## Provenance

The project uses the pinned upstream PIE and bio-sae code through Git submodules, with PIE-specific extraction, independent layer-wise training, analyses and reporting in the local notebooks and `src/`. The original checkpoint and 14 SAE checkpoint hashes were verified against saved manifests during internal audit (15/15 files; both submodules at expected clean commits). Large model files and biological data are not distributed in this repository.

See the notebooks and [`reports/Preliminary_Findings.md`](reports/Preliminary_Findings.md) for implementation details and carefully qualified conclusions.
