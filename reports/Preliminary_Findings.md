# Preliminary Findings: Layer-wise Sparse Autoencoder Analysis of PIE

## Experimental scope

Four encoder, six processor and four decoder blocks of a pretrained PIE model were analyzed with 14 independently initialized TopK sparse autoencoders. The training and validation splits comprise 128 and 32 perturbation conditions, respectively. Each SAE maps a 768-dimensional representation to 3,072 features with TopK=32.

We measured dictionary reconstruction, feature activation frequencies and geometry, 50-axis SVD alignment, output-gene associations, functional enrichment, and correspondence between adjacent blocks within each architectural component. PIE-specific extraction and analysis were implemented around the pinned upstream `bio-sae` TopK SAE and trainer.

## Quantitative findings

| Analysis | Observation |
|---|---|
| Held-out reconstruction | SAE R² = 0.8119–0.9821 across the 14 block-specific models |
| SVD alignment | 200/43,008 SAE directions have absolute cosine >0.5 with an individual leading-50 SVD direction; 20/43,008 exceed 0.7 |
| Decoder annotation coverage | 1,694 / 1,616 / 1,603 / 1,612 annotated features out of 3,072 in decoder blocks 0–3 |
| Matched co-activation | 0.339–0.518 versus 0.009–0.015 mean under feature-reassignment controls across 11 adjacent transitions |
| Decoder Top-20 gene intersection | 4.437 / 4.577 / 4.894 shared genes on held-out conditions across the three adjacent decoder transitions |
| Decoder functional terms | 3.873 / 4.103 / 4.435 shared enriched terms across the corresponding decoder transitions |

Matched features were identified using signed cosine similarity and mutual nearest-neighbor matching between independently trained SAE decoder dictionaries. The co-activation control uses 200 permutations for each transition; all 11 observed values lie beyond the sampled null distribution, yielding empirical unadjusted p=1/201. With 11 comparisons, the minimum achievable Bonferroni-adjusted p is about 0.0547 at this permutation resolution. Decoder gene and annotation correspondences were also evaluated against feature-reassignment controls.

The functional annotation workflow uses the PIE output-gene vocabulary as background and tests enrichment in GO Biological Process, KEGG and Reactome gene sets with one-sided Fisher tests and Benjamini–Hochberg correction. STRING and TRRUST interaction analyses are reported separately where applicable.

## Interpretation and limitations

Observed feature correspondences exceed the selected random matching baselines and exhibit partial stability in decoder gene associations. The co-activation, gene intersection and annotation results are complementary descriptions of related structures, not three statistically independent confirmations. In particular, shared output-gene identity may contribute to decoder correspondence.

The observations come from a fixed set of perturbation conditions, with many latent or gene-position vectors per condition. These vectors are not independent biological replicates. No repeated same-layer SAE refits or SAE-targeted causal interventions were performed. The matching rule also differs from the thresholded absolute-cosine feature-chain approach in the upstream `bio-sae` work.

The original SAE implementation normalizes rows of its decoder weight matrix, a behavior retained in the reported checkpoints. This detail affects the interpretation of dictionary weights and is documented in the repository README.

## Reproducible summaries

- [`sae_training_summary.csv`](tables/sae_training_summary.csv)
- [`decoder_annotation_summary.csv`](tables/decoder_annotation_summary.csv)
- [`cross_layer_summary.csv`](tables/cross_layer_summary.csv)

Full notebook outputs, the implementation of each statistical procedure, and additional qualitative analyses are in [`../notebooks/`](../notebooks/). Large data and trained checkpoints are stored separately from the source repository.
