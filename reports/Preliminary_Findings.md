# Preliminary Findings — Layer-wise SAE Analysis of PIE

**Status:** exploratory analysis for initial scientific review; not a peer-reviewed manuscript or causal validation.

## Scope and methods

We inspected 14 blocks of the pretrained PIE K562 checkpoint (4 encoder, 6 processor, 4 decoder), extracted train and held-out activations using 128/32 perturbation conditions, and trained 14 separately initialized TopK sparse autoencoders using the pinned upstream `bio-sae` model and trainer. Each SAE maps 768-dimensional activations to 3072 features (TopK=32). We then analyzed feature frequencies, concentration, dictionary geometry, 50-axis SVD comparisons, decoder gene-association catalogs, pathway annotations and within-component adjacent-layer correspondence.

## Quantitative results

- **Reconstruction.** Held-out SAE R² ranges from 0.8119 to 0.9821 across the 14 block-specific models. Comparisons across blocks should account for different activation distributions and position semantics.
- **Dictionary geometry.** Median best-direction cosine for successive encoder pairs increases from 0.410 to 0.554, remains 0.620–0.640 for processor pairs and is 0.241–0.269 for decoder pairs. Signed cosine and mutual nearest neighbors were used without an absolute-cosine threshold. MNN counts alone are not evidence of stable feature identity.
- **Matched feature co-activation.** On the held-out fixed-position subsample, observed conditional co-activation is about 0.34–0.52, versus 0.009–0.015 under randomly reassigned matched pairs. Each of 11 pairs has the Monte-Carlo floor p=1/201 after 200 permutations; the results must not be reported as family-wise significant at 0.05 after Bonferroni correction (smallest adjusted p approximately 0.0547).
- **Decoder associations.** Adjacent decoder MNN feature pairs share on average 4.44, 4.58 and 4.89 validation Top-20 genes, respectively, and about 3.87, 4.10 and 4.43 enriched functional terms. Such agreement can reflect shared gene-aligned positions and correlated annotation databases; it does not prove conserved mechanism.
- **Annotation coverage.** The reported decoder annotation workflow associates 1603–1694 of 3072 features per decoder layer with at least one curated functional term. Check notebook 04 for the precise statistical background, term filters and null design.

Relevant machine-readable summaries: `reports/tables/sae_training_summary.csv`, `reports/tables/decoder_annotation_summary.csv`, and `reports/tables/cross_layer_summary.csv`.

## Major limitations and follow-up

1. **Restricted study size:** 128 training and 32 validation conditions; multiple positions within a condition are not independent biological samples.
2. **Training protocol:** row-wise decoder weight normalization inherited from upstream `bio-sae` (`dim=1` for a 768×3072 matrix), not per-feature unit norms. Existing checkpoints should be interpreted under the implementation actually used.
3. **Different matching criteria:** the feature matching scheme differs from the upstream 0.7 absolute-cosine, thresholded-chain methodology. It is an adaptation, not a direct full replication.
4. **Control design:** shuffling feature correspondence is informative but does not address all sources of shared activity; annotations can be correlated, and multiple comparisons matter.
5. **No causality:** neither pathway enrichment nor co-activation demonstrates regulatory circuitry. Causal patching and repeated-seed stability analysis remain future directions.
6. **Artifact availability:** large original data, activation HDF5 and trained checkpoints are omitted from the public repository; full third-party reproduction requires acquiring the upstream inputs and regenerating intermediates.

## Interpretation

Across the tested PIE checkpoint and sample, SAE dictionaries show nonrandom-looking cross-block co-activation relative to the chosen feature-reassignment control, while the decoder analysis exhibits shared gene-association and annotation structure despite relatively low geometric cosine. These are **associational, exploratory observations**. Strong claims about specific biological mechanisms, feature identity across layers, or causal computation require additional tests.
