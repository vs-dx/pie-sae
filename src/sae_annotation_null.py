"""Gene-label permutation controls for PIE decoder SAE functional annotations.

Preserve each feature's associated-gene set structure and gene popularity by
applying a single random bijection of PIE gene identities to all feature lists.
Recompute Fisher/hypergeometric tests and BH correction exactly as implemented
in sae_biological_annotation.py. First reproduce the published observed counts
before accepting the permutation control.

This tests whether gene IDs carry more functional coherence than chance under
this specific null, not causality or perturbation-dependent biological logic.
"""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy.sparse import csr_matrix
from scipy.stats import hypergeom

from sae_biological_annotation import (
    _bh_adjust_candidates,
    _load_networks,
    _prepare_terms,
)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as file:
        for chunk in iter(lambda: file.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def _reference_matrix(ref_dir: Path, names: list[str], min_size: int, max_size: int):
    names_upper = [name.upper() for name in names]
    if len(names_upper) != len(set(names_upper)):
        raise ValueError('Gene symbols in PIE vocabulary are not unique after uppercasing')
    gene_index = {name: i for i, name in enumerate(names_upper)}
    terms, _ = _prepare_terms(ref_dir, set(names_upper), min_size, max_size)
    _, _, tf_names, _ = _load_networks(ref_dir, set(names_upper))
    if len(tf_names) >= min_size:
        terms.append({'source': 'TRRUST_TF', 'term': 'curated_transcription_factors', 'genes': tf_names})
    if not terms:
        raise ValueError('No annotation terms')
    rows, cols = [], []
    sizes = np.empty(len(terms), dtype=np.int32)
    for j, term in enumerate(terms):
        sizes[j] = len(term['genes'])
        for name in term['genes']:
            rows.append(gene_index[name])
            cols.append(j)
    term_matrix = csr_matrix(
        (np.ones(len(rows), dtype=np.int16), (rows, cols)),
        shape=(len(names_upper), len(terms)), dtype=np.int16,
    )
    return gene_index, term_matrix, sizes, len(terms)


def _catalog_matrix(catalog: dict, gene_index: dict[str,int], n_genes: int):
    features = catalog['features']
    row_idx, col_idx = [], []
    for i, row in enumerate(features):
        if row['feature_idx'] != i:
            raise ValueError('Unexpected feature index')
        genes = {str(item['gene_name']).upper() for item in row['top_genes']}
        if len(genes) > 20 or not genes.issubset(gene_index):
            raise ValueError(f'Invalid feature-associated gene list: {i}')
        row_idx.extend([i] * len(genes))
        col_idx.extend(gene_index[g] for g in genes)
    return csr_matrix(
        (np.ones(len(col_idx), dtype=np.int16), (row_idx, col_idx)),
        shape=(len(features), n_genes), dtype=np.int16,
    )


def _enrichment_counts(
    features_by_genes: csr_matrix,
    genes_by_terms: csr_matrix,
    term_sizes: np.ndarray,
    n_genes: int,
    alpha: float,
):
    overlap = (features_by_genes @ genes_by_terms).tocoo()
    mask = overlap.data >= 2
    rows = overlap.row[mask]
    term_cols = overlap.col[mask]
    counts = overlap.data[mask]
    n_features, n_terms = features_by_genes.shape[0], genes_by_terms.shape[1]
    m = n_features * n_terms
    feature_sizes = np.diff(features_by_genes.indptr)
    if len(rows):
        pvals = hypergeom.sf(counts - 1, n_genes, term_sizes[term_cols], feature_sizes[rows])
        qvals = _bh_adjust_candidates(pvals, m)
        significant_rows = rows[qvals < alpha]
        n_assoc = int(len(significant_rows))
        n_annotated = int(np.unique(significant_rows).size)
    else:
        n_assoc = 0
        n_annotated = 0
    return {
        'annotated_features': n_annotated,
        'significant_associations': n_assoc,
        'tested_overlap_ge2': int(len(rows)),
    }


def run_annotation_permutation_control(
    project_root: str | Path,
    layer_index: int = 3,
    *,
    n_permutations: int = 25,
    seed: int = 42,
    alpha: float = 0.05,
    min_term_size: int = 5,
    max_term_size: int = 500,
):
    """Run or reuse validated null control for a decoder layer.

    Returns dict with observed and per-permutation null counts. Writes one JSON
    in `results/sae_feature_analysis_independent/annotation_controls`.
    """
    if layer_index not in range(4) or not 1 <= n_permutations <= 5000:
        raise ValueError('Invalid decoder layer or permutation count')
    root = Path(project_root).resolve()
    analysis_root = root / 'results' / 'sae_feature_analysis_independent'
    catalog_path = analysis_root / 'feature_catalogs' / f'decoder_{layer_index}_train_catalog.json'
    annotation_path = analysis_root / 'biological_annotations' / f'decoder_{layer_index}_annotations.json'
    vocab_path = root / 'results' / 'gene_reference' / 'gene_index_to_symbol.json'
    ref_dir = root / 'data' / 'biological_references'
    ref_names = ('go_bp_gene_sets.json', 'kegg_gene_sets.json',
                 'reactome_gene_sets.json', 'trrust_human.tsv', 'string_ppi_edges.json')
    paths = [catalog_path, annotation_path, vocab_path] + [ref_dir / x for x in ref_names]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)

    input_hashes = {str(path.relative_to(root)): _sha256(path) for path in paths}
    config = {
        'layer': f'decoder.{layer_index}', 'n_permutations': n_permutations,
        'seed': seed, 'alpha': alpha,
        'min_term_size': min_term_size, 'max_term_size': max_term_size,
        'null': 'one global random bijection of 6642 gene labels per replicate; feature co-occurrence and gene popularity preserved',
    }
    output_dir = analysis_root / 'annotation_controls'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f'decoder_{layer_index}_gene_permutation_{n_permutations}_seed{seed}.json'
    if output_path.exists():
        previous = json.loads(output_path.read_text(encoding='utf-8'))
        if previous.get('config') != config or previous.get('input_sha256') != input_hashes:
            raise ValueError(f'Existing null control has changed provenance: {output_path}')
        print(f'Reusing verified null control: {output_path}')
        return {'path': str(output_path), **previous}

    catalog = json.loads(catalog_path.read_text(encoding='utf-8'))
    annotation = json.loads(annotation_path.read_text(encoding='utf-8'))
    vocab = json.loads(vocab_path.read_text(encoding='utf-8'))
    gene_names = vocab['genes']
    n_genes = len(gene_names)
    if n_genes != 6642 or catalog['provenance']['layer'] != config['layer']:
        raise ValueError('Incorrect PIE catalog / gene vocabulary')
    if catalog['provenance']['split'] != 'train':
        raise ValueError('Expected training catalog')
    gene_index, ref_matrix, term_sizes, n_terms = _reference_matrix(
        ref_dir, gene_names, min_term_size, max_term_size
    )
    X = _catalog_matrix(catalog, gene_index, n_genes)
    if X.shape != (3072, 6642):
        raise ValueError('Unexpected catalog dimensions')
    if n_terms != annotation['summary']['n_ontology_terms']:
        raise ValueError('Annotation term count differs from original calculation')
    method = annotation['method']
    if method['fdr_alpha'] != alpha or method['min_overlap_for_explicit_test'] != 2:
        raise ValueError('Null parameters differ from original annotation')
    if method['gene_background_size'] != n_genes:
        raise ValueError('Unexpected original background')

    print(f'Calibrating decoder.{layer_index} against saved enrichment ...', flush=True)
    observed = _enrichment_counts(X, ref_matrix, term_sizes, n_genes, alpha)
    original_summary = annotation['summary']
    expected = (
        original_summary['n_features_with_significant_annotation'],
        original_summary['n_significant_feature_term_associations'],
        original_summary['n_explicit_gene_set_tests_with_overlap_ge_2'],
    )
    actual = (
        observed['annotated_features'], observed['significant_associations'],
        observed['tested_overlap_ge2'],
    )
    if actual != expected:
        raise AssertionError(f'Original annotation not reproduced! calculated={actual}, saved={expected}')
    print('Observed enrichment reproduced exactly:', actual, flush=True)

    rng = np.random.default_rng(seed)
    null_counts = []
    for r in range(n_permutations):
        perm = rng.permutation(n_genes)
        # Globally relabel each gene ID, preserving sizes and overlap patterns.
        permuted = csr_matrix(
            (X.data.copy(), perm[X.indices], X.indptr.copy()),
            shape=X.shape, dtype=np.int16,
        )
        stats = _enrichment_counts(permuted, ref_matrix, term_sizes, n_genes, alpha)
        null_counts.append(stats)
        if (r+1) % max(1, min(10, n_permutations // 5)) == 0 or r+1 == n_permutations:
            print(f'  Permutations {r+1}/{n_permutations} | annotated={stats["annotated_features"]}', flush=True)

    annotated = np.array([item['annotated_features'] for item in null_counts])
    associations = np.array([item['significant_associations'] for item in null_counts])
    report = {
        'config': config,
        'input_sha256': input_hashes,
        'n_reference_terms': n_terms,
        'n_features': int(X.shape[0]),
        'observed': observed,
        'null_replicates': null_counts,
        'summary': {
            'null_mean_annotated': float(np.mean(annotated)),
            'null_max_annotated': int(np.max(annotated)),
            'null_mean_associations': float(np.mean(associations)),
            'null_max_associations': int(np.max(associations)),
            'empirical_p_annotated': float((1 + np.sum(annotated >= observed['annotated_features'])) / (n_permutations+1)),
            'empirical_p_associations': float((1 + np.sum(associations >= observed['significant_associations'])) / (n_permutations+1)),
        },
        'limitations': [
            'Gene-label permutation breaks mapping to biology while preserving gene-popularity/co-occurrence structure; it is not a perturbation intervention.',
            'Gene set enrichment and network edge counts do not demonstrate causal regulatory computation.',
            'P-value resolution is limited by the number of permutations.',
        ],
        'generated_utc': datetime.now(timezone.utc).isoformat(),
    }
    partial = output_path.with_suffix('.partial.json')
    partial.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    partial.replace(output_path)
    print(f'Control saved: {output_path}')
    return {'path': str(output_path), **report}
