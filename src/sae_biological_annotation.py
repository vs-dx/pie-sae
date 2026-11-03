"""GO/KEGG/Reactome enrichment and descriptive STRING/TRRUST network summaries.

Adaptation of Biodyn-AI/bio-sae/src/04_annotate_features.py to PIE decoder
feature--gene catalogs. Uses the 6642 assayed output genes as the test universe.
The original script BH-corrects only p-values for terms with >=2 overlapping
members. Here omitted tests are retained in the multiple-testing denominator,
with p=1 (a conservative approximation), and the full prespecified
hypothesis count is used for BH correction.

Network links are described, NOT labelled statistically significant/enriched
without a separate network null model. TRRUST links remain directed.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy.stats import hypergeom

REFERENCE_FILES = {
    "GO_BP": "go_bp_gene_sets.json",
    "KEGG": "kegg_gene_sets.json",
    "Reactome": "reactome_gene_sets.json",
}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _bh_adjust_candidates(candidate_p: np.ndarray, total_hypotheses: int) -> np.ndarray:
    """BH adjusted p-values with omitted small-overlap tests set to p=1.

    This is conservative for the omitted tests, not mathematically identical
    to exact Fisher p-values when a term overlaps exactly one gene.
    The full hypothesis count is used, avoiding FDR inflation from filtering.
    """
    p = np.asarray(candidate_p, dtype=np.float64)
    if len(p) == 0:
        return p.copy()
    if total_hypotheses < len(p) or np.any((p < 0) | (p > 1)):
        raise ValueError("Invalid p-values or total hypothesis count")
    order = np.argsort(p, kind="stable")
    sorted_p = p[order]
    ranks = np.arange(1, len(p) + 1, dtype=np.float64)
    sorted_q = np.minimum.accumulate(
        (sorted_p * (total_hypotheses / ranks))[::-1]
    )[::-1]
    output = np.empty_like(p)
    output[order] = np.minimum(1.0, sorted_q)
    return output


def _prepare_terms(reference_dir: Path, vocabulary: set[str], min_size: int, max_size: int):
    # Each item is a dict with: source, term, genes.
    terms = []
    refs = {}
    for source, filename in REFERENCE_FILES.items():
        path = reference_dir / filename
        if not path.is_file():
            raise FileNotFoundError(path)
        raw = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(raw, dict):
            raise ValueError(f"Invalid gene set structure in {filename}")
        refs[source] = {"path": str(path), "sha256": _sha256(path)}
        for name, members in sorted(raw.items()):
            # Filter after intersection with the 6642 assayed genes.
            annotated = {str(g).upper() for g in members} & vocabulary
            if min_size <= len(annotated) <= max_size:
                terms.append({"source": source, "term": name, "genes": annotated})
    return terms, refs


def _load_networks(reference_dir: Path, vocabulary: set[str]):
    string_path = reference_dir / "string_ppi_edges.json"
    trrust_path = reference_dir / "trrust_human.tsv"
    for path in (string_path, trrust_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    raw = json.loads(string_path.read_text(encoding="utf-8"))
    pairs = {}
    for threshold in (700, 900):
        key = f"pairs_{threshold}"
        if key not in raw:
            raise ValueError(f"Missing STRING network tier {key}")
        links = set()
        for a, b in raw[key]:
            a, b = str(a).upper(), str(b).upper()
            if a != b and a in vocabulary and b in vocabulary:
                links.add(tuple(sorted((a, b))))
        pairs[threshold] = links
    if not pairs[900].issubset(pairs[700]):
        raise ValueError("STRING score-900 links must be in score-700 tier")

    directed = set()
    all_tfs = set()
    for line in trrust_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 4:
            raise ValueError(f"Invalid TRRUST record: {line[:100]}")
        tf, target = fields[0].upper(), fields[1].upper()
        if tf in vocabulary:
            all_tfs.add(tf)
        if tf in vocabulary and target in vocabulary and tf != target:
            directed.add((tf, target))
    refs = {
        "STRING": {"path": str(string_path), "sha256": _sha256(string_path)},
        "TRRUST": {"path": str(trrust_path), "sha256": _sha256(trrust_path)},
    }
    return pairs, directed, all_tfs, refs


def annotate_decoder_layer(
    project_root: str | Path,
    layer_index: int = 3,
    *,
    alpha: float = 0.05,
    min_term_size: int = 5,
    max_term_size: int = 500,
    min_overlap: int = 2,
) -> dict:
    """Annotate one PIE decoder layer; never runs PIE inference or SAE training.

    Tests GO/KEGG/Reactome and TF membership with one-sided hypergeometric
    (identical to one-sided Fisher exact with the fixed contingency margins).
    BH q-values computed across all terms x all features for this decoder layer;
    terms with fewer than two overlapping genes count as p=1 hypotheses.
    STRING and directed TRRUST interactions are descriptive, not p-tested.
    """
    if layer_index not in range(4):
        raise ValueError("Expected a decoder layer in [0, 3]")
    if not 0 < alpha < 1 or min_overlap < 1 or min_term_size < 1:
        raise ValueError("Invalid enrichment configuration")
    root = Path(project_root).resolve()
    ref_dir = root / "data" / "biological_references"
    catalog_path = (root / "results" / "sae_feature_analysis_independent"
                    / "feature_catalogs" / f"decoder_{layer_index}_train_catalog.json")
    vocab_path = root / "results" / "gene_reference" / "gene_index_to_symbol.json"
    for path in (catalog_path, vocab_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    vocab = json.loads(vocab_path.read_text(encoding="utf-8"))
    vocabulary = {str(g).upper() for g in vocab["genes"]}
    N = len(vocabulary)
    if N != 6642 or len(catalog["features"]) != 3072:
        raise ValueError("Incorrect PIE gene universe or feature count")
    if catalog["provenance"]["layer"] != f"decoder.{layer_index}":
        raise ValueError("Catalog layer mismatch")
    if catalog["provenance"]["split"] != "train":
        raise ValueError("Expected training catalog; validation is held out")

    terms, gene_set_refs = _prepare_terms(ref_dir, vocabulary, min_term_size, max_term_size)
    string_pairs, trrust_edges, trrust_tfs, network_refs = _load_networks(ref_dir, vocabulary)
    if len(trrust_tfs) >= min_term_size:
        terms.append({"source": "TRRUST_TF", "term": "curated_transcription_factors", "genes": trrust_tfs})
    n_terms = len(terms)
    if n_terms == 0:
        raise ValueError("No valid annotation terms after PIE vocabulary filtering")

    inverted = defaultdict(list)
    term_sizes = np.empty(n_terms, dtype=np.int32)
    for term_id, term in enumerate(terms):
        term_sizes[term_id] = len(term["genes"])
        for gene in term["genes"]:
            inverted[gene].append(term_id)

    features = catalog["features"]
    feature_gene_sets = []
    candidates_fi = []
    candidates_term = []
    candidates_overlap = []
    candidates_p = []
    network_stats = []

    print(f"Annotating decoder.{layer_index}: {len(features)} features, "
          f"{n_terms} eligible annotation terms, background={N} PIE genes")
    for i, feature in enumerate(features):
        if int(feature["feature_idx"]) != i:
            raise ValueError("Nonsequential SAE feature IDs")
        top = {str(item["gene_name"]).upper() for item in feature["top_genes"]}
        if not top.issubset(vocabulary) or len(top) > 20:
            raise ValueError(f"Invalid feature top genes at index {i}")
        feature_gene_sets.append(top)
        if len(top) >= min_overlap:
            counts = Counter(t for g in top for t in inverted[g])
            ids = np.fromiter((t for t, c in counts.items() if c >= min_overlap), dtype=np.int32)
            if ids.size:
                overlaps = np.array([counts[int(t)] for t in ids], dtype=np.int32)
                # Fisher's exact test (alternative='greater'), vectorized via hypergeom.sf.
                pvals = hypergeom.sf(overlaps - 1, N, term_sizes[ids], len(top))
                candidates_fi.extend([i] * len(ids))
                candidates_term.extend(ids.tolist())
                candidates_overlap.extend(overlaps.tolist())
                candidates_p.extend(np.asarray(pvals, dtype=np.float64).tolist())

        # STRING is an *undirected* functional association network.
        within_700 = 0
        within_900 = 0
        for a, b in combinations(sorted(top), 2):
            pair = (a, b)
            within_700 += pair in string_pairs[700]
            within_900 += pair in string_pairs[900]
        # TRRUST links preserve TF -> target direction.
        directed_links = sorted((tf, target) for tf, target in trrust_edges
                                if tf in top and target in top)
        network_stats.append({
            "feature_idx": i,
            "n_genes": len(top),
            "n_possible_undirected_pairs": len(top) * (len(top) - 1) // 2,
            "string_pairs_700": int(within_700),
            "string_pairs_900": int(within_900),
            "trrust_directed_pairs": len(directed_links),
            "trrust_directed_edges": [list(p) for p in directed_links],
        })
        if (i + 1) % 512 == 0:
            print(f"  Features: {i + 1}/{len(features)}", flush=True)

    # ALL terms x features are included in the BH family, even overlap=0/1
    # terms that are not explicitly evaluated (conservatively treated as p=1).
    m = n_terms * len(features)
    qvals = _bh_adjust_candidates(np.array(candidates_p), m)
    significant = []
    annotation_by_feature = defaultdict(list)
    significant_by_source = Counter()
    for fi, term_id, overlap_count, pval, q in zip(
        candidates_fi, candidates_term, candidates_overlap, candidates_p, qvals
    ):
        if q >= alpha:
            continue
        term = terms[term_id]
        hit_genes = sorted(feature_gene_sets[fi] & term["genes"])
        record = {
            "feature_idx": int(fi),
            "source": term["source"],
            "term": term["term"],
            "p_value": float(pval),
            "q_value": float(q),
            "overlap_count": int(overlap_count),
            "overlap_genes": hit_genes,
            "term_size_in_PIE": int(term_sizes[term_id]),
            "feature_gene_count": len(feature_gene_sets[fi]),
        }
        significant.append(record)
        annotation_by_feature[str(fi)].append(record)
        significant_by_source[term["source"]] += 1
    significant.sort(key=lambda r: (r["q_value"], r["p_value"], r["feature_idx"]))

    method = {
        "gene_background": "full PIE decoder output vocabulary",
        "gene_background_size": N,
        "term_filter": f"retain {min_term_size}..{max_term_size} member genes after intersecting with PIE vocabulary",
        "min_overlap_for_explicit_test": min_overlap,
        "omitted_gene_set_tests": "assigned p=1 (conservative; exact Fisher p for one-gene overlap is not evaluated)",
        "statistical_test": "one-sided Fisher exact; vectorized equivalent hypergeometric survival function",
        "multiple_testing": "Benjamini-Hochberg across ALL gene-set/TF-term tests x 3072 features, per decoder layer",
        "fdr_alpha": alpha,
        "network_metrics": "descriptive STRING 700/900 undirected links and directed TRRUST edges; no network p-values",
        "same_SAE_feature_train_validation": "train catalog used; validation withheld for independent assessment",
    }
    input_hashes = {
        "catalog_sha256": _sha256(catalog_path),
        "vocabulary_sha256": _sha256(vocab_path),
        **{f"{k}_sha256": v["sha256"] for k, v in {**gene_set_refs, **network_refs}.items()},
    }
    summary = {
        "layer": f"decoder.{layer_index}",
        "n_features": len(features),
        "n_ontology_terms": n_terms,
        "n_hypotheses_for_BH": m,
        "n_explicit_gene_set_tests_with_overlap_ge_2": len(candidates_p),
        "n_significant_feature_term_associations": len(significant),
        "n_features_with_significant_annotation": len(annotation_by_feature),
        "significant_by_source": dict(significant_by_source),
        "n_features_with_STRING_700_edges": sum(r["string_pairs_700"] > 0 for r in network_stats),
        "n_features_with_STRING_900_edges": sum(r["string_pairs_900"] > 0 for r in network_stats),
        "n_features_with_directed_TRRUST_edges": sum(r["trrust_directed_pairs"] > 0 for r in network_stats),
    }
    output = {
        "method": method,
        "inputs": input_hashes,
        "summary": summary,
        "significant_enrichments": significant,
        "feature_annotations": annotation_by_feature,
        "network_statistics": network_stats,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }
    out_dir = (root / "results" / "sae_feature_analysis_independent"
               / "biological_annotations")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"decoder_{layer_index}_annotations.json"
    if out_path.is_file():
        existing = json.loads(out_path.read_text(encoding="utf-8"))
        if existing.get("inputs") != input_hashes or existing.get("method") != method:
            raise FileExistsError(f"Existing annotation results have different provenance: {out_path}")
        print(f"Verified existing annotations; skipping recomputation: {out_path}")
        return {"path": str(out_path), "summary": existing["summary"], "reused": True}
    tmp = out_path.with_suffix(".partial.json")
    if tmp.exists():
        raise RuntimeError(f"Leftover incomplete annotation file: {tmp}")
    tmp.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp.replace(out_path)
    print(f"Saved annotations: {out_path}")
    return {"path": str(out_path), "summary": summary, "reused": False}
