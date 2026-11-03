"""PIE decoder SAE feature-to-gene catalogs.

Adaptation of Biodyn-AI/bio-sae/src/03_analyze_features.py:
* Rank genes by mean positive SAE feature activation (conditional on firing).
* Record feature-level activation frequency and mean activation across all positions.
* Record feature-gene firing counts, which equal supporting condition counts in PIE
  (one decoder position per gene per perturbation condition).
* Read checked HDF5 activations in batches, use original upstream TopKSAE,
  and resolve checkpoints from the independent-series manifest.

This module does NOT perform biological pathway enrichment or causal inference.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import torch

D_MODEL = 768
N_FEATURES = 3072
N_GENES = 6642
N_LAYERS = 4
TOPK = 32


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _accumulate_sparse(
    gene_ids: np.ndarray,
    topk_features: np.ndarray,
    topk_values: np.ndarray,
    gene_sum: np.ndarray,
    gene_count: np.ndarray,
    feature_sum: np.ndarray,
    feature_count: np.ndarray,
) -> int:
    """Accumulate strictly positive TopK activations; no full sparse codes on CPU."""
    if topk_features.shape != topk_values.shape:
        raise ValueError("TopK shape mismatch")
    if topk_features.ndim != 2 or len(gene_ids) != topk_features.shape[0]:
        raise ValueError("Position and feature indices do not match")

    active = topk_values > 0
    features = topk_features[active].astype(np.intp, copy=False)
    values = topk_values[active].astype(np.float64, copy=False)
    genes = np.broadcast_to(gene_ids[:, None], topk_features.shape)[active].astype(np.intp, copy=False)

    if not len(features):
        return 0
    if (features.min() < 0 or features.max() >= len(feature_count) or
            genes.min() < 0 or genes.max() >= gene_count.shape[1]):
        raise ValueError("Unexpected gene or feature index")

    feature_sum += np.bincount(features, weights=values, minlength=len(feature_count))
    feature_count += np.bincount(features, minlength=len(feature_count))
    # Each decoder gene occurs once per condition. Thus, for a gene-feature pair,
    # a firing count is also the number of experimental conditions supporting it.
    np.add.at(gene_sum, (features, genes), values)
    np.add.at(gene_count, (features, genes), 1)
    return len(features)


def _create_catalog(
    gene_sum: np.ndarray,
    gene_count: np.ndarray,
    feature_sum: np.ndarray,
    feature_count: np.ndarray,
    gene_symbols: list[str],
    n_positions: int,
    n_conditions: int,
    top_n: int,
) -> tuple[list[dict], dict]:
    """Match Kendiukhov's conditional-on-firing gene ranking."""
    feature_catalog = []
    for feature_idx in range(len(feature_count)):
        hits = int(feature_count[feature_idx])
        gene_idx = np.flatnonzero(gene_count[feature_idx])
        if len(gene_idx):
            conditional_mean = gene_sum[feature_idx, gene_idx] / gene_count[feature_idx, gene_idx]
            # Stable deterministic ordering: descending conditional mean, then gene index.
            order = np.lexsort((gene_idx, -conditional_mean))[:top_n]
            ranked = gene_idx[order]
            top_genes = [
                {
                    "gene_id": int(gid),
                    "gene_name": gene_symbols[int(gid)],
                    "mean_activation": float(gene_sum[feature_idx, gid] / gene_count[feature_idx, gid]),
                    "fire_count": int(gene_count[feature_idx, gid]),
                    "supporting_conditions": int(gene_count[feature_idx, gid]),
                }
                for gid in ranked
            ]
        else:
            top_genes = []

        feature_catalog.append({
            "feature_idx": feature_idx,
            "is_dead": hits == 0,  # only on selected data split
            "activation_freq": float(hits / n_positions),
            "mean_activation": float(feature_sum[feature_idx] / n_positions),
            "mean_activation_when_firing": float(feature_sum[feature_idx] / hits) if hits else 0.0,
            "fire_count": hits,
            "n_unique_genes": int(len(gene_idx)),
            "top_genes": top_genes,
        })
    summary = {
        "n_features": int(len(feature_count)),
        "n_alive": int(np.count_nonzero(feature_count)),
        "n_inactive_on_split": int(np.count_nonzero(feature_count == 0)),
        "n_conditions": n_conditions,
        "n_genes": len(gene_symbols),
        "n_positions": n_positions,
        "mean_l0": float(feature_count.sum() / n_positions),
        "top_genes_per_feature": top_n,
        "top_gene_ranking": "mean positive feature activation conditional on feature firing",
        "gene_feature_fire_count_meaning": "number of perturbation conditions where feature fired at that output gene",
    }
    return feature_catalog, summary


def build_decoder_catalog(
    project_root: str | Path,
    layer_index: int = 3,
    *,
    split: str = "train",
    batch_size: int = 256,
    top_genes: int = 20,
) -> dict:
    """Build or verify a catalog for a decoder layer on one HDF5 split.

    Saves JSON in results/sae_feature_analysis_independent/feature_catalogs.
    Existing matching outputs are reused; mismatches raise instead of overwrite.
    """
    if layer_index not in range(N_LAYERS):
        raise ValueError(f"Invalid decoder layer index: {layer_index}")
    if split not in {"train", "validation"}:
        raise ValueError("split must be 'train' or 'validation'")
    if batch_size <= 0 or not 1 <= top_genes <= N_GENES:
        raise ValueError("Invalid batch_size or top_genes")

    root = Path(project_root).resolve()
    manifest_path = root / "results" / "sae_training" / "independent_sae_manifest.json"
    gene_path = root / "results" / "gene_reference" / "gene_index_to_symbol.json"
    data_name, expected_conditions = ("train_128.h5", 128) if split == "train" else ("val_32.h5", 32)
    data_path = root / "data" / "activations" / "expanded_128_32" / data_name
    for path in (manifest_path, gene_path, data_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    layer_name = f"decoder.{layer_index}"
    matches = [row for row in manifest["checkpoints"] if row["layer"] == layer_name]
    if len(matches) != 1:
        raise ValueError(f"Expected one exact checkpoint for {layer_name}")
    record = matches[0]
    checkpoint_path = root / record["checkpoint"]
    current_sha256 = _sha256(checkpoint_path)
    if current_sha256 != record["checkpoint_sha256"]:
        raise ValueError(f"Checkpoint SHA-256 mismatch: {layer_name}")

    with gene_path.open(encoding="utf-8") as file:
        gene_symbols = json.load(file)["genes"]
    if len(gene_symbols) != N_GENES or len(set(gene_symbols)) != N_GENES:
        raise ValueError("Invalid ordered gene vocabulary")

    source_record = {
        "checkpoint": record["checkpoint"],
        "checkpoint_sha256": current_sha256,
        "bio_sae_commit": manifest["bio_sae_commit"],
        "gene_vocabulary_file": str(gene_path.relative_to(root)).replace("\\", "/"),
        "data_file": str(data_path.relative_to(root)).replace("\\", "/"),
        "data_file_size_bytes": data_path.stat().st_size,
        "data_file_mtime_ns": data_path.stat().st_mtime_ns,
        "layer": layer_name,
        "split": split,
        "batch_size": batch_size,
        "top_genes": top_genes,
    }

    output_dir = root / "results" / "sae_feature_analysis_independent" / "feature_catalogs"
    output_path = output_dir / f"decoder_{layer_index}_{split}_catalog.json"
    if output_path.exists():
        existing = json.loads(output_path.read_text(encoding="utf-8"))
        if existing.get("provenance") != source_record or len(existing.get("features", [])) != N_FEATURES:
            raise RuntimeError(f"Existing catalog does not match this configuration: {output_path}")
        print(f"Verified existing feature catalog; skipping recomputation: {output_path}")
        return {"path": str(output_path), "summary": existing["summary"], "reused": True}

    src = root / "external" / "bio-sae" / "src"
    if not (src / "sae_model.py").is_file():
        raise FileNotFoundError(src / "sae_model.py")
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    from sae_model import TopKSAE

    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    config = checkpoint["config"]
    if config["layer"] != layer_name or config["n_features"] != N_FEATURES or config["topk"] != TOPK:
        raise ValueError("Checkpoint config mismatch")
    model = TopKSAE(d_model=D_MODEL, n_features=N_FEATURES, k=TOPK)
    model.load_state_dict(checkpoint["model_state_dict"], strict=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device).eval()
    mean = np.load(checkpoint_path.parent / "activation_mean.npy")
    if mean.shape != (D_MODEL,) or not np.isfinite(mean).all():
        raise ValueError("Invalid training centering mean")
    mean_gpu = torch.from_numpy(mean).to(device=device, dtype=torch.float32)

    # Dense gene x feature accumulators, but no full activation dataset or SAE-code matrix.
    # About 245 MB for both (float64 sums, int32 counts) at 6642 x 3072.
    gene_sum = np.zeros((N_FEATURES, N_GENES), dtype=np.float64)
    gene_count = np.zeros((N_FEATURES, N_GENES), dtype=np.int32)
    feature_sum = np.zeros(N_FEATURES, dtype=np.float64)
    feature_count = np.zeros(N_FEATURES, dtype=np.int64)
    n_positions = 0
    started = time.perf_counter()

    print(f"Analyzing {layer_name}, split={split}, device={device}")
    with h5py.File(data_path, "r") as f:
        if f.attrs["split"] != split or int(f.attrs["completed_count"]) != expected_conditions:
            raise ValueError("Incorrect input HDF5 metadata")
        data = f["decoder"][str(layer_index)]
        if data.shape != (expected_conditions, N_GENES, D_MODEL) or data.dtype != np.float32:
            raise ValueError(f"Unexpected decoder activation shape: {data.shape}")
        gene_ids_ds = f["gene_ids"]
        if gene_ids_ds.shape != (expected_conditions, N_GENES):
            raise ValueError("Unexpected gene_id shape")

        with torch.inference_mode():
            for condition in range(expected_conditions):
                ids = np.asarray(gene_ids_ds[condition], dtype=np.int64)
                if not np.array_equal(ids, np.arange(N_GENES, dtype=np.int64)):
                    raise ValueError(f"Gene indices do not match gene vocabulary at condition {condition}")

                for start in range(0, N_GENES, batch_size):
                    end = min(start + batch_size, N_GENES)
                    raw = np.asarray(data[condition, start:end, :], dtype=np.float32)
                    if not np.isfinite(raw).all():
                        raise ValueError(f"Nonfinite activations at condition {condition}, rows {start}:{end}")
                    tensor = torch.from_numpy(raw.copy()).to(device)
                    sparse, indices = model.encode(tensor - mean_gpu)
                    top_values = sparse.gather(1, indices).cpu().numpy()
                    top_indices = indices.cpu().numpy()
                    _accumulate_sparse(
                        ids[start:end], top_indices, top_values,
                        gene_sum, gene_count, feature_sum, feature_count,
                    )
                    n_positions += end - start

                if (condition + 1) % 16 == 0 or condition + 1 == expected_conditions:
                    print(f"  Conditions: {condition + 1}/{expected_conditions}")

    expected_positions = expected_conditions * N_GENES
    if n_positions != expected_positions or int(feature_count.sum()) > TOPK * n_positions:
        raise RuntimeError("Incorrect total feature activation count")
    if not np.array_equal(gene_count.sum(axis=1), feature_count):
        raise RuntimeError("Gene-level and global feature counts disagree")
    if not np.allclose(gene_sum.sum(axis=1), feature_sum, rtol=1e-5, atol=1e-5):
        raise RuntimeError("Gene-level and global feature activations disagree")
    if np.any(gene_count > expected_conditions):
        raise RuntimeError("Gene-feature pair fires more often than available conditions")

    features, summary = _create_catalog(
        gene_sum, gene_count, feature_sum, feature_count,
        gene_symbols, n_positions, expected_conditions, top_genes,
    )
    summary["elapsed_seconds"] = time.perf_counter() - started
    output = {"provenance": source_record, "summary": summary, "features": features}

    output_dir.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(".partial.json")
    if temporary_path.exists():
        raise RuntimeError(f"Incomplete previous catalog found: {temporary_path}")
    temporary_path.write_text(json.dumps(output, indent=2, ensure_ascii=False), encoding="utf-8")
    temporary_path.replace(output_path)
    print(f"Catalog saved: {output_path}")
    print(f"Features: {summary['n_features']} | alive: {summary['n_alive']} | "
          f"mean L0: {summary['mean_l0']:.4f} | elapsed: {summary['elapsed_seconds']:.1f}s")
    return {"path": str(output_path), "summary": summary, "reused": False}
