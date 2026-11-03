"""Layer-wise validation statistics for the independent PIE-SAE series.

Adapted from the activation-frequency computations in the project's archived
04_sae_feature_statistics notebook, using the pinned bio-sae TopKSAE class.

Inputs: independent_sae_manifest.json, val_32.h5, saved checkpoint and mean.
Outputs: one compressed NPZ per layer and a 14-layer JSON summary.
No model training and no PIE inference are performed.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import torch


LAYER_NAMES = (
    [f"encoder.{i}" for i in range(4)]
    + [f"processor.{i}" for i in range(6)]
    + [f"decoder.{i}" for i in range(4)]
)


def gini_coefficient(values: np.ndarray) -> float:
    """Gini coefficient of nonnegative per-feature firing frequencies."""
    x = np.sort(np.asarray(values, dtype=np.float64))
    if x.ndim != 1 or not np.isfinite(x).all() or (x < 0).any():
        raise ValueError("Expected finite nonnegative feature frequencies")
    total = x.sum()
    if total == 0:
        return 0.0
    n = x.size
    return float((2 * np.dot(np.arange(1, n + 1), x) / (n * total)) - (n + 1) / n)


def _sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def analyze_feature_statistics(
    project_root: str | Path,
    batch_size: int = 256,
) -> dict:
    """Compute and save held-out feature usage statistics for all 14 layers."""
    root = Path(project_root).resolve()
    manifest_file = root / "results/sae_training/independent_sae_manifest.json"
    validation_file = root / "data/activations/expanded_128_32/val_32.h5"
    upstream_file = root / "external/bio-sae/src/sae_model.py"
    output_dir = root / "results/sae_feature_analysis_independent/feature_statistics"
    for path in (manifest_file, validation_file, upstream_file):
        if not path.is_file():
            raise FileNotFoundError(path)
    if not isinstance(batch_size, int) or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")

    upstream_path = str(upstream_file.parent)
    if upstream_path not in sys.path:
        sys.path.insert(0, upstream_path)
    from sae_model import TopKSAE

    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    assert manifest["training_config"]["n_features"] == 3072
    assert manifest["training_config"]["topk"] == 32
    entries = manifest["checkpoints"]
    if len(entries) != 14:
        raise ValueError("Manifest must contain 14 checkpoints")
    registry = {entry["layer"]: entry for entry in entries}
    if set(registry) != set(LAYER_NAMES):
        raise ValueError("Manifest layer registry mismatch")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    results = []
    per_layer_output = {}

    with h5py.File(validation_file, "r") as h5:
        n_conditions = int(h5.attrs["n_conditions"])
        assert h5.attrs["split"] == "validation"
        assert n_conditions == 32 and h5.attrs["completed_count"] == 32

        for layer in LAYER_NAMES:
            component, layer_index = layer.split(".")
            record = registry[layer]
            checkpoint_path = root / record["checkpoint"]
            mean_path = checkpoint_path.parent / "activation_mean.npy"
            if not checkpoint_path.is_file() or not mean_path.is_file():
                raise FileNotFoundError(checkpoint_path)
            if _sha256(checkpoint_path) != record["checkpoint_sha256"]:
                raise ValueError(f"Checkpoint checksum mismatch for {layer}")
            checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
            config = checkpoint["config"]
            assert config["layer"] == layer
            assert config["init_seed"] == record["init_seed"]
            assert config["d_model"] == 768
            assert config["n_features"] == 3072
            assert config["topk"] == 32

            model = TopKSAE(d_model=768, n_features=3072, k=32)
            model.load_state_dict(checkpoint["model_state_dict"], strict=True)
            model = model.to(device).eval()
            mean_array = np.load(mean_path)
            assert mean_array.shape == (768,) and np.isfinite(mean_array).all()
            mean = torch.from_numpy(mean_array.astype(np.float32)).to(device)

            data = h5[component][layer_index]
            expected_positions = 6642 if component == "decoder" else 512
            assert data.shape == (32, expected_positions, 768)
            assert data.dtype == np.float32
            counts = np.zeros(3072, dtype=np.int64)
            positive_sums = np.zeros(3072, dtype=np.float64)
            n_vectors = 0
            positive_activations = 0

            with torch.inference_mode():
                for condition in range(32):
                    for start in range(0, expected_positions, batch_size):
                        x_np = np.asarray(
                            data[condition, start:min(start + batch_size, expected_positions), :],
                            dtype=np.float32,
                        )
                        if not np.isfinite(x_np).all():
                            raise ValueError(f"Nonfinite activation: {layer}, condition {condition}")
                        x = torch.from_numpy(x_np.copy()).to(device) - mean
                        z, _ = model.encode(x)
                        active = z > 0
                        counts += active.sum(dim=0).cpu().numpy().astype(np.int64)
                        positive_sums += z.sum(dim=0).cpu().numpy().astype(np.float64)
                        n_vectors += x.shape[0]
                        positive_activations += int(active.sum().item())

            frequencies = counts / n_vectors
            mean_active_amplitude = np.divide(
                positive_sums,
                counts,
                out=np.zeros(3072, dtype=np.float64),
                where=counts > 0,
            )
            alive = mean_active_amplitude[counts > 0]
            row = {
                "layer": layer,
                "init_seed": int(record["init_seed"]),
                "n_vectors": n_vectors,
                "n_features": 3072,
                "inactive_features": int((counts == 0).sum()),
                "mean_l0": float(positive_activations / n_vectors),
                "mean_frequency": float(frequencies.mean()),
                "median_frequency": float(np.median(frequencies)),
                "p95_frequency": float(np.percentile(frequencies, 95)),
                "max_frequency": float(frequencies.max()),
                "frequency_gini": gini_coefficient(frequencies),
                "median_feature_amplitude": float(np.median(alive)) if alive.size else None,
                "p95_feature_amplitude": float(np.percentile(alive, 95)) if alive.size else None,
            }
            assert abs(row["mean_l0"] - 32.0) < 1e-5, (layer, row["mean_l0"])
            assert row["inactive_features"] >= 0
            per_layer_output[layer] = {
                "feature_counts": counts,
                "activation_frequency": frequencies,
                "mean_active_amplitude": mean_active_amplitude,
                "total_vectors": np.array(n_vectors, dtype=np.int64),
            }
            results.append(row)
            print(
                f"{layer:<13} | N={n_vectors:>6} | "
                f"median freq={row['median_frequency']:.5f} | "
                f"Gini={row['frequency_gini']:.4f} | "
                f"median amplitude={row['median_feature_amplitude']:.3f} | "
                f"inactive={row['inactive_features']}"
            )
            del model

    # Commit output only after all 14 layers completed successfully.
    output_dir.mkdir(parents=True, exist_ok=True)
    for layer, arrays in per_layer_output.items():
        np.savez_compressed(output_dir / f"{layer.replace('.', '_')}.npz", **arrays)
    result = {
        "series": "independent",
        "source_manifest": str(manifest_file.relative_to(root)).replace("\\", "/"),
        "validation_hdf5": str(validation_file.relative_to(root)).replace("\\", "/"),
        "batch_size": batch_size,
        "frequency_definition": "fraction of validation vectors with SAE activation > 0",
        "amplitude_definition": "mean positive activation per feature, conditioned on firing",
        "gini_definition": "Gini coefficient of 3072 feature firing frequencies",
        "results": results,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"\nVerified all {len(results)} layers. Summary: {summary_path}")
    return result
