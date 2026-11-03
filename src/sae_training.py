"""Reproducible PIE layer-wise SAE training using Kendiukhov's bio-sae classes.

Adaptations to PIE: read HDF5 activations; separate held-out conditions;
independent initialization and minibatch-order seeds; idempotent run output.
This module does not modify the upstream TopKSAE or SAETrainer.
"""

from __future__ import annotations

import gc
import json
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import torch


D_MODEL = 768
EXPANSION = 4
N_FEATURES = D_MODEL * EXPANSION
TOPK = 32
_LAYER_COUNTS = {"encoder": 4, "processor": 6, "decoder": 4}
_POSITIONS = {"encoder": 512, "processor": 512, "decoder": 6642}


def _load_layer(path: Path, component: str, index: int, expected_split: str,
                expected_conditions: int) -> np.ndarray:
    with h5py.File(path, "r") as f:
        if f.attrs["split"] != expected_split:
            raise ValueError(f"Unexpected split in {path}")
        if int(f.attrs["n_conditions"]) != expected_conditions:
            raise ValueError(f"Unexpected condition count in {path}")
        if int(f.attrs["completed_count"]) != expected_conditions:
            raise ValueError(f"Incomplete extraction: {path}")
        ds = f[component][str(index)]
        expected_shape = (expected_conditions, _POSITIONS[component], D_MODEL)
        if ds.shape != expected_shape or ds.dtype != np.float32:
            raise ValueError(f"Unexpected dataset: {component}.{index}: {ds.shape}, {ds.dtype}")
        x = ds[:].reshape(-1, D_MODEL)

    if not np.isfinite(x).all():
        raise ValueError(f"Non-finite activations: {component}.{index}, {path}")
    return np.ascontiguousarray(x)


@torch.inference_mode()
def _evaluate(sae, activations: np.ndarray, batch_size: int,
              device: torch.device) -> dict:
    sae.eval()
    n = len(activations)
    d = activations.shape[1]
    n_features = sae.n_features

    squared_error = 0.0
    input_sum = np.zeros(d, dtype=np.float64)
    input_squares = 0.0
    residual_sum = np.zeros(d, dtype=np.float64)
    residual_squares = 0.0
    feature_counts = np.zeros(n_features, dtype=np.int64)
    active_total = 0

    for start in range(0, n, batch_size):
        batch_np = activations[start:start + batch_size]
        x = torch.from_numpy(batch_np).to(device)
        reconstructed, codes, _ = sae(x)
        residual = (x - reconstructed).cpu().numpy().astype(np.float64)
        original = batch_np.astype(np.float64)
        active = (codes > 0)

        squared_error += float(np.square(residual).sum())
        input_sum += original.sum(axis=0)
        input_squares += float(np.square(original).sum())
        residual_sum += residual.sum(axis=0)
        residual_squares += float(np.square(residual).sum())
        feature_counts += active.sum(dim=0).cpu().numpy().astype(np.int64)
        active_total += int(active.sum().item())

    total_ss = input_squares - float(np.dot(input_sum, input_sum) / n)
    residual_ss = residual_squares - float(np.dot(residual_sum, residual_sum) / n)
    if total_ss <= 0:
        raise ValueError("Nonpositive validation variance")

    return {
        "mse": squared_error / (n * d),
        "r2": 1.0 - squared_error / total_ss,
        "explained_variance": 1.0 - residual_ss / total_ss,
        "l0": active_total / n,
        "inactive_features": int((feature_counts == 0).sum()),
        "n_vectors": n,
    }


def train_sae_layer(
    project_root: str | Path,
    component: str,
    layer_index: int,
    *,
    init_seed: int,
    order_seed: int = 42,
    epochs: int = 20,
    batch_size: int = 256,
    learning_rate: float = 3e-4,
    run_group: str = "sae_protocol_tests",
) -> dict:
    """Train one TopK SAE, or load an already completed matching run.

    A pre-existing incomplete or mismatched experiment raises an error.
    No previous checkpoints are overwritten.
    """
    root = Path(project_root).expanduser().resolve()
    if component not in _LAYER_COUNTS or not 0 <= layer_index < _LAYER_COUNTS[component]:
        raise ValueError(f"Invalid PIE layer: {component}.{layer_index}")
    if epochs <= 0 or batch_size <= 0 or learning_rate <= 0:
        raise ValueError("epochs, batch_size, learning_rate must be positive")

    upstream_src = root / "external" / "bio-sae" / "src"
    if not (upstream_src / "sae_model.py").is_file():
        raise FileNotFoundError(f"Official bio-sae implementation missing: {upstream_src}")
    if str(upstream_src) not in sys.path:
        sys.path.insert(0, str(upstream_src))
    from sae_model import SAETrainer, TopKSAE

    layer = f"{component}.{layer_index}"
    run_name = (f"{component}_{layer_index}_init{init_seed}_order{order_seed}"
                f"_ep{epochs}_bs{batch_size}_lr{learning_rate:.0e}")
    output_dir = root / "runs" / run_group / run_name
    config = {
        "layer": layer,
        "dataset": "expanded_128_32",
        "normalization": "original_dim1",
        "d_model": D_MODEL,
        "n_features": N_FEATURES,
        "topk": TOPK,
        "init_seed": int(init_seed),
        "order_seed": int(order_seed),
        "epochs": int(epochs),
        "batch_size": int(batch_size),
        "learning_rate": float(learning_rate),
    }

    # Idempotency: a completed, matching experiment is read, not trained again.
    if output_dir.exists():
        results_path = output_dir / "results.json"
        checkpoint_path = output_dir / "sae_final.pt"
        mean_path = output_dir / "activation_mean.npy"
        if not all(p.is_file() and p.stat().st_size > 0
                   for p in (results_path, checkpoint_path, mean_path)):
            raise RuntimeError(f"Incomplete existing run: {output_dir}")
        results = json.loads(results_path.read_text(encoding="utf-8"))
        if results.get("config") != config or not results.get("complete"):
            raise RuntimeError(f"Existing run configuration mismatch: {output_dir}")
        saved = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
        if (saved.get("config") != config or
                saved["model_state_dict"]["W_dec.weight"].shape !=
                (D_MODEL, N_FEATURES)):
            raise RuntimeError(f"Checkpoint validation failed: {output_dir}")
        saved_mean = np.load(mean_path)
        if saved_mean.shape != (D_MODEL,) or not np.isfinite(saved_mean).all():
            raise RuntimeError(f"Invalid saved mean: {mean_path}")
        print(f"Verified existing SAE, skipping training: {layer}")
        print(f"Validation R²: {results['final_validation']['r2']:.4f}")
        return results

    data_dir = root / "data" / "activations" / "expanded_128_32"
    x_train = _load_layer(data_dir / "train_128.h5", component, layer_index,
                          "train", 128)
    x_val = _load_layer(data_dir / "val_32.h5", component, layer_index,
                        "validation", 32)
    activation_mean = x_train.mean(axis=0, dtype=np.float64).astype(np.float32)
    x_train -= activation_mean
    x_val -= activation_mean

    # Independent seeds for weights and shuffling.
    torch.manual_seed(init_seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(init_seed)
    rng = np.random.default_rng(order_seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    sae = TopKSAE(d_model=D_MODEL, n_features=N_FEATURES, k=TOPK)
    trainer = SAETrainer(sae, lr=learning_rate, device=device)

    print(f"Training {layer} | device={device} | train={len(x_train)} | "
          f"val={len(x_val)} | init={init_seed} | order={order_seed}")
    history = []
    started = time.perf_counter()
    for epoch in range(1, epochs + 1):
        order = rng.permutation(len(x_train))
        losses = []
        for start in range(0, len(order), batch_size):
            batch = torch.from_numpy(x_train[order[start:start + batch_size]])
            losses.append(trainer.train_step(batch))

        row = {"epoch": epoch, "train_mse": float(np.mean(losses)),
               "optimizer_steps": int(trainer.step)}
        if epoch == 1 or epoch % 5 == 0 or epoch == epochs:
            valid = _evaluate(trainer.sae, x_val, batch_size, device)
            row.update({"val_mse": valid["mse"], "val_r2": valid["r2"]})
            print(f"Epoch {epoch:2d}/{epochs} | train MSE={row['train_mse']:.4f} | "
                  f"val MSE={valid['mse']:.4f} | R²={valid['r2']:.4f} | "
                  f"steps={trainer.step}")
        history.append(row)

    final = _evaluate(trainer.sae, x_val, batch_size, device)
    elapsed = time.perf_counter() - started
    upstream_commit = subprocess.check_output(
        ["git", "-C", str(root / "external" / "bio-sae"), "rev-parse", "HEAD"],
        text=True,
    ).strip()

    results = {
        "complete": True,
        "config": config,
        "train_vectors": len(x_train),
        "validation_vectors": len(x_val),
        "optimizer_steps": int(trainer.step),
        "bio_sae_commit": upstream_commit,
        "history": history,
        "final_validation": final,
        "elapsed_seconds": elapsed,
        "checkpoint": str((output_dir / "sae_final.pt").relative_to(root)),
    }

    # Save only after successful training, and write results last.
    output_dir.mkdir(parents=True, exist_ok=False)
    np.save(output_dir / "activation_mean.npy", activation_mean)
    temporary_ckpt = output_dir / "sae_final.tmp.pt"
    torch.save({"config": config,
                "model_state_dict": trainer.sae.state_dict(),
                "optimizer_state_dict": trainer.optimizer.state_dict()},
               temporary_ckpt)
    temporary_ckpt.replace(output_dir / "sae_final.pt")
    temporary_json = output_dir / "results.tmp.json"
    temporary_json.write_text(json.dumps(results, indent=2), encoding="utf-8")
    temporary_json.replace(output_dir / "results.json")
    print(f"Saved {layer}: R²={final['r2']:.4f} | L0={final['l0']:.2f} | "
          f"inactive={final['inactive_features']} | {output_dir}")

    del x_train, x_val, trainer, sae
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return results
