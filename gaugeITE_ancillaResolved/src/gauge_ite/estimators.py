"""Class conditioning, prescribed reweighting, and sampling-cost diagnostics."""

from __future__ import annotations

import numpy as np
import pandas as pd


def effective_sample_size(
    total_shots: int,
    native_probabilities: np.ndarray,
    target_probabilities: np.ndarray,
) -> float:
    native = np.asarray(native_probabilities, dtype=float)
    target = np.asarray(target_probabilities, dtype=float)
    if native.shape != target.shape:
        raise ValueError("native and target probabilities must have the same shape.")
    if np.any((target > 0) & (native <= 0)):
        return 0.0
    denominator = float(np.sum(np.divide(target**2, native, where=native > 0)))
    return float(total_shots / denominator) if denominator > 0 else 0.0


def reweight_class_observables(
    class_observables: pd.DataFrame,
    native_classes: pd.DataFrame,
    target_prior: pd.DataFrame,
    observable_column: str = "energy_estimate",
    target_column: str = "target_probability",
    total_shots: int | None = None,
) -> dict:
    merged = (
        class_observables[["class_id", observable_column]]
        .merge(native_classes[["class_id", "probability"]], on="class_id")
        .merge(target_prior[["class_id", target_column]], on="class_id")
    )
    coverage = float(merged[target_column].sum())
    estimate = (
        float(np.sum(merged[target_column] * merged[observable_column]))
        if np.isclose(coverage, 1.0)
        else np.nan
    )
    result = {
        "estimate": estimate,
        "target_coverage": coverage,
        "weights": merged,
    }
    if total_shots is not None:
        result["effective_sample_size"] = effective_sample_size(
            total_shots,
            merged["probability"].to_numpy(),
            merged[target_column].to_numpy(),
        )
    return result

