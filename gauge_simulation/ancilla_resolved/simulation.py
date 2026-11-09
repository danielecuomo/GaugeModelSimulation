"""Small, direct calculations for first-time package users."""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from .models import TFIMSpec
from .thermodynamics import thermal_quantities_from_eigenvalues


def clean_spectrum(spec: TFIMSpec) -> np.ndarray:
    """Return the sorted spectrum of the clean TFIM Hamiltonian."""

    return np.linalg.eigvalsh(spec.clean_hamiltonian())


def clean_ground_energy(spec: TFIMSpec) -> float:
    """Return the exact ground-state energy of the clean model."""

    return float(clean_spectrum(spec)[0])


def clean_thermodynamics(
    spec: TFIMSpec,
    beta_values: Sequence[float],
) -> pd.DataFrame:
    """Return exact clean-model thermodynamics as a plain table."""

    betas = np.asarray(tuple(beta_values), dtype=float)
    if betas.ndim != 1 or len(betas) == 0:
        raise ValueError("beta_values must be a non-empty one-dimensional sequence.")
    if np.any(betas < 0):
        raise ValueError("beta_values must be non-negative.")
    spectrum = clean_spectrum(spec)
    return pd.DataFrame(
        [
            {
                "beta": float(beta),
                **thermal_quantities_from_eigenvalues(spectrum, float(beta)),
            }
            for beta in betas
        ]
    )


__all__ = ["clean_ground_energy", "clean_spectrum", "clean_thermodynamics"]
