"""Exact small-system thermodynamics for gauge-flux classes."""

from __future__ import annotations

from itertools import product
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.integrate import simpson, trapezoid
from scipy.special import logsumexp

from .labels import publication_label
from .models import TFIMSpec


def thermal_quantities_from_eigenvalues(
    eigenvalues: Sequence[float], beta: float
) -> dict:
    eigenvalues = np.asarray(eigenvalues, dtype=float)
    beta = float(beta)
    shifted = np.exp(-beta * (eigenvalues - eigenvalues[0]))
    normalized = shifted / shifted.sum()
    energy = float(np.sum(normalized * eigenvalues))
    energy_sq = float(np.sum(normalized * eigenvalues**2))
    log_partition = float(-beta * eigenvalues[0] + np.log(shifted.sum()))
    variance = max(energy_sq - energy**2, 0.0)
    return {
        "ground_energy": float(eigenvalues[0]),
        "thermal_energy": energy,
        "thermal_variance": variance,
        "heat_capacity": beta**2 * variance,
        "log_partition_function": log_partition,
        "partition_function": float(np.exp(log_partition)),
        "free_energy": -log_partition / beta if beta > 0 else np.nan,
    }


def branch_thermodynamics(
    spec: TFIMSpec,
    bond_signs: Sequence[int],
    field_signs: Sequence[int],
    beta: float,
) -> dict:
    hamiltonian = spec.branch_hamiltonian(bond_signs, field_signs)
    eigenvalues = np.linalg.eigvalsh(hamiltonian)
    return {
        **thermal_quantities_from_eigenvalues(eigenvalues, beta),
        "eigenvalues": eigenvalues,
    }


def class_representatives(
    spec: TFIMSpec,
    max_bonds: int = 16,
) -> pd.DataFrame:
    """Enumerate one effective bond representative for every cycle class."""

    if spec.num_bond_terms > max_bonds:
        raise ValueError(
            f"Class enumeration is limited to {max_bonds} bonds; "
            f"this model has {spec.num_bond_terms}."
        )
    classes: dict[tuple[int, ...], dict] = {}
    for bond_signs in product((1, -1), repeat=spec.num_bond_terms):
        fluxes = spec.lattice.cycle_fluxes(bond_signs)
        entry = classes.setdefault(
            fluxes,
            {
                "class_id": spec.lattice.class_id(bond_signs),
                "cycle_fluxes": fluxes,
                "plaquette_fluxes": spec.lattice.plaquette_fluxes(bond_signs),
                "flux_label": spec.lattice.format_fluxes(fluxes),
                "representative_bond_signs": tuple(bond_signs),
                "bond_multiplicity": 0,
            },
        )
        entry["bond_multiplicity"] += 1

    table = pd.DataFrame(classes.values()).sort_values("class_id").reset_index(drop=True)
    table["representative_field_signs"] = [
        (1,) * spec.num_system_qubits
    ] * len(table)
    table["full_sign_multiplicity"] = (
        table["bond_multiplicity"] * 2**spec.num_system_qubits
    )
    return table


def exact_class_thermodynamics(
    spec: TFIMSpec,
    beta_values: Sequence[float],
    max_system_qubits: int = 10,
) -> pd.DataFrame:
    """Compute exact Gibbs quantities for every gauge-flux class."""

    beta_values = np.asarray(beta_values, dtype=float)
    if beta_values.ndim != 1 or len(beta_values) == 0:
        raise ValueError("beta_values must be a non-empty one-dimensional sequence.")
    if np.any(beta_values < 0):
        raise ValueError("beta_values must be non-negative.")
    if spec.num_system_qubits > max_system_qubits:
        raise ValueError(
            f"Dense exact thermodynamics are limited to {max_system_qubits} "
            f"system qubits; received {spec.num_system_qubits}."
        )

    representatives = class_representatives(spec)
    rows = []
    for class_row in representatives.to_dict("records"):
        hamiltonian = spec.branch_hamiltonian(
            class_row["representative_bond_signs"],
            class_row["representative_field_signs"],
        )
        eigenvalues = np.linalg.eigvalsh(hamiltonian)
        for beta in beta_values:
            rows.append(
                {
                    **class_row,
                    "beta": float(beta),
                    **thermal_quantities_from_eigenvalues(eigenvalues, float(beta)),
                }
            )

    table = pd.DataFrame(rows)
    table["multiplicity_weight"] = (
        table["full_sign_multiplicity"]
        / table.groupby("beta")["full_sign_multiplicity"].transform("sum")
    )
    table["_log_annealed"] = (
        np.log(table["full_sign_multiplicity"].astype(float))
        + table["log_partition_function"]
    )
    table["_max_log"] = table.groupby("beta")["_log_annealed"].transform("max")
    table["_weight"] = np.exp(table["_log_annealed"] - table["_max_log"])
    table["annealed_class_weight"] = (
        table["_weight"] / table.groupby("beta")["_weight"].transform("sum")
    )
    return (
        table.drop(columns=["_log_annealed", "_max_log", "_weight"])
        .sort_values(["beta", "class_id"])
        .reset_index(drop=True)
    )


def within_class_invariance_report(
    spec: TFIMSpec,
    beta: float,
    samples_per_class: int = 4,
) -> pd.DataFrame:
    """Numerically verify that ideal Hamiltonians in one class are equivalent."""

    representatives = class_representatives(spec)
    candidates: dict[int, list[tuple[int, ...]]] = {
        int(row.class_id): [] for row in representatives.itertuples()
    }
    for bonds in product((1, -1), repeat=spec.num_bond_terms):
        class_id = spec.lattice.class_id(bonds)
        if len(candidates[class_id]) < samples_per_class:
            candidates[class_id].append(tuple(bonds))
        if all(len(values) >= samples_per_class for values in candidates.values()):
            break

    rows = []
    fields = (1,) * spec.num_system_qubits
    for class_row in representatives.itertuples():
        spectra = []
        energies = []
        logs = []
        for bonds in candidates[int(class_row.class_id)]:
            values = branch_thermodynamics(spec, bonds, fields, beta)
            spectra.append(values["eigenvalues"])
            energies.append(values["thermal_energy"])
            logs.append(values["log_partition_function"])
        reference = spectra[0]
        rows.append(
            {
                "class_id": int(class_row.class_id),
                "cycle_fluxes": class_row.cycle_fluxes,
                "flux_label": class_row.flux_label,
                "tested_branches": len(spectra),
                "max_spectral_deviation": max(
                    float(np.max(np.abs(values - reference))) for values in spectra
                ),
                "thermal_energy_spread": float(max(energies) - min(energies)),
                "log_partition_spread": float(max(logs) - min(logs)),
            }
        )
    return pd.DataFrame(rows)


def class_prior_from_random_bonds(
    spec: TFIMSpec,
    p_negative: float = 0.5,
) -> pd.DataFrame:
    """Class prior for independently drawn quenched bond signs."""

    if not 0 <= p_negative <= 1:
        raise ValueError("p_negative must lie in [0, 1].")
    weights: dict[int, float] = {}
    metadata: dict[int, dict] = {}
    for bonds in product((1, -1), repeat=spec.num_bond_terms):
        negatives = sum(sign == -1 for sign in bonds)
        probability = p_negative**negatives * (1 - p_negative) ** (
            spec.num_bond_terms - negatives
        )
        class_id = spec.lattice.class_id(bonds)
        fluxes = spec.lattice.cycle_fluxes(bonds)
        weights[class_id] = weights.get(class_id, 0.0) + probability
        metadata[class_id] = {
            "class_id": class_id,
            "cycle_fluxes": fluxes,
            "flux_label": spec.lattice.format_fluxes(fluxes),
        }
    rows = [
        {**metadata[class_id], "quenched_class_weight": probability}
        for class_id, probability in sorted(weights.items())
    ]
    table = pd.DataFrame(rows)
    table["quenched_class_weight"] /= table["quenched_class_weight"].sum()
    return table


def reconstruction_table(
    exact_classes: pd.DataFrame,
    native_class_probabilities: pd.DataFrame | None = None,
    quenched_prior: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Compare clean, annealed, quenched, and optional native energies."""

    rows = []
    for beta, group in exact_classes.groupby("beta", sort=True):
        clean = group.loc[group["class_id"] == 0].iloc[0]
        rows.append(
            {
                "beta": beta,
                "estimator_id": "ed_clean_sector",
                "estimator": publication_label("ed_clean_sector"),
                "energy": clean.thermal_energy,
            }
        )
        rows.append(
            {
                "beta": beta,
                "estimator_id": "ed_annealed",
                "estimator": publication_label("ed_annealed"),
                "energy": float(np.sum(group["annealed_class_weight"] * group["thermal_energy"])),
            }
        )
        if quenched_prior is not None:
            merged = group.merge(quenched_prior, on="class_id")
            rows.append(
                {
                    "beta": beta,
                    "estimator_id": "ed_uniform_quenched",
                    "estimator": publication_label("ed_uniform_quenched"),
                    "energy": float(
                        np.sum(merged["quenched_class_weight"] * merged["thermal_energy"])
                    ),
                }
            )
        if native_class_probabilities is not None:
            native = native_class_probabilities.loc[
                np.isclose(native_class_probabilities["beta"], beta)
            ]
            merged = group.merge(native[["class_id", "probability"]], on="class_id")
            rows.append(
                {
                    "beta": beta,
                    "estimator_id": "ed_circuit_weighted",
                    "estimator": publication_label("ed_circuit_weighted"),
                    "energy": float(np.sum(merged["probability"] * merged["thermal_energy"])),
                }
            )
    return pd.DataFrame(rows)


def _integral(values: np.ndarray, beta: np.ndarray, method: str) -> float:
    """Integrate one prefix, including the two-point Simpson fallback."""

    if len(beta) < 2:
        return 0.0
    if method == "trapezoid" or len(beta) == 2:
        return float(trapezoid(values, x=beta))
    if method == "simpson":
        return float(simpson(values, x=beta))
    raise ValueError("method must be 'simpson' or 'trapezoid'.")


def _integration_weights(beta: np.ndarray, method: str) -> np.ndarray:
    """Return the linear quadrature weights for every cumulative endpoint."""

    weights = np.zeros((len(beta), len(beta)), dtype=float)
    for endpoint in range(1, len(beta)):
        prefix = beta[: endpoint + 1]
        for index in range(endpoint + 1):
            basis = np.zeros(endpoint + 1, dtype=float)
            basis[index] = 1.0
            weights[endpoint, index] = _integral(basis, prefix, method)
    return weights


def reconstruct_class_log_partitions(
    class_energies: pd.DataFrame,
    *,
    dimension: int,
    energy_column: str,
    energy_stderr_column: str | None = None,
    method: str = "simpson",
) -> pd.DataFrame:
    """Reconstruct class log partition functions by thermodynamic integration.

    The input must contain one energy for every class on a common increasing
    beta grid beginning at zero.  Standard errors are propagated under the
    explicit assumption that estimates at different beta values are
    statistically independent.
    """

    if dimension <= 0:
        raise ValueError("dimension must be positive.")
    required = {"beta", "class_id", "flux_label", energy_column}
    missing = required.difference(class_energies.columns)
    if missing:
        raise ValueError(f"class_energies is missing columns: {sorted(missing)}")

    beta = np.sort(class_energies["beta"].astype(float).unique())
    if len(beta) == 0 or not np.isclose(beta[0], 0.0):
        raise ValueError("The thermodynamic-integration grid must begin at beta=0.")
    if np.any(np.diff(beta) <= 0):
        raise ValueError("The beta grid must be strictly increasing.")
    weights = _integration_weights(beta, method)
    log_dimension = float(np.log(dimension))
    rows = []
    for class_id, group in class_energies.groupby("class_id", sort=True):
        group = group.sort_values("beta").reset_index(drop=True)
        if not np.allclose(group["beta"].to_numpy(dtype=float), beta):
            raise ValueError("Every class must use the same beta grid.")
        energy = group[energy_column].to_numpy(dtype=float)
        integrals = weights @ energy
        if energy_stderr_column is None:
            integral_stderr = np.zeros_like(integrals)
        else:
            stderr = group[energy_stderr_column].to_numpy(dtype=float)
            if np.any(~np.isfinite(stderr)):
                raise ValueError(
                    f"{energy_stderr_column} must be finite for uncertainty propagation."
                )
            integral_stderr = np.sqrt((weights**2) @ (stderr**2))
        for index, source in group.iterrows():
            beta_value = float(source.beta)
            log_partition = log_dimension - float(integrals[index])
            rows.append(
                {
                    "beta": beta_value,
                    "class_id": int(class_id),
                    "cycle_fluxes": source.get("cycle_fluxes"),
                    "flux_label": source.flux_label,
                    "integrated_energy": float(integrals[index]),
                    "log_partition_reconstructed": log_partition,
                    "log_partition_stderr": float(integral_stderr[index]),
                    "free_energy_reconstructed": (
                        -log_partition / beta_value if beta_value > 0 else np.nan
                    ),
                    "free_energy_stderr": (
                        float(integral_stderr[index]) / beta_value
                        if beta_value > 0
                        else np.nan
                    ),
                    "integration_method": method,
                }
            )
    return pd.DataFrame(rows).sort_values(["beta", "class_id"]).reset_index(
        drop=True
    )


def quenched_free_energy_summary(
    exact_classes: pd.DataFrame,
    operational_classes: pd.DataFrame,
    quenched_prior: pd.DataFrame,
    *,
    dimension: int,
    method: str = "simpson",
) -> dict[str, pd.DataFrame]:
    """Reconstruct quenched thermodynamics and separate its error sources."""

    operational_logz = reconstruct_class_log_partitions(
        operational_classes,
        dimension=dimension,
        energy_column="energy_estimate",
        energy_stderr_column=(
            "energy_stderr"
            if "energy_stderr" in operational_classes
            and operational_classes["energy_stderr"].notna().all()
            else None
        ),
        method=method,
    )
    ed_grid_logz = reconstruct_class_log_partitions(
        exact_classes,
        dimension=dimension,
        energy_column="thermal_energy",
        method=method,
    )
    class_table = (
        exact_classes[
            [
                "beta",
                "class_id",
                "flux_label",
                "thermal_energy",
                "ground_energy",
                "log_partition_function",
            ]
        ]
        .merge(
            operational_classes[
                [
                    "beta",
                    "class_id",
                    "energy_estimate",
                    "energy_stderr",
                    "probability",
                ]
            ],
            on=["beta", "class_id"],
        )
        .merge(
            operational_logz[
                [
                    "beta",
                    "class_id",
                    "log_partition_reconstructed",
                    "log_partition_stderr",
                ]
            ].rename(
                columns={
                    "log_partition_reconstructed": "operational_ti_log_partition",
                    "log_partition_stderr": "operational_ti_log_partition_stderr",
                }
            ),
            on=["beta", "class_id"],
        )
        .merge(
            ed_grid_logz[
                ["beta", "class_id", "log_partition_reconstructed"]
            ].rename(
                columns={
                    "log_partition_reconstructed": "ed_grid_ti_log_partition"
                }
            ),
            on=["beta", "class_id"],
        )
        .merge(
            quenched_prior[["class_id", "quenched_class_weight"]],
            on="class_id",
        )
        .sort_values(["beta", "class_id"])
        .reset_index(drop=True)
    )
    class_table["operational_ti_log_partition_error"] = (
        class_table["operational_ti_log_partition"]
        - class_table["log_partition_function"]
    )
    class_table["ed_grid_ti_log_partition_error"] = (
        class_table["ed_grid_ti_log_partition"]
        - class_table["log_partition_function"]
    )

    rows = []
    log_dimension = float(np.log(dimension))
    for beta, group in class_table.groupby("beta", sort=True):
        beta = float(beta)
        prior = group["quenched_class_weight"].to_numpy(dtype=float)
        prior = prior / prior.sum()
        exact_logz = group["log_partition_function"].to_numpy(dtype=float)
        operational_logz_values = group[
            "operational_ti_log_partition"
        ].to_numpy(dtype=float)
        ed_grid_logz_values = group["ed_grid_ti_log_partition"].to_numpy(
            dtype=float
        )
        log_prior = np.log(prior)
        ed_quenched_energy = float(
            np.sum(prior * group["thermal_energy"].to_numpy(dtype=float))
        )
        operational_quenched_energy = float(
            np.sum(prior * group["energy_estimate"].to_numpy(dtype=float))
        )
        operational_energy_stderr = float(
            np.sqrt(
                np.sum(
                    (
                        prior
                        * group["energy_stderr"].to_numpy(dtype=float)
                    )
                    ** 2
                )
            )
        )
        ground_limit = float(
            np.sum(prior * group["ground_energy"].to_numpy(dtype=float))
        )
        if beta > 0:
            ed_fq = -float(np.sum(prior * exact_logz)) / beta
            operational_fq = -float(np.sum(prior * operational_logz_values)) / beta
            ed_grid_fq = -float(np.sum(prior * ed_grid_logz_values)) / beta
            ed_fann = -float(logsumexp(log_prior + exact_logz)) / beta
            operational_fann = -float(
                logsumexp(log_prior + operational_logz_values)
            ) / beta
            operational_fq_stderr = float(
                np.sqrt(
                    np.sum(
                        (
                            prior
                            * group[
                                "operational_ti_log_partition_stderr"
                            ].to_numpy(dtype=float)
                        )
                        ** 2
                    )
                )
                / beta
            )
        else:
            ed_fq = operational_fq = ed_grid_fq = np.nan
            ed_fann = operational_fann = np.nan
            operational_fq_stderr = np.nan
        rows.append(
            {
                "beta": beta,
                "ed_uniform_quenched_energy": ed_quenched_energy,
                "operational_uniform_quenched_energy": operational_quenched_energy,
                "operational_uniform_quenched_energy_stderr": (
                    operational_energy_stderr
                ),
                "ed_quenched_ground_energy": ground_limit,
                "ed_uniform_quenched_free_energy": ed_fq,
                "operational_ti_uniform_quenched_free_energy": operational_fq,
                "operational_ti_uniform_quenched_free_energy_stderr": (
                    operational_fq_stderr
                ),
                "ed_grid_ti_uniform_quenched_free_energy": ed_grid_fq,
                "ed_annealed_free_energy": ed_fann,
                "operational_ti_annealed_free_energy": operational_fann,
                "ed_quenched_minus_annealed": ed_fq - ed_fann,
                "operational_ti_quenched_error": operational_fq - ed_fq,
                "ed_grid_quadrature_error": ed_grid_fq - ed_fq,
                "ed_reduced_quenched_free_energy": (
                    beta * ed_fq + log_dimension if beta > 0 else 0.0
                ),
                "operational_ti_reduced_quenched_free_energy": (
                    beta * operational_fq + log_dimension if beta > 0 else 0.0
                ),
                "ed_reduced_annealed_free_energy": (
                    beta * ed_fann + log_dimension if beta > 0 else 0.0
                ),
            }
        )
    return {
        "operational_log_partitions": operational_logz,
        "ed_grid_log_partitions": ed_grid_logz,
        "class_reconstruction": class_table,
        "summary": pd.DataFrame(rows),
    }


def quadrature_convergence_report(
    spec: TFIMSpec,
    beta_max: float,
    *,
    point_counts: Sequence[int] = (5, 9, 15, 29),
    method: str = "simpson",
) -> pd.DataFrame:
    """Measure integration-only error using ED energies and exact ED log Z."""

    if beta_max <= 0:
        raise ValueError("beta_max must be positive.")
    rows = []
    dimension = 2**spec.num_system_qubits
    for points in point_counts:
        if int(points) < 2:
            raise ValueError("Every point count must be at least two.")
        beta_values = np.linspace(0.0, float(beta_max), int(points))
        exact = exact_class_thermodynamics(spec, beta_values)
        reconstructed = reconstruct_class_log_partitions(
            exact,
            dimension=dimension,
            energy_column="thermal_energy",
            method=method,
        )
        comparison = exact.merge(
            reconstructed[
                ["beta", "class_id", "log_partition_reconstructed"]
            ],
            on=["beta", "class_id"],
        )
        error = (
            comparison["log_partition_reconstructed"]
            - comparison["log_partition_function"]
        )
        endpoint = comparison.loc[np.isclose(comparison["beta"], beta_max)]
        endpoint_error = (
            endpoint["log_partition_reconstructed"]
            - endpoint["log_partition_function"]
        )
        rows.append(
            {
                "beta_max": float(beta_max),
                "points": int(points),
                "beta_spacing": float(beta_max) / (int(points) - 1),
                "integration_method": method,
                "max_abs_log_partition_error_all_beta": float(
                    np.max(np.abs(error))
                ),
                "max_abs_log_partition_error_at_beta_max": float(
                    np.max(np.abs(endpoint_error))
                ),
            }
        )
    return pd.DataFrame(rows)
