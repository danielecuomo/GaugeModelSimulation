"""Aer execution helpers that preserve every measured branch."""

from __future__ import annotations

import warnings
from math import sqrt
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from qiskit import transpile
from qiskit_aer import AerSimulator
from scipy.sparse import SparseEfficiencyWarning

from ..branches import record_from_bits
from ..circuits import CircuitBundle, add_final_measurements, build_tfim_circuit
from ..models import TFIMSpec
from ..protocols import ITEProtocol


def _key_to_int(key) -> int:
    if isinstance(key, str):
        compact = key.replace(" ", "")
        if compact.startswith(("0x", "0X")):
            return int(compact, 16)
        if compact and set(compact) <= {"0", "1"}:
            return int(compact, 2)
        return int(compact, 0)
    return int(key)


def exact_probability_dict(
    bundle: CircuitBundle,
    basis: str = "ancilla",
    backend=None,
) -> dict[int, float]:
    """Return an exact ideal probability dictionary from Aer MPS."""

    basis = basis.lower()
    if basis not in {"ancilla", "z", "x"}:
        raise ValueError("basis must be 'ancilla', 'z', or 'x'.")
    if bundle.protocol.architecture == "measure_reset_trajectory":
        raise ValueError("Exact deferred probabilities do not support mid-circuit measurement.")

    circuit = bundle.circuit.copy()
    qubits = list(bundle.layout.record_qubits)
    if basis in {"z", "x"}:
        if basis == "x":
            for qubit in bundle.layout.system_qubits:
                circuit.h(qubit)
        qubits.extend(bundle.layout.system_qubits)
    label = f"{basis}_probabilities"
    circuit.save_probabilities_dict(qubits=qubits, label=label)
    backend = backend or AerSimulator(method="matrix_product_state")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)
        compiled = transpile(circuit, backend, optimization_level=0)
        raw = backend.run(compiled).result().data(0)[label]
    probabilities = {
        _key_to_int(key): float(probability)
        for key, probability in raw.items()
        if float(probability) > 0
    }
    total = sum(probabilities.values())
    return {key: value / total for key, value in probabilities.items()}


def sampled_counts(
    bundle: CircuitBundle,
    basis: str = "ancilla",
    shots: int = 4096,
    seed: int = 2026,
    backend=None,
) -> dict[str, int]:
    """Run a finite-shot Aer experiment, including dynamic trajectories."""

    if shots <= 0:
        raise ValueError("shots must be positive.")
    backend = backend or AerSimulator(method="matrix_product_state")
    measured = add_final_measurements(bundle, basis=basis)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)
        compiled = transpile(measured, backend, optimization_level=0)
        counts = backend.run(
            compiled,
            shots=shots,
            seed_simulator=seed,
        ).result().get_counts(0)
    return {str(key): int(value) for key, value in counts.items()}


def sampled_record_probabilities(
    bundle: CircuitBundle,
    shots: int = 4096,
    seed: int = 2026,
    backend=None,
) -> dict:
    """Decode finite-shot records using the canonical public table format."""

    counts = sampled_counts(
        bundle,
        basis="ancilla",
        shots=shots,
        seed=seed,
        backend=backend,
    )
    # Local import avoids a module-level cycle: the user-facing execution
    # module depends on the circuit package but not on this backend module.
    from ..execution import decode_counts

    return decode_counts(bundle, counts)


def _record_bits_from_index(index: int, num_bits: int) -> tuple[int, ...]:
    return tuple((index >> bit) & 1 for bit in range(num_bits))


def exact_class_probabilities(
    bundle: CircuitBundle,
    backend=None,
) -> dict:
    probabilities = exact_probability_dict(bundle, basis="ancilla", backend=backend)
    rows = []
    for record_index, probability in probabilities.items():
        bits = _record_bits_from_index(record_index, bundle.layout.num_record_bits)
        record = record_from_bits(
            bundle.spec,
            bits,
            steps=bundle.protocol.steps,
            architecture=bundle.protocol.architecture,
        )
        rows.append(
            {
                "beta": bundle.protocol.beta,
                "ancilla_bitstring": record.ancilla_bitstring,
                "class_id": record.class_id,
                "cycle_fluxes": record.cycle_fluxes,
                "plaquette_fluxes": record.plaquette_fluxes,
                "flux_label": bundle.spec.lattice.format_fluxes(record.cycle_fluxes),
                "persistent_signs": record.persistent_signs,
                "probability": probability,
            }
        )
    branches = pd.DataFrame(rows).sort_values("probability", ascending=False).reset_index(drop=True)
    classes = (
        branches.groupby(
            ["beta", "class_id", "cycle_fluxes", "plaquette_fluxes", "flux_label"],
            as_index=False,
            sort=True,
        )
        .agg(
            probability=("probability", "sum"),
            num_branches=("ancilla_bitstring", "count"),
            persistent_probability=(
                "probability",
                lambda values: float(
                    values[
                        branches.loc[values.index, "persistent_signs"].to_numpy(dtype=bool)
                    ].sum()
                ),
            ),
        )
        .sort_values("class_id")
        .reset_index(drop=True)
    )
    if bundle.protocol.architecture == "coherent_reuse":
        # A final coherent record is not a time-resolved trajectory, so the
        # notion of sign persistence does not apply.
        classes["persistent_probability"] = np.nan
    return {"branches": branches, "classes": classes, "probabilities": probabilities}


def probability_sweep(
    spec: TFIMSpec,
    beta_values: Sequence[float],
    steps: int = 1,
    angle_mapping: str = "exact_kraus",
    architecture: str = "coherent_reuse",
    input_mode: str = "purification",
    backend=None,
    return_branches: bool = False,
) -> dict:
    class_tables = []
    branch_tables = []
    bundles = []
    backend = backend or AerSimulator(method="matrix_product_state")
    for beta in beta_values:
        protocol = ITEProtocol(
            beta=float(beta),
            steps=steps,
            angle_mapping=angle_mapping,
            architecture=architecture,
            input_mode=input_mode,
        )
        bundle = build_tfim_circuit(spec, protocol)
        result = exact_class_probabilities(bundle, backend=backend)
        bundles.append(bundle)
        class_tables.append(result["classes"])
        if return_branches:
            branch_tables.append(result["branches"])
    output = {
        "classes": pd.concat(class_tables, ignore_index=True),
        "bundles": bundles,
    }
    if return_branches:
        output["branches"] = pd.concat(branch_tables, ignore_index=True)
    return output


def probability_diagnostics(
    exact_thermodynamics: pd.DataFrame,
    circuit_classes: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    reference = exact_thermodynamics[
        [
            "beta",
            "class_id",
            "flux_label",
            "annealed_class_weight",
            "multiplicity_weight",
        ]
    ]
    comparison = circuit_classes.merge(
        reference,
        on=["beta", "class_id", "flux_label"],
        how="inner",
    )
    comparison["annealed_probability_error"] = (
        comparison["probability"] - comparison["annealed_class_weight"]
    )
    comparison["uniform_probability_error"] = (
        comparison["probability"] - comparison["multiplicity_weight"]
    )
    rows = []
    for beta, group in comparison.groupby("beta", sort=True):
        rows.append(
            {
                "beta": float(beta),
                "total_variation_to_annealed": 0.5
                * float(np.abs(group["annealed_probability_error"]).sum()),
                "total_variation_to_uniform": 0.5
                * float(np.abs(group["uniform_probability_error"]).sum()),
                "max_abs_error_to_annealed": float(
                    np.abs(group["annealed_probability_error"]).max()
                ),
                "max_abs_error_to_uniform": float(
                    np.abs(group["uniform_probability_error"]).max()
                ),
            }
        )
    return comparison, pd.DataFrame(rows)


def sampled_class_probabilities(
    exact_classes: pd.DataFrame,
    shots: int = 4096,
    seed: int = 2026,
) -> pd.DataFrame:
    """Draw reproducible finite-shot class counts from ideal probabilities."""

    if shots <= 0:
        raise ValueError("shots must be positive.")
    rng = np.random.default_rng(seed)
    rows = []
    for beta, group in exact_classes.groupby("beta", sort=True):
        group = group.sort_values("class_id").reset_index(drop=True)
        probabilities = group["probability"].to_numpy(dtype=float).copy()
        probabilities = probabilities / probabilities.sum()
        counts = rng.multinomial(shots, probabilities)
        for row, count in zip(group.itertuples(), counts):
            probability = count / shots
            rows.append(
                {
                    "beta": float(beta),
                    "class_id": int(row.class_id),
                    "cycle_fluxes": row.cycle_fluxes,
                    "flux_label": row.flux_label,
                    "shots": int(count),
                    "probability": probability,
                    "probability_stderr": sqrt(
                        max(probability * (1 - probability), 0.0) / shots
                    ),
                }
            )
    return pd.DataFrame(rows)


def _joint_rows(
    bundle: CircuitBundle,
    probabilities: Mapping[int, float],
    basis: str,
) -> pd.DataFrame:
    n_record = bundle.layout.num_record_bits
    n_system = bundle.spec.num_system_qubits
    mask = (1 << n_record) - 1
    rows = []
    for joint_index, probability in probabilities.items():
        record_index = joint_index & mask
        system_index = joint_index >> n_record
        record_bits = _record_bits_from_index(record_index, n_record)
        system_bits = _record_bits_from_index(system_index, n_system)
        record = record_from_bits(
            bundle.spec,
            record_bits,
            steps=bundle.protocol.steps,
            architecture=bundle.protocol.architecture,
        )
        eigenvalues = tuple(1 if bit == 0 else -1 for bit in system_bits)
        if basis == "z":
            value = sum(
                coupling * sign * eigenvalues[u] * eigenvalues[v]
                for coupling, sign, (u, v) in zip(
                    bundle.spec.gamma,
                    record.bond_signs,
                    bundle.spec.lattice.edges,
                )
            )
        else:
            value = sum(
                coupling * sign * eigenvalues[vertex]
                for vertex, (coupling, sign) in enumerate(
                    zip(bundle.spec.eta, record.field_signs)
                )
            )
        rows.append(
            {
                "class_id": record.class_id,
                "cycle_fluxes": record.cycle_fluxes,
                "flux_label": bundle.spec.lattice.format_fluxes(record.cycle_fluxes),
                "record_index": record_index,
                "system_index": system_index,
                "probability": probability,
                "energy_sample": float(value),
            }
        )
    return pd.DataFrame(rows)


def exact_conditioned_energies(
    bundle: CircuitBundle,
    backend=None,
) -> dict:
    """Exact ideal class-conditioned energy using joint Z/X probabilities."""

    if bundle.protocol.steps != 1:
        raise ValueError(
            "Static branch-energy interpretation is enabled only for a one-step circuit."
        )
    backend = backend or AerSimulator(method="matrix_product_state")
    z_probabilities = exact_probability_dict(bundle, "z", backend=backend)
    x_probabilities = exact_probability_dict(bundle, "x", backend=backend)
    z_rows = _joint_rows(bundle, z_probabilities, "z")
    x_rows = _joint_rows(bundle, x_probabilities, "x")

    def summarize(table: pd.DataFrame, value_name: str) -> pd.DataFrame:
        rows = []
        for keys, group in table.groupby(
            ["class_id", "cycle_fluxes", "flux_label"], sort=True
        ):
            probability = float(group["probability"].sum())
            rows.append(
                {
                    "class_id": keys[0],
                    "cycle_fluxes": keys[1],
                    "flux_label": keys[2],
                    f"{value_name}_probability": probability,
                    value_name: float(
                        np.sum(group["probability"] * group["energy_sample"])
                        / probability
                    ),
                }
            )
        return pd.DataFrame(rows)

    z_summary = summarize(z_rows, "zz_energy")
    x_summary = summarize(x_rows, "x_energy")
    classes = z_summary.merge(
        x_summary,
        on=["class_id", "cycle_fluxes", "flux_label"],
    )
    classes["probability"] = (
        classes["zz_energy_probability"] + classes["x_energy_probability"]
    ) / 2
    classes["energy_estimate"] = classes["zz_energy"] + classes["x_energy"]
    direct_energy = float(
        np.sum(z_rows["probability"] * z_rows["energy_sample"])
        + np.sum(x_rows["probability"] * x_rows["energy_sample"])
    )
    reconstructed = float(np.sum(classes["probability"] * classes["energy_estimate"]))
    return {
        "class_energies": classes,
        "z_joint": z_rows,
        "x_joint": x_rows,
        "direct_unconditioned_energy": direct_energy,
        "class_reconstructed_energy": reconstructed,
        "absolute_closure_error": abs(direct_energy - reconstructed),
    }


def _sample_joint_table(
    table: pd.DataFrame,
    shots: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    probabilities = table["probability"].to_numpy(dtype=float)
    counts = rng.multinomial(shots, probabilities / probabilities.sum())
    sampled = table.copy()
    sampled["count"] = counts
    return sampled.loc[sampled["count"] > 0].reset_index(drop=True)


def sampled_conditioned_energies(
    exact_conditioned: dict,
    shots: int = 4096,
    seed: int = 2026,
) -> dict:
    """Draw ideal finite-shot samples from exact joint probabilities."""

    if shots <= 0:
        raise ValueError("shots must be positive.")
    rng = np.random.default_rng(seed)
    z = _sample_joint_table(exact_conditioned["z_joint"], shots, rng)
    x = _sample_joint_table(exact_conditioned["x_joint"], shots, rng)

    def summarize(table: pd.DataFrame, value_name: str) -> pd.DataFrame:
        rows = []
        for keys, group in table.groupby(
            ["class_id", "cycle_fluxes", "flux_label"], sort=True
        ):
            values = group["energy_sample"].to_numpy(dtype=float)
            weights = group["count"].to_numpy(dtype=float)
            n = int(weights.sum())
            mean = float(np.average(values, weights=weights))
            variance = (
                float(np.sum(weights * (values - mean) ** 2) / (n - 1))
                if n > 1
                else np.nan
            )
            rows.append(
                {
                    "class_id": keys[0],
                    "cycle_fluxes": keys[1],
                    "flux_label": keys[2],
                    f"shots_{value_name}": n,
                    value_name: mean,
                    f"{value_name}_stderr": sqrt(max(variance, 0.0) / n)
                    if n > 1
                    else np.nan,
                }
            )
        return pd.DataFrame(rows)

    z_summary = summarize(z, "zz_energy")
    x_summary = summarize(x, "x_energy")
    classes = z_summary.merge(
        x_summary,
        on=["class_id", "cycle_fluxes", "flux_label"],
        how="outer",
    )
    classes["energy_estimate"] = classes["zz_energy"] + classes["x_energy"]
    classes["energy_stderr"] = np.sqrt(
        classes["zz_energy_stderr"] ** 2 + classes["x_energy_stderr"] ** 2
    )
    return {"class_energies": classes, "z_samples": z, "x_samples": x}
