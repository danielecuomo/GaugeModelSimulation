"""Kraus operators and operational Gibbs-state certification."""

from __future__ import annotations

from itertools import product
from math import sqrt
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from qiskit.quantum_info import Statevector
from scipy.linalg import expm

from .branches import BranchRecord, record_from_bits, record_from_signs
from .circuits import build_tfim_circuit
from .models import TFIMSpec, local_pauli_unitary, pauli_term_matrix
from .protocols import ITEProtocol
from .circuits.primitive import ite_angle
from .thermodynamics import class_representatives


def elementary_kraus(
    spec: TFIMSpec,
    term_index: int,
    measured_sign: int,
    protocol: ITEProtocol,
) -> np.ndarray:
    """Exact conditional map of one ancilla-Pauli interaction."""

    if measured_sign not in (-1, 1):
        raise ValueError("measured_sign must be +1 or -1.")
    term = spec.terms[term_index]
    angle = ite_angle(
        protocol.beta,
        protocol.steps,
        term.signed_coefficient,
        protocol.angle_mapping,
    )
    dimension = 2**spec.num_system_qubits
    pauli = pauli_term_matrix(term.kind, term.sites, spec.num_system_qubits)
    return (
        np.cos(angle) * np.eye(dimension)
        + measured_sign * np.sin(angle) * pauli
    ) / sqrt(2)


def trajectory_kraus(
    spec: TFIMSpec,
    protocol: ITEProtocol,
    record: BranchRecord,
) -> np.ndarray:
    """Kraus map for one layer or a fully recorded fresh-ancilla trajectory."""

    if protocol.architecture == "coherent_reuse" and protocol.steps > 1:
        raise ValueError(
            "A coherent-reuse final record contains hidden histories. Use "
            "extract_all_coherent_kraus instead."
        )
    if len(record.measured_signs_by_step) != protocol.steps:
        if not (protocol.steps == 1 and len(record.measured_signs_by_step) == 1):
            raise ValueError("The record does not match the protocol step count.")
    dimension = 2**spec.num_system_qubits
    result = np.eye(dimension, dtype=complex)
    for signs in record.measured_signs_by_step:
        for term_index, sign in enumerate(signs):
            result = elementary_kraus(spec, term_index, sign, protocol) @ result
    return result


def operational_class_energy_sweep(
    spec: TFIMSpec,
    beta_values: Sequence[float],
    *,
    angle_mapping: str = "exact_kraus",
    shots_per_basis: int | None = None,
) -> pd.DataFrame:
    """Evaluate one noiseless one-step operational branch per gauge class.

    Gauge covariance makes the energy and branch probability identical for
    all microscopic records in one class after transformation to a common
    representative frame.  Evaluating one representative is therefore
    equivalent to pooling the full noiseless circuit distribution, while
    avoiding a separate large-register simulation at every beta value.

    ``shots_per_basis`` does not draw random samples.  When supplied, the
    returned standard error estimates use the exact conditional variances and
    the expected number of class records in independent Z- and X-basis runs.
    The beta-zero value is used as the known thermodynamic-integration anchor
    and is assigned zero sampling error.
    """

    beta_values = np.asarray(beta_values, dtype=float)
    if beta_values.ndim != 1 or len(beta_values) == 0:
        raise ValueError("beta_values must be a non-empty one-dimensional sequence.")
    if np.any(beta_values < 0):
        raise ValueError("beta_values must be non-negative.")
    if len(np.unique(beta_values)) != len(beta_values):
        raise ValueError("beta_values must not contain duplicates.")
    if shots_per_basis is not None and shots_per_basis <= 0:
        raise ValueError("shots_per_basis must be positive when supplied.")

    representatives = class_representatives(spec)
    dimension = 2**spec.num_system_qubits
    rows = []
    for beta in np.sort(beta_values):
        protocol = ITEProtocol(
            beta=float(beta),
            steps=1,
            angle_mapping=angle_mapping,
            architecture="coherent_reuse",
            input_mode="zero",
        )
        for class_row in representatives.itertuples():
            measured_bonds = tuple(
                effective * compiled
                for effective, compiled in zip(
                    class_row.representative_bond_signs,
                    spec.compiled_bond_signs,
                )
            )
            measured_signs = measured_bonds + tuple(
                class_row.representative_field_signs
            )
            record = record_from_signs(spec, (measured_signs,))
            kraus = trajectory_kraus(spec, protocol, record)
            omega = kraus @ kraus.conj().T
            omega = (omega + omega.conj().T) / 2
            omega_trace = float(np.real(np.trace(omega)))
            state = omega / omega_trace

            zz_hamiltonian = np.zeros((dimension, dimension), dtype=complex)
            for coefficient, sign, edge in zip(
                spec.gamma,
                record.bond_signs,
                spec.lattice.edges,
            ):
                zz_hamiltonian += (
                    coefficient
                    * sign
                    * pauli_term_matrix("ZZ", edge, spec.num_system_qubits)
                )
            x_hamiltonian = np.zeros((dimension, dimension), dtype=complex)
            for vertex, (coefficient, sign) in enumerate(
                zip(spec.eta, record.field_signs)
            ):
                x_hamiltonian += (
                    coefficient
                    * sign
                    * pauli_term_matrix("X", (vertex,), spec.num_system_qubits)
                )
            hamiltonian = zz_hamiltonian + x_hamiltonian

            def expectation(operator: np.ndarray) -> float:
                return float(np.real(np.trace(operator @ state)))

            zz_energy = expectation(zz_hamiltonian)
            x_energy = expectation(x_hamiltonian)
            zz_variance = max(
                expectation(zz_hamiltonian @ zz_hamiltonian) - zz_energy**2,
                0.0,
            )
            x_variance = max(
                expectation(x_hamiltonian @ x_hamiltonian) - x_energy**2,
                0.0,
            )
            rows.append(
                {
                    "beta": float(beta),
                    "class_id": int(class_row.class_id),
                    "cycle_fluxes": class_row.cycle_fluxes,
                    "flux_label": class_row.flux_label,
                    "representative_bond_signs": (
                        class_row.representative_bond_signs
                    ),
                    "full_sign_multiplicity": int(
                        class_row.full_sign_multiplicity
                    ),
                    "representative_branch_probability": omega_trace / dimension,
                    "unnormalized_class_probability": (
                        class_row.full_sign_multiplicity * omega_trace / dimension
                    ),
                    "zz_energy": zz_energy,
                    "x_energy": x_energy,
                    "energy_estimate": zz_energy + x_energy,
                    "zz_energy_variance": zz_variance,
                    "x_energy_variance": x_variance,
                    "ground_energy": float(
                        np.linalg.eigvalsh(hamiltonian)[0]
                    ),
                }
            )

    table = pd.DataFrame(rows).sort_values(["beta", "class_id"]).reset_index(
        drop=True
    )
    table["probability"] = table["unnormalized_class_probability"] / table.groupby(
        "beta"
    )["unnormalized_class_probability"].transform("sum")
    table["expected_class_shots_per_basis"] = (
        np.nan
        if shots_per_basis is None
        else shots_per_basis * table["probability"]
    )
    if shots_per_basis is None:
        table["energy_stderr"] = np.nan
    else:
        table["energy_stderr"] = np.sqrt(
            (table["zz_energy_variance"] + table["x_energy_variance"])
            / table["expected_class_shots_per_basis"]
        )
        table.loc[np.isclose(table["beta"], 0.0), "energy_stderr"] = 0.0
    return table


def all_one_step_kraus(
    spec: TFIMSpec,
    beta: float,
    angle_mapping: str = "exact_kraus",
) -> dict[BranchRecord, np.ndarray]:
    protocol = ITEProtocol(
        beta=beta,
        steps=1,
        angle_mapping=angle_mapping,
        architecture="coherent_reuse",
        input_mode="zero",
    )
    output = {}
    for bits in product((0, 1), repeat=spec.num_terms):
        record = record_from_bits(spec, bits, steps=1, architecture="coherent_reuse")
        output[record] = trajectory_kraus(spec, protocol, record)
    return output


def persistent_trajectory_kraus(
    spec: TFIMSpec,
    beta: float,
    steps: int,
    angle_mapping: str = "exact_kraus",
) -> dict[BranchRecord, np.ndarray]:
    """Return the subset of fresh trajectories whose signs persist."""

    protocol = ITEProtocol(
        beta=beta,
        steps=steps,
        angle_mapping=angle_mapping,
        architecture="fresh_trajectory",
        input_mode="zero",
    )
    output = {}
    for signs in product((1, -1), repeat=spec.num_terms):
        record = record_from_signs(spec, (signs,) * steps)
        output[record] = trajectory_kraus(spec, protocol, record)
    return output


def extract_all_coherent_kraus(
    spec: TFIMSpec,
    beta: float,
    steps: int,
    angle_mapping: str = "exact_kraus",
    max_evolution_qubits: int = 16,
) -> dict[BranchRecord, np.ndarray]:
    """Extract every final-record Kraus map of the coherent legacy circuit.

    The extraction simulates each system computational-basis input once and
    slices the resulting statevector by the final ancilla record. It is meant
    for small certification lattices such as 2x2.
    """

    evolution_qubits = spec.num_system_qubits + spec.num_terms
    if evolution_qubits > max_evolution_qubits:
        raise ValueError(
            f"Coherent Kraus extraction is limited to {max_evolution_qubits} "
            f"system-plus-ancilla qubits; received {evolution_qubits}."
        )
    n_system = spec.num_system_qubits
    n_record = spec.num_terms
    system_dimension = 2**n_system
    record_dimension = 2**n_record
    operators = np.zeros(
        (record_dimension, system_dimension, system_dimension), dtype=complex
    )

    for input_index in range(system_dimension):
        bits = tuple((input_index >> qubit) & 1 for qubit in range(n_system))
        protocol = ITEProtocol(
            beta=beta,
            steps=steps,
            angle_mapping=angle_mapping,
            architecture="coherent_reuse",
            input_mode="basis",
            basis_bits=bits,
        )
        circuit = build_tfim_circuit(spec, protocol).circuit
        state = np.asarray(Statevector.from_instruction(circuit).data)
        block = state.reshape((record_dimension, system_dimension))
        operators[:, :, input_index] = block

    output = {}
    for record_index in range(record_dimension):
        bits = tuple((record_index >> bit) & 1 for bit in range(n_record))
        record = record_from_bits(
            spec, bits, steps=steps, architecture="coherent_reuse"
        )
        output[record] = operators[record_index]
    return output


def _trace_norm_hermitian(matrix: np.ndarray) -> float:
    matrix = (matrix + matrix.conj().T) / 2
    return float(np.sum(np.abs(np.linalg.eigvalsh(matrix))))


def _operational_hamiltonian(omega: np.ndarray, beta: float) -> np.ndarray:
    omega = (omega + omega.conj().T) / 2
    eigenvalues, eigenvectors = np.linalg.eigh(omega)
    scale = max(float(np.max(eigenvalues)), 1.0)
    clipped = np.clip(eigenvalues, np.finfo(float).eps * scale, None)
    return eigenvectors @ np.diag(-np.log(clipped) / beta) @ eigenvectors.conj().T


def _procrustes_relative_error(
    kraus: np.ndarray,
    target_half: np.ndarray,
    gamma: float,
) -> float:
    scaled = sqrt(max(gamma, 0.0)) * target_half
    overlap = scaled.conj().T @ kraus
    nuclear = float(np.sum(np.linalg.svd(overlap, compute_uv=False)))
    minimum_sq = max(
        np.linalg.norm(kraus, "fro") ** 2
        + np.linalg.norm(scaled, "fro") ** 2
        - 2 * nuclear,
        0.0,
    )
    return sqrt(minimum_sq) / max(float(np.linalg.norm(kraus, "fro")), 1e-15)


def certify_kraus_map(
    spec: TFIMSpec,
    beta: float,
    kraus_map: Mapping[BranchRecord, np.ndarray],
    label: str,
    require_complete: bool = True,
) -> dict:
    """Compare exact operational branches with intended branch Gibbs states."""

    if beta <= 0:
        raise ValueError("Operational Hamiltonian certification requires beta > 0.")
    dimension = 2**spec.num_system_qubits
    rows = []
    completeness = np.zeros((dimension, dimension), dtype=complex)
    omega_by_record: dict[BranchRecord, np.ndarray] = {}

    for record, kraus in kraus_map.items():
        omega = kraus @ kraus.conj().T
        omega = (omega + omega.conj().T) / 2
        omega_by_record[record] = omega
        completeness += kraus.conj().T @ kraus
        hamiltonian = spec.branch_hamiltonian(record.bond_signs, record.field_signs)
        target_omega = expm(-beta * hamiltonian)
        target_half = expm(-beta * hamiltonian / 2)
        partition = float(np.real(np.trace(target_omega)))
        omega_trace = float(np.real(np.trace(omega)))
        gamma = omega_trace / partition
        delta = _trace_norm_hermitian(omega - gamma * target_omega)
        normalized_distance = 0.5 * _trace_norm_hermitian(
            omega / omega_trace - target_omega / partition
        )
        hop = _operational_hamiltonian(omega, beta)
        scalar_shift = float(np.real(np.trace(hop - hamiltonian)) / dimension)
        hop_difference = hop - hamiltonian - scalar_shift * np.eye(dimension)
        hop_relative = float(
            np.linalg.norm(hop_difference, "fro")
            / max(np.linalg.norm(hamiltonian, "fro"), 1e-15)
        )
        rows.append(
            {
                "label": label,
                "ancilla_bitstring": record.ancilla_bitstring,
                "class_id": record.class_id,
                "cycle_fluxes": record.cycle_fluxes,
                "flux_label": spec.lattice.format_fluxes(record.cycle_fluxes),
                "bond_signs": record.bond_signs,
                "field_signs": record.field_signs,
                "persistent_signs": record.persistent_signs,
                "branch_probability": omega_trace / dimension,
                "partition_function": partition,
                "gamma_prefactor": gamma,
                "alpha_magnitude": sqrt(max(gamma, 0.0)),
                "delta_trace_norm": delta,
                "state_trace_distance": normalized_distance,
                "operational_H_relative_error": hop_relative,
                "operational_scalar_shift": scalar_shift,
                "kraus_right_unitary_error": _procrustes_relative_error(
                    kraus, target_half, gamma
                ),
            }
        )

    branches = pd.DataFrame(rows)
    coverage = float(branches["branch_probability"].sum())
    branches["normalized_branch_probability"] = (
        branches["branch_probability"] / coverage
    )
    target_total = float(branches["partition_function"].sum())
    branches["ideal_normalized_probability"] = (
        branches["partition_function"] / target_total
    )
    tv_distance = 0.5 * float(
        np.abs(
            branches["normalized_branch_probability"]
            - branches["ideal_normalized_probability"]
        ).sum()
    )
    completeness_error = float(
        np.linalg.norm(completeness - np.eye(dimension), ord=2)
    )

    class_table = (
        branches.groupby(
            ["class_id", "cycle_fluxes", "flux_label"],
            as_index=False,
            sort=True,
        )
        .agg(
            probability=("normalized_branch_probability", "sum"),
            unnormalized_probability=("branch_probability", "sum"),
            ideal_probability=("ideal_normalized_probability", "sum"),
            num_branches=("ancilla_bitstring", "count"),
            mean_state_trace_distance=("state_trace_distance", "mean"),
            max_state_trace_distance=("state_trace_distance", "max"),
        )
    )
    summary = pd.DataFrame(
        [
            {
                "label": label,
                "num_branches": len(branches),
                "probability_coverage": coverage,
                "completeness_error": completeness_error,
                "completeness_expected": bool(require_complete),
                "total_variation_to_intended_Gibbs": tv_distance,
                "max_state_trace_distance": float(
                    branches["state_trace_distance"].max()
                ),
                "mean_state_trace_distance": float(
                    branches["state_trace_distance"].mean()
                ),
                "max_operational_H_relative_error": float(
                    branches["operational_H_relative_error"].max()
                ),
                "max_kraus_right_unitary_error": float(
                    branches["kraus_right_unitary_error"].max()
                ),
                "relative_alpha_spread": float(
                    branches["alpha_magnitude"].std(ddof=0)
                    / max(branches["alpha_magnitude"].mean(), 1e-15)
                ),
            }
        ]
    )
    return {
        "branches": branches,
        "classes": class_table,
        "summary": summary,
        "kraus": dict(kraus_map),
        "omega": omega_by_record,
    }


def operational_gauge_covariance_report(
    spec: TFIMSpec,
    certification: dict,
) -> pd.DataFrame:
    """Compare gauge-corrected operational states within every class."""

    branches = certification["branches"]
    omega_map: Mapping[BranchRecord, np.ndarray] = certification["omega"]
    record_by_bits = {record.ancilla_bitstring: record for record in omega_map}
    rows = []
    for class_id, group in branches.groupby("class_id", sort=True):
        representative_row = group.sort_values("ancilla_bitstring").iloc[0]
        representative = record_by_bits[representative_row.ancilla_bitstring]
        rep_omega = omega_map[representative]
        rep_state = rep_omega / np.trace(rep_omega)
        rep_probability = float(representative_row.branch_probability)
        for row in group.itertuples():
            record = record_by_bits[row.ancilla_bitstring]
            vertex_gauge = spec.lattice.find_vertex_gauge(
                record.bond_signs, representative.bond_signs
            )
            field_gauge = tuple(
                source * target
                for source, target in zip(
                    record.field_signs, representative.field_signs
                )
            )
            unitary = local_pauli_unitary(vertex_gauge, field_gauge)
            omega = omega_map[record]
            corrected = unitary @ (omega / np.trace(omega)) @ unitary.conj().T
            rows.append(
                {
                    "class_id": int(class_id),
                    "flux_label": representative_row.flux_label,
                    "representative_bitstring": representative.ancilla_bitstring,
                    "ancilla_bitstring": record.ancilla_bitstring,
                    "gauge_corrected_trace_distance": 0.5
                    * _trace_norm_hermitian(corrected - rep_state),
                    "relative_probability_difference": abs(
                        float(row.branch_probability) - rep_probability
                    )
                    / max(rep_probability, 1e-15),
                }
            )
    return pd.DataFrame(rows)
