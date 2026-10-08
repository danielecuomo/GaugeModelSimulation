"""Execute circuits and decode every measured gauge record.

The lower-level backend helpers in :mod:`gauge_ite.backends` remain useful
for numerical studies.  This module provides the user-facing execution
object: it keeps the circuit that was built, the circuit that was measured,
the circuit that was actually sent to the simulator, and the resulting
counts in one place.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np
import pandas as pd
from qiskit import QuantumCircuit

from .branches import record_from_bits
from .circuits import CircuitBundle


@dataclass(frozen=True)
class CircuitExecution:
    """Complete provenance for one finite-shot circuit execution."""

    bundle: CircuitBundle
    basis: str
    shots: int
    seed: int | None
    measured_circuit: QuantumCircuit
    transpiled_circuit: QuantumCircuit | None
    counts: dict[str, int]
    backend_name: str
    provider: str = "aer"
    job_id: str | None = None
    register_counts: dict[str, dict[str, int]] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def logical_circuit(self) -> QuantumCircuit:
        return self.bundle.circuit


def _record_bits_from_displayed_key(
    displayed_key: str,
    num_record_bits: int,
) -> tuple[int, ...]:
    """Recover little-endian record bits from a Qiskit count key.

    Record bits occupy the least-significant classical positions in every
    circuit built by this package.  Qiskit displays those bits on the right
    and reverses their order for human-readable count strings.
    """

    compact = str(displayed_key).replace(" ", "").replace("_", "")
    if compact.startswith(("0x", "0X")):
        displayed_record = format(int(compact, 16), f"0{num_record_bits}b")[-num_record_bits:]
    else:
        if not compact or set(compact) - {"0", "1"}:
            raise ValueError(f"Unsupported count key: {displayed_key!r}.")
        if len(compact) < num_record_bits:
            compact = compact.zfill(num_record_bits)
        displayed_record = compact[-num_record_bits:]
    return tuple(int(bit) for bit in reversed(displayed_record))


def decode_counts(
    bundle: CircuitBundle,
    counts: Mapping[str, int],
) -> dict[str, object]:
    """Convert raw counts into branch, gauge-class, and summary tables."""

    normalized_counts = {str(key): int(value) for key, value in counts.items()}
    total_shots = int(sum(normalized_counts.values()))
    if total_shots <= 0:
        raise ValueError("counts must contain at least one shot.")

    rows: list[dict[str, object]] = []
    for displayed_key, count in normalized_counts.items():
        bits = _record_bits_from_displayed_key(
            displayed_key,
            bundle.layout.num_record_bits,
        )
        record = record_from_bits(
            bundle.spec,
            bits,
            steps=bundle.protocol.steps,
            architecture=bundle.protocol.architecture,
        )
        rows.append(
            {
                "beta": bundle.protocol.beta,
                "architecture": bundle.protocol.architecture,
                "steps": bundle.protocol.steps,
                "displayed_key": displayed_key,
                "ancilla_bitstring": record.ancilla_bitstring,
                "measured_signs_by_step": record.measured_signs_by_step,
                "effective_bond_signs_by_step": record.effective_bond_signs_by_step,
                "field_signs_by_step": record.field_signs_by_step,
                "class_ids_by_step": record.class_ids_by_step,
                "class_id": record.class_id,
                "cycle_fluxes_by_step": record.cycle_fluxes_by_step,
                "cycle_fluxes": record.cycle_fluxes,
                "plaquette_fluxes_by_step": record.plaquette_fluxes_by_step,
                "plaquette_fluxes": record.plaquette_fluxes,
                "flux_label": bundle.spec.lattice.format_fluxes(record.cycle_fluxes),
                "persistent_signs": (
                    record.persistent_signs
                    if bundle.protocol.architecture != "coherent_reuse"
                    else None
                ),
                "count": count,
                "probability": count / total_shots,
            }
        )

    branches = (
        pd.DataFrame(rows)
        .sort_values(["count", "ancilla_bitstring"], ascending=[False, True])
        .reset_index(drop=True)
    )
    grouping = [
        "beta",
        "architecture",
        "steps",
        "class_id",
        "cycle_fluxes",
        "plaquette_fluxes",
        "flux_label",
    ]
    classes = (
        branches.groupby(grouping, as_index=False, sort=True)
        .agg(
            count=("count", "sum"),
            probability=("probability", "sum"),
            observed_records=("ancilla_bitstring", "count"),
        )
        .sort_values("class_id")
        .reset_index(drop=True)
    )
    classes["probability_stderr"] = np.sqrt(
        classes["probability"] * (1.0 - classes["probability"]) / total_shots
    )

    if bundle.protocol.architecture == "coherent_reuse":
        persistent_probability = None
    else:
        persistent_probability = float(
            branches.loc[
                branches["persistent_signs"].fillna(False).astype(bool),
                "probability",
            ].sum()
        )

    summary = {
        "shots": total_shots,
        "unique_records": int(len(branches)),
        "observed_classes": int(len(classes)),
        "probability_sum": float(branches["probability"].sum()),
        "persistent_probability": persistent_probability,
        "persistence_is_defined": bundle.protocol.architecture != "coherent_reuse",
    }
    return {
        "counts": normalized_counts,
        "branches": branches,
        "classes": classes,
        "summary": summary,
    }


def run_aer_experiment(
    bundle: CircuitBundle,
    basis: str = "ancilla",
    shots: int = 4096,
    seed: int | None = 2026,
    backend=None,
    optimization_level: int = 0,
) -> CircuitExecution:
    """Backward-compatible convenience wrapper around :class:`AerExecutor`.

    New code can use :class:`gauge_ite.executors.AerExecutor` directly.  The
    import is local so importing the core package never opens a provider
    session and keeps backend dependencies isolated from the physics code.
    """

    from .executors import AerExecutor

    return AerExecutor(
        shots=shots,
        seed=seed,
        backend=backend,
        optimization_level=optimization_level,
    ).run(bundle, basis=basis)


__all__ = ["CircuitExecution", "decode_counts", "run_aer_experiment"]
