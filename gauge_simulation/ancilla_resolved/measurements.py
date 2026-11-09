"""Provider-neutral reduction of finite-shot measurement counts."""

from __future__ import annotations

from math import sqrt
from typing import Mapping

import numpy as np
import pandas as pd

from .branches import record_from_bits
from .circuits import CircuitBundle


def displayed_bits_little_endian(displayed_key: str, total_bits: int) -> tuple[int, ...]:
    """Decode a Qiskit count key into classical-bit order ``c[0], c[1], ...``."""

    compact = str(displayed_key).replace(" ", "").replace("_", "")
    if compact.startswith(("0x", "0X")):
        compact = format(int(compact, 16), f"0{total_bits}b")
    elif not compact or set(compact) - {"0", "1"}:
        raise ValueError(f"Unsupported count key: {displayed_key!r}.")
    if len(compact) < total_bits:
        compact = compact.zfill(total_bits)
    if len(compact) > total_bits:
        raise ValueError(
            f"Count key {displayed_key!r} contains more than {total_bits} bits."
        )
    return tuple(int(bit) for bit in reversed(compact))


def joint_count_table(
    bundle: CircuitBundle,
    counts: Mapping[str, int],
    basis: str,
) -> pd.DataFrame:
    """Convert joint system-record counts into per-outcome energy samples."""

    basis = basis.lower()
    if basis not in {"z", "x"}:
        raise ValueError("basis must be 'z' or 'x'.")
    n_record = bundle.layout.num_record_bits
    n_system = bundle.spec.num_system_qubits
    total_bits = n_record + n_system
    rows = []
    for displayed_key, raw_count in counts.items():
        count = int(raw_count)
        if count <= 0:
            continue
        bits = displayed_bits_little_endian(displayed_key, total_bits)
        record_bits = bits[:n_record]
        system_bits = bits[n_record:]
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
                "displayed_key": str(displayed_key),
                "record_bits": record_bits,
                "system_bits": system_bits,
                "count": count,
                "energy_sample": float(value),
            }
        )
    if not rows:
        raise ValueError("counts must contain at least one positive count.")
    return pd.DataFrame(rows)


def _summarize_basis(table: pd.DataFrame, value_name: str) -> pd.DataFrame:
    rows = []
    total_shots = int(table["count"].sum())
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
                "class_id": int(keys[0]),
                "cycle_fluxes": keys[1],
                "flux_label": keys[2],
                f"shots_{value_name}": n,
                f"probability_{value_name}": n / total_shots,
                value_name: mean,
                f"{value_name}_stderr": (
                    sqrt(max(variance, 0.0) / n) if n > 1 else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def conditioned_energy_from_counts(
    bundle: CircuitBundle,
    z_counts: Mapping[str, int],
    x_counts: Mapping[str, int],
) -> dict[str, object]:
    """Estimate class-resolved energy from paired Z- and X-basis runs.

    For a multi-step coherent-reuse circuit, this is explicitly an
    operational energy conditioned on the final record; it is not assigned to
    a static Hamiltonian repeated at every step.
    """

    z_rows = joint_count_table(bundle, z_counts, "z")
    x_rows = joint_count_table(bundle, x_counts, "x")
    z_summary = _summarize_basis(z_rows, "zz_energy")
    x_summary = _summarize_basis(x_rows, "x_energy")
    classes = z_summary.merge(
        x_summary,
        on=["class_id", "cycle_fluxes", "flux_label"],
        how="outer",
    ).sort_values("class_id").reset_index(drop=True)
    classes["probability"] = (
        classes["probability_zz_energy"] + classes["probability_x_energy"]
    ) / 2
    classes["energy_estimate"] = classes["zz_energy"] + classes["x_energy"]
    classes["energy_stderr"] = np.sqrt(
        classes["zz_energy_stderr"] ** 2 + classes["x_energy_stderr"] ** 2
    )

    z_total = int(z_rows["count"].sum())
    x_total = int(x_rows["count"].sum())
    direct = float(
        np.sum(z_rows["count"] * z_rows["energy_sample"]) / z_total
        + np.sum(x_rows["count"] * x_rows["energy_sample"]) / x_total
    )
    valid = classes.dropna(subset=["energy_estimate", "probability"])
    reconstructed = float(np.sum(valid["probability"] * valid["energy_estimate"]))
    return {
        "class_energies": classes,
        "z_samples": z_rows,
        "x_samples": x_rows,
        "direct_unconditioned_energy": direct,
        "class_reconstructed_energy": reconstructed,
        "absolute_closure_error": abs(direct - reconstructed),
        "interpretation": bundle.protocol.interpretation,
    }


__all__ = [
    "conditioned_energy_from_counts",
    "displayed_bits_little_endian",
    "joint_count_table",
]
