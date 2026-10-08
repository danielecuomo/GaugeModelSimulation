"""Save readable and machine-readable circuit artifacts."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

from qiskit import QuantumCircuit, qpy


def circuit_statistics(circuit: QuantumCircuit) -> dict[str, object]:
    """Return a JSON-friendly summary of one circuit."""

    return {
        "num_qubits": int(circuit.num_qubits),
        "num_clbits": int(circuit.num_clbits),
        "depth": int(circuit.depth() or 0),
        "size": int(circuit.size()),
        "operations": {
            str(name): int(count) for name, count in circuit.count_ops().items()
        },
    }


def save_circuit(
    circuit: QuantumCircuit,
    stem: str | Path,
    *,
    save_mpl: bool = True,
    fold: int = 160,
) -> dict[str, Path | str]:
    """Save text, QPY, and optional PNG versions of a circuit."""

    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    text_path = stem.with_suffix(".txt")
    qpy_path = stem.with_suffix(".qpy")
    text_path.write_text(
        str(circuit.draw(output="text", fold=fold)),
        encoding="utf-8",
    )
    with qpy_path.open("wb") as handle:
        qpy.dump(circuit, handle)

    result: dict[str, Path | str] = {"text": text_path, "qpy": qpy_path}
    if save_mpl:
        try:
            import matplotlib.pyplot as plt

            figure = circuit.draw(output="mpl", fold=fold, scale=0.65)
            png_path = stem.with_suffix(".png")
            figure.savefig(png_path, dpi=220, bbox_inches="tight")
            plt.close(figure)
            result["png"] = png_path
        except Exception as exc:  # pragma: no cover - depends on optional renderer
            result["png_warning"] = f"{type(exc).__name__}: {exc}"
    return result


def save_circuit_set(
    circuits: Mapping[str, QuantumCircuit],
    output_dir: str | Path,
    *,
    save_mpl: bool = True,
) -> dict[str, dict[str, Path | str]]:
    """Save several named circuits using a consistent file convention."""

    output = Path(output_dir)
    return {
        name: save_circuit(circuit, output / name, save_mpl=save_mpl)
        for name, circuit in circuits.items()
    }


__all__ = ["circuit_statistics", "save_circuit", "save_circuit_set"]
