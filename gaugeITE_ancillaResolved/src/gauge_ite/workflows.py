"""High-level workflows used by the command line and examples."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib import metadata
import json
from pathlib import Path
import platform
from typing import Sequence

import numpy as np
import pandas as pd
from qiskit import transpile
from qiskit_aer import AerSimulator

from .backends import exact_class_probabilities, exact_conditioned_energies
from .circuit_io import circuit_statistics, save_circuit_set
from .circuits import CircuitBundle, add_final_measurements, build_tfim_circuit, circuit_resource_row
from .execution import CircuitExecution, decode_counts, run_aer_experiment
from .experiments import ExperimentResult
from .models import TFIMSpec
from .protocols import ITEProtocol
from .simulation import clean_ground_energy, clean_thermodynamics


@dataclass(frozen=True)
class RunArtifacts:
    """Paths and decoded data produced by :func:`save_experiment`."""

    output_dir: Path
    execution: CircuitExecution
    branches: pd.DataFrame
    classes: pd.DataFrame
    summary: dict[str, object]


def _package_version(name: str) -> str | None:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


def environment_report(backend_name: str | None = None) -> dict[str, object]:
    """Return the minimum software information needed to repeat a run."""

    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "backend": backend_name,
        "packages": {
            name: _package_version(name)
            for name in (
                "gauge-ite",
                "numpy",
                "scipy",
                "pandas",
                "matplotlib",
                "qiskit",
                "qiskit-aer",
                "qiskit-ibm-runtime",
            )
        },
    }


def experiment_config(execution: CircuitExecution) -> dict[str, object]:
    """Return a JSON-friendly description of an execution."""

    spec = execution.bundle.spec
    return {
        "model": "tfim",
        "lattice": {
            "rows": spec.lattice.rows,
            "cols": spec.lattice.cols,
            "boundary": spec.lattice.boundary,
            "edges": [list(edge) for edge in spec.lattice.edges],
        },
        "couplings": {
            "gamma": list(spec.gamma),
            "eta": list(spec.eta),
            "compiled_bond_signs": list(spec.compiled_bond_signs),
        },
        "protocol": asdict(execution.bundle.protocol),
        "measurement": {
            "basis": execution.basis,
            "shots": execution.shots,
            "seed": execution.seed,
        },
        "execution": {
            "provider": execution.provider,
            "backend": execution.backend_name,
            "job_id": execution.job_id,
        },
    }


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def save_experiment(
    execution: CircuitExecution,
    output_root: str | Path = "runs",
    run_name: str | None = None,
    *,
    save_mpl: bool = True,
) -> RunArtifacts:
    """Save one execution as a self-contained, inspectable run directory."""

    if run_name is None:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        spec = execution.bundle.spec
        protocol = execution.bundle.protocol
        beta_label = f"{protocol.beta:g}".replace(".", "p")
        run_name = (
            f"{stamp}_{spec.lattice.rows}x{spec.lattice.cols}_"
            f"{protocol.architecture}_s{protocol.steps}_beta{beta_label}"
        )
    output = Path(output_root) / run_name
    if output.exists():
        raise FileExistsError(f"Run directory already exists: {output}")
    output.mkdir(parents=True)

    decoded = decode_counts(execution.bundle, execution.counts)
    _write_json(output / "config.json", experiment_config(execution))
    _write_json(output / "environment.json", environment_report(execution.backend_name))
    _write_json(output / "counts.json", decoded["counts"])
    if execution.register_counts:
        _write_json(output / "register_counts.json", execution.register_counts)
    if execution.metadata:
        _write_json(output / "provider_metadata.json", execution.metadata)
    _write_json(output / "summary.json", decoded["summary"])

    decoded["branches"].to_csv(output / "branches.csv", index=False)
    decoded["classes"].to_csv(output / "classes.csv", index=False)
    pd.DataFrame([circuit_resource_row(execution.bundle.spec, execution.bundle.protocol)]).to_csv(
        output / "resource_summary.csv",
        index=False,
    )

    circuits = {
        "logical_circuit": execution.logical_circuit,
        "measured_circuit": execution.measured_circuit,
    }
    if execution.transpiled_circuit is not None:
        circuits["transpiled_circuit"] = execution.transpiled_circuit
    diagram_paths = save_circuit_set(circuits, output, save_mpl=save_mpl)
    circuit_report = {
        name: {
            "statistics": circuit_statistics(circuit),
            "files": {key: str(value) for key, value in diagram_paths[name].items()},
        }
        for name, circuit in circuits.items()
    }
    _write_json(output / "circuits.json", circuit_report)

    return RunArtifacts(
        output_dir=output,
        execution=execution,
        branches=decoded["branches"],
        classes=decoded["classes"],
        summary=decoded["summary"],
    )


def save_experiment_result(
    result: ExperimentResult,
    output_root: str | Path = "runs",
    run_name: str | None = None,
    *,
    save_mpl: bool = True,
) -> RunArtifacts | dict[str, object]:
    """Save a high-level result, while keeping the default path write-free.

    A one-basis class/trajectory run retains the original flat run layout.
    Energy estimation contains two independently submitted circuits, so its
    Z- and X-basis provenance is stored in sibling subdirectories.
    """

    if len(result.executions) == 1:
        return save_experiment(
            next(iter(result.executions.values())),
            output_root=output_root,
            run_name=run_name,
            save_mpl=save_mpl,
        )

    if run_name is None:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        protocol = result.bundle.protocol
        beta_label = f"{protocol.beta:g}".replace(".", "p")
        run_name = (
            f"{stamp}_{result.bundle.spec.lattice.rows}x"
            f"{result.bundle.spec.lattice.cols}_{protocol.architecture}_"
            f"s{protocol.steps}_beta{beta_label}_energy"
        )
    output = Path(output_root) / run_name
    if output.exists():
        raise FileExistsError(f"Run directory already exists: {output}")
    output.mkdir(parents=True)
    basis_artifacts = {
        basis: save_experiment(
            execution,
            output_root=output,
            run_name=basis,
            save_mpl=save_mpl,
        )
        for basis, execution in result.executions.items()
    }
    class_energies = result.analysis.get("class_energies")
    if class_energies is not None:
        class_energies.to_csv(output / "class_energies.csv", index=False)
    for key in ("z_samples", "x_samples"):
        table = result.analysis.get(key)
        if table is not None:
            table.to_csv(output / f"{key}.csv", index=False)
    summary = {
        "observable": result.config.observable,
        "interpretation": result.interpretation,
        "warnings": list(result.warnings),
        "job_ids": result.job_ids,
        "backend_names": result.backend_names,
        "direct_unconditioned_energy": result.analysis.get(
            "direct_unconditioned_energy"
        ),
        "class_reconstructed_energy": result.analysis.get(
            "class_reconstructed_energy"
        ),
        "absolute_closure_error": result.analysis.get("absolute_closure_error"),
    }
    _write_json(output / "summary.json", summary)
    return {
        "output_dir": output,
        "basis_artifacts": basis_artifacts,
        "result": result,
        "summary": summary,
    }


def run_and_save(
    spec: TFIMSpec,
    protocol: ITEProtocol,
    *,
    basis: str = "ancilla",
    shots: int = 4096,
    seed: int | None = 2026,
    optimization_level: int = 0,
    output_root: str | Path = "runs",
    run_name: str | None = None,
    save_mpl: bool = True,
) -> RunArtifacts:
    """Build, execute, decode, and save one complete Aer experiment."""

    bundle = build_tfim_circuit(spec, protocol)
    execution = run_aer_experiment(
        bundle,
        basis=basis,
        shots=shots,
        seed=seed,
        optimization_level=optimization_level,
    )
    return save_experiment(
        execution,
        output_root=output_root,
        run_name=run_name,
        save_mpl=save_mpl,
    )


def inspect_architectures(
    spec: TFIMSpec,
    *,
    beta: float = 0.75,
    steps: int = 2,
    angle_mapping: str = "exact_kraus",
    input_mode: str = "purification",
    basis_bits: tuple[int, ...] | None = None,
    output_dir: str | Path = "artifacts/circuits",
    optimization_level: int = 0,
    save_mpl: bool = True,
) -> pd.DataFrame:
    """Save side-by-side logical, measured, and Aer-transpiled architectures."""

    backend = AerSimulator(method="matrix_product_state")
    rows = []
    for architecture in (
        "coherent_reuse",
        "fresh_trajectory",
        "measure_reset_trajectory",
    ):
        protocol = ITEProtocol(
            beta=beta,
            steps=steps,
            angle_mapping=angle_mapping,
            architecture=architecture,
            input_mode=input_mode,
            basis_bits=basis_bits,
        )
        bundle = build_tfim_circuit(spec, protocol)
        measured = add_final_measurements(bundle, basis="ancilla")
        compiled = transpile(measured, backend=backend, optimization_level=optimization_level)
        architecture_dir = Path(output_dir) / architecture
        save_circuit_set(
            {
                "logical_circuit": bundle.circuit,
                "measured_circuit": measured,
                "transpiled_for_aer": compiled,
            },
            architecture_dir,
            save_mpl=save_mpl,
        )
        rows.append(
            {
                **circuit_resource_row(spec, protocol),
                "logical_depth": bundle.circuit.depth(),
                "transpiled_depth": compiled.depth(),
                "transpiled_two_qubit_gates": sum(
                    count
                    for name, count in compiled.count_ops().items()
                    if name in {"cx", "cz", "ecr", "rzz", "swap"}
                ),
                "artifact_directory": str(architecture_dir),
            }
        )
    table = pd.DataFrame(rows)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    table.to_csv(output / "architecture_resources.csv", index=False)
    return table


def quickstart_report(
    spec: TFIMSpec,
    protocol: ITEProtocol,
    beta_values: Sequence[float],
) -> dict[str, object]:
    """Return the simple energy summary that replaces the legacy driver."""

    thermodynamics = clean_thermodynamics(spec, beta_values)
    bundle = build_tfim_circuit(spec, protocol)
    probabilities = exact_class_probabilities(bundle)
    output: dict[str, object] = {
        "ground_energy": clean_ground_energy(spec),
        "clean_thermodynamics": thermodynamics,
        "class_probabilities": probabilities["classes"],
        "resources": pd.DataFrame([circuit_resource_row(spec, protocol)]),
    }
    if protocol.architecture != "measure_reset_trajectory":
        conditioned = exact_conditioned_energies(bundle)
        output["circuit_energy"] = conditioned["direct_unconditioned_energy"]
        output["class_reconstructed_energy"] = conditioned["class_reconstructed_energy"]
        output["energy_closure_error"] = conditioned["absolute_closure_error"]
    return output


__all__ = [
    "RunArtifacts",
    "environment_report",
    "experiment_config",
    "inspect_architectures",
    "quickstart_report",
    "run_and_save",
    "save_experiment",
    "save_experiment_result",
]
