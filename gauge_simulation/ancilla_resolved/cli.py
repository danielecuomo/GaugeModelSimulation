"""Command-line entry point for common Gauge ITE tasks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import numpy as np

from .config import ExecutionConfig, ExperimentConfig, OBSERVABLES
from .circuits import build_tfim_circuit
from .execution import decode_counts
from .experiments import ExperimentResult, GaugeITEExperiment
from .models import TFIMSpec
from .protocols import ANGLE_MAPPINGS, ARCHITECTURES, INPUT_MODES, ITEProtocol
from .workflows import (
    inspect_architectures,
    quickstart_report,
)


def _comma_values(value: str, cast=float) -> tuple:
    try:
        return tuple(cast(part.strip()) for part in value.split(",") if part.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def _sign_values(value: str) -> tuple[int, ...]:
    signs = _comma_values(value, int)
    if any(sign not in (-1, 1) for sign in signs):
        raise argparse.ArgumentTypeError("signs must contain only +1 and -1")
    return signs


def _bit_values(value: str) -> tuple[int, ...]:
    compact = value.replace(",", "").replace(" ", "")
    if not compact or set(compact) - {"0", "1"}:
        raise argparse.ArgumentTypeError("basis bits must contain only 0 and 1")
    return tuple(int(bit) for bit in compact)


def _add_model_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--rows", type=int, default=2)
    parser.add_argument("--cols", type=int, default=2)
    parser.add_argument("--boundary", choices=("open", "periodic"), default="open")
    parser.add_argument("--gamma", type=float, default=float(np.pi / 4))
    parser.add_argument("--eta", type=float, default=float(np.pi / 4))
    parser.add_argument(
        "--compiled-bond-signs",
        type=_sign_values,
        help="Comma-separated fixed bond signs, for example 1,-1,1,1.",
    )


def _add_protocol_arguments(
    parser: argparse.ArgumentParser,
    *,
    include_architecture: bool = True,
) -> None:
    parser.add_argument("--beta", type=float, default=0.75)
    parser.add_argument("--steps", type=int, default=1)
    if include_architecture:
        parser.add_argument(
            "--architecture",
            choices=tuple(sorted(ARCHITECTURES)),
            default="coherent_reuse",
        )
    parser.add_argument(
        "--angle-mapping",
        choices=tuple(sorted(ANGLE_MAPPINGS)),
        default="exact_kraus",
    )
    parser.add_argument(
        "--input-mode",
        choices=tuple(sorted(INPUT_MODES)),
        default="purification",
    )
    parser.add_argument("--basis-bits", type=_bit_values)


def _spec_from_args(args) -> TFIMSpec:
    return TFIMSpec.rectangular(
        rows=args.rows,
        cols=args.cols,
        gamma=args.gamma,
        eta=args.eta,
        compiled_bond_signs=args.compiled_bond_signs,
        boundary=args.boundary,
    )


def _protocol_from_args(args) -> ITEProtocol:
    return ITEProtocol(
        beta=args.beta,
        steps=args.steps,
        architecture=args.architecture,
        angle_mapping=args.angle_mapping,
        input_mode=args.input_mode,
        basis_bits=args.basis_bits,
    )


def _command_quickstart(args) -> int:
    spec = _spec_from_args(args)
    protocol = _protocol_from_args(args)
    beta_values = args.beta_values or (0.0, 0.25, 0.5, args.beta)
    report = quickstart_report(spec, protocol, beta_values)
    print(f"Exact clean ground-state energy: {report['ground_energy']:.8f}")
    if "circuit_energy" in report:
        print(f"Circuit energy at beta={args.beta:g}: {report['circuit_energy']:.8f}")
        print(f"Class reconstruction error: {report['energy_closure_error']:.3e}")
    print("\nExact clean thermodynamics")
    print(report["clean_thermodynamics"].round(8).to_string(index=False))
    print("\nCircuit gauge-class probabilities")
    print(report["class_probabilities"].round(8).to_string(index=False))
    if args.output:
        output = Path(args.output)
        output.mkdir(parents=True, exist_ok=True)
        report["clean_thermodynamics"].to_csv(output / "clean_thermodynamics.csv", index=False)
        report["class_probabilities"].to_csv(output / "class_probabilities.csv", index=False)
        report["resources"].to_csv(output / "resources.csv", index=False)
        print(f"\nSaved tables to {output.resolve()}")
    return 0


def _command_run(args) -> int:
    spec = _spec_from_args(args)
    protocol = _protocol_from_args(args)
    execution_config = ExecutionConfig(
        provider=args.provider,
        backend=args.backend,
        shots=args.shots,
        seed=args.seed,
        optimization_level=args.optimization_level,
        profile=args.profile,
        instance=args.instance,
        submit=args.submit,
    )
    experiment = GaugeITEExperiment(spec, protocol, execution=execution_config)
    if args.observable == "energy":
        result = experiment.run("energy")
    elif args.basis == "ancilla":
        result = experiment.run(args.observable)
    else:
        raw_execution = experiment.execute(args.basis)
        decoded = decode_counts(experiment.bundle, raw_execution.counts)
        result = ExperimentResult(
            config=ExperimentConfig(
                spec=spec,
                protocol=protocol,
                observable=args.observable,
                execution=execution_config,
            ),
            bundle=experiment.bundle,
            executions={args.basis: raw_execution},
            analysis={**decoded, "interpretation": protocol.interpretation},
            warnings=experiment.interpretation_notes(args.observable),
        )

    if args.save:
        artifacts = result.save(
            output_root=args.output_root,
            run_name=args.run_name,
            save_mpl=not args.no_png,
        )
        output_dir = artifacts.output_dir if hasattr(artifacts, "output_dir") else artifacts["output_dir"]
        print(f"Saved complete run to {output_dir.resolve()}")
    else:
        print("Run completed in memory; no output files were written.")
        print("Add --save if you want circuits, counts, and tables on disk.")

    for basis, execution in result.executions.items():
        print(
            f"{basis.upper()} run: provider={execution.provider}, "
            f"backend={execution.backend_name}, job_id={execution.job_id or 'local'}"
        )
    for note in result.warnings:
        print(f"Interpretation: {note}")
    if args.observable == "energy":
        summary = {
            key: result.analysis[key]
            for key in (
                "direct_unconditioned_energy",
                "class_reconstructed_energy",
                "absolute_closure_error",
            )
        }
        classes = result.analysis["class_energies"]
        print(json.dumps(summary, indent=2))
        print("\nClass-conditioned energy table")
    else:
        summary = result.analysis["summary"]
        classes = result.analysis["classes"]
        print(json.dumps(summary, indent=2))
        print("\nGauge-class table")
    print(classes.round(8).to_string(index=False))
    return 0


def _command_configure_ibm(args) -> int:
    from .accounts import configure_ibm_account

    result = configure_ibm_account(
        profile=args.profile,
        instance=args.instance,
        channel=args.channel,
        overwrite=args.overwrite,
        set_as_default=not args.no_default,
    )
    print(json.dumps(result, indent=2))
    print("The token was saved by Qiskit Runtime and was not written to this repository.")
    return 0


def _command_account_check(args) -> int:
    from .accounts import ibm_account_status

    result = ibm_account_status(profile=args.profile)
    print(json.dumps(result, indent=2))
    return 0 if result["configured"] else 1


def _command_draw(args) -> int:
    table = inspect_architectures(
        _spec_from_args(args),
        beta=args.beta,
        steps=args.steps,
        angle_mapping=args.angle_mapping,
        input_mode=args.input_mode,
        basis_bits=args.basis_bits,
        output_dir=args.output,
        optimization_level=args.optimization_level,
        save_mpl=not args.no_png,
    )
    print(table.to_string(index=False))
    print(f"\nSaved circuit files to {Path(args.output).resolve()}")
    print(
        "Persistent records are selected from fresh or measure-reset data; "
        "they are not a separate circuit architecture."
    )
    return 0


def _command_classify(args) -> int:
    counts = json.loads(Path(args.counts).read_text(encoding="utf-8"))
    bundle = build_tfim_circuit(_spec_from_args(args), _protocol_from_args(args))
    result = decode_counts(bundle, counts)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    result["branches"].to_csv(output / "branches.csv", index=False)
    result["classes"].to_csv(output / "classes.csv", index=False)
    (output / "summary.json").write_text(
        json.dumps(result["summary"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(result["classes"].round(8).to_string(index=False))
    print(f"\nSaved decoded tables to {output.resolve()}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gauge-ite",
        description="Ancilla-resolved imaginary-time and gauge-class simulation.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    quickstart = subparsers.add_parser(
        "quickstart",
        help="Print ground-state, thermal, circuit-energy, and class information.",
    )
    _add_model_arguments(quickstart)
    _add_protocol_arguments(quickstart)
    quickstart.add_argument("--beta-values", type=lambda value: _comma_values(value, float))
    quickstart.add_argument("--output")
    quickstart.set_defaults(func=_command_quickstart)

    run = subparsers.add_parser(
        "run",
        help="Run on Aer (default) or explicitly submit to IBM hardware.",
    )
    _add_model_arguments(run)
    _add_protocol_arguments(run)
    run.add_argument("--basis", choices=("ancilla", "z", "x"), default="ancilla")
    run.add_argument(
        "--observable",
        choices=tuple(sorted(OBSERVABLES)),
        default="class_probabilities",
        help="Energy automatically runs both Z and X measurement circuits.",
    )
    run.add_argument("--provider", choices=("aer", "ibm"), default="aer")
    run.add_argument(
        "--backend",
        default="auto",
        help="Aer uses its default simulator; IBM 'auto' chooses the least busy compatible device.",
    )
    run.add_argument("--profile", help="Name of a locally saved IBM Runtime profile.")
    run.add_argument("--instance", help="Optional IBM Cloud CRN or service instance.")
    run.add_argument(
        "--submit",
        action="store_true",
        help="Required safety confirmation before a real IBM hardware job is submitted.",
    )
    run.add_argument("--shots", type=int, default=4096)
    run.add_argument("--seed", type=int, default=2026)
    run.add_argument(
        "--optimization-level",
        type=int,
        choices=range(4),
        help="Defaults to 0 for Aer and 3 for IBM hardware.",
    )
    run.add_argument(
        "--save",
        action="store_true",
        help="Write circuits, counts, tables, and metadata to a run directory.",
    )
    run.add_argument("--output-root", default="runs")
    run.add_argument("--run-name")
    run.add_argument("--no-png", action="store_true")
    run.set_defaults(func=_command_run)

    draw = subparsers.add_parser(
        "draw",
        help="Compare coherent, fresh, and measure-reset circuit architectures.",
    )
    _add_model_arguments(draw)
    _add_protocol_arguments(draw, include_architecture=False)
    draw.add_argument("--optimization-level", type=int, choices=range(4), default=0)
    draw.add_argument("--output", default="artifacts/circuits")
    draw.add_argument("--no-png", action="store_true")
    draw.set_defaults(func=_command_draw)

    classify = subparsers.add_parser(
        "classify",
        help="Decode an existing Qiskit counts JSON file into gauge classes.",
    )
    _add_model_arguments(classify)
    _add_protocol_arguments(classify)
    classify.add_argument("counts")
    classify.add_argument("--output", default="decoded_counts")
    classify.set_defaults(func=_command_classify)

    configure = subparsers.add_parser(
        "configure",
        help="Configure a provider account without placing credentials in code.",
    )
    configure_subparsers = configure.add_subparsers(dest="provider_name", required=True)
    configure_ibm = configure_subparsers.add_parser("ibm")
    configure_ibm.add_argument("--profile", default="default")
    configure_ibm.add_argument("--instance")
    configure_ibm.add_argument("--channel", default="ibm_quantum_platform")
    configure_ibm.add_argument("--overwrite", action="store_true")
    configure_ibm.add_argument("--no-default", action="store_true")
    configure_ibm.set_defaults(func=_command_configure_ibm)

    account = subparsers.add_parser("account", help="Inspect saved provider profiles.")
    account_subparsers = account.add_subparsers(dest="account_command", required=True)
    account_check = account_subparsers.add_parser("check")
    account_check.add_argument("--profile")
    account_check.set_defaults(func=_command_account_check)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line interface."""

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (
        ValueError,
        FileNotFoundError,
        FileExistsError,
        ImportError,
        PermissionError,
        RuntimeError,
    ) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
