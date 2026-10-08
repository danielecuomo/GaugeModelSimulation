"""Small object-oriented facade over the unchanged physics and circuit code."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

from .circuits import CircuitBundle, build_tfim_circuit
from .config import ExecutionConfig, ExperimentConfig, OBSERVABLES
from .execution import CircuitExecution, decode_counts
from .executors import BackendExecutor, ExecutionJob, create_executor
from .measurements import conditioned_energy_from_counts
from .models import TFIMSpec
from .protocols import ITEProtocol


@dataclass(frozen=True)
class ExperimentResult:
    """Provider-neutral results plus their physical interpretation."""

    config: ExperimentConfig
    bundle: CircuitBundle
    executions: dict[str, CircuitExecution]
    analysis: dict[str, Any]
    warnings: tuple[str, ...] = field(default_factory=tuple)

    @property
    def interpretation(self) -> str:
        return self.bundle.interpretation

    @property
    def job_ids(self) -> dict[str, str | None]:
        return {basis: execution.job_id for basis, execution in self.executions.items()}

    @property
    def backend_names(self) -> dict[str, str]:
        return {
            basis: execution.backend_name
            for basis, execution in self.executions.items()
        }

    @property
    def classes(self):
        return self.analysis.get("class_energies", self.analysis.get("classes"))

    def save(
        self,
        output_root: str | Path = "runs",
        run_name: str | None = None,
        *,
        save_mpl: bool = True,
    ):
        """Save only when explicitly called."""

        from .workflows import save_experiment_result

        return save_experiment_result(
            self,
            output_root=output_root,
            run_name=run_name,
            save_mpl=save_mpl,
        )


class GaugeITEExperiment:
    """Configure once, then run on Aer or a supported hardware executor."""

    def __init__(
        self,
        spec: TFIMSpec,
        protocol: ITEProtocol,
        *,
        executor: BackendExecutor | None = None,
        execution: ExecutionConfig | None = None,
    ) -> None:
        if executor is not None and execution is not None:
            raise ValueError("Pass either executor or execution, not both.")
        self.spec = spec
        self.protocol = protocol
        self.execution_config = execution or ExecutionConfig()
        self.executor = executor or create_executor(self.execution_config)
        self._bundle: CircuitBundle | None = None

    @property
    def bundle(self) -> CircuitBundle:
        if self._bundle is None:
            self._bundle = build_tfim_circuit(self.spec, self.protocol)
        return self._bundle

    def interpretation_notes(self, observable: str = "class_probabilities") -> tuple[str, ...]:
        notes = [f"Circuit interpretation: {self.protocol.interpretation}."]
        if self.protocol.architecture == "coherent_reuse" and self.protocol.steps > 1:
            notes.append(
                "The final coherent record is not a step-by-step disorder history; "
                "do not interpret it as one fixed Hamiltonian reused at every step."
            )
        if observable == "energy" and self.protocol.steps > 1:
            notes.append(
                "The reported multi-step energy is conditioned on the final "
                "operational record, not a static-Hamiltonian Gibbs class."
            )
        return tuple(notes)

    def submit(self, basis: str = "ancilla") -> ExecutionJob:
        """Submit a raw basis measurement without waiting for its result."""

        return self.executor.submit(self.bundle, basis=basis)

    def execute(self, basis: str = "ancilla") -> CircuitExecution:
        """Submit and wait for one raw basis measurement."""

        return self.executor.run(self.bundle, basis=basis)

    def run(self, observable: str = "class_probabilities") -> ExperimentResult:
        """Run a high-level observable and return a common result object."""

        observable = observable.lower()
        if observable not in OBSERVABLES:
            raise ValueError(f"observable must be one of {sorted(OBSERVABLES)}.")
        config = ExperimentConfig(
            spec=self.spec,
            protocol=self.protocol,
            observable=observable,
            execution=self.execution_config,
        )
        if observable in {"class_probabilities", "trajectory_records"}:
            execution = self.execute("ancilla")
            decoded = decode_counts(self.bundle, execution.counts)
            analysis = {
                **decoded,
                "interpretation": self.protocol.interpretation,
            }
            executions = {"ancilla": execution}
        else:
            z_execution = self.execute("z")
            x_execution = self.execute("x")
            analysis = conditioned_energy_from_counts(
                self.bundle,
                z_execution.counts,
                x_execution.counts,
            )
            executions = {"z": z_execution, "x": x_execution}
        return ExperimentResult(
            config=config,
            bundle=self.bundle,
            executions=executions,
            analysis=analysis,
            warnings=self.interpretation_notes(observable),
        )

    def with_protocol(self, **changes) -> "GaugeITEExperiment":
        """Create a related experiment without mutating this one."""

        return GaugeITEExperiment(
            self.spec,
            replace(self.protocol, **changes),
            executor=self.executor,
        )


__all__ = ["ExperimentResult", "GaugeITEExperiment"]
