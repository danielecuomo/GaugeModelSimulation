"""IBM Quantum Runtime executor with safe account and submission defaults."""

from __future__ import annotations

from collections import Counter
from typing import Any

from ..circuits import CircuitBundle, add_final_measurements
from ..execution import CircuitExecution
from .base import BackendExecutor, ExecutionJob, backend_name, raw_job_id


def supports_measure_reset(backend: Any) -> bool:
    """Return whether a backend exposes measurement and reset operations.

    A measure-reset trajectory contains mid-circuit measurements.  IBM
    backends evolve, so capability is checked on the selected backend rather
    than encoded as a list of device names.
    """

    target = getattr(backend, "target", None)
    operations = set(getattr(target, "operation_names", ()) or ())
    if not operations:
        configuration = getattr(backend, "configuration", None)
        configuration = configuration() if callable(configuration) else configuration
        operations = set(getattr(configuration, "basis_gates", ()) or ())
    return {"measure", "reset"}.issubset(operations)


def _counts(bit_array: Any) -> dict[str, int]:
    return {str(key): int(value) for key, value in bit_array.get_counts().items()}


def _joint_register_counts(data: Any) -> dict[str, int]:
    """Combine system and record registers shot by shot.

    Marginal count dictionaries cannot recover correlations.  Energy
    estimation needs the system outcome paired with its record outcome, so we
    use the Runtime BitArray shot strings and form the canonical Qiskit key
    ``"system record"``.
    """

    system = list(data.system.get_bitstrings())
    record = list(data.record.get_bitstrings())
    if len(system) != len(record):
        raise RuntimeError("IBM result registers contain different shot counts.")
    return dict(Counter(f"{system_bits} {record_bits}" for system_bits, record_bits in zip(system, record)))


class IBMRuntimeExecutor(BackendExecutor):
    """Submit circuits through IBM Runtime SamplerV2.

    ``submit=False`` is the default safety lock.  The service loads a locally
    saved account profile; this class never accepts or stores an API token.
    """

    provider = "ibm"

    def __init__(
        self,
        *,
        shots: int = 4096,
        backend: str | Any = "auto",
        profile: str | None = None,
        instance: str | None = None,
        optimization_level: int = 3,
        submit: bool = False,
        service: Any = None,
    ) -> None:
        if shots <= 0:
            raise ValueError("shots must be positive.")
        if optimization_level not in range(4):
            raise ValueError("optimization_level must be 0, 1, 2, or 3.")
        self.shots = int(shots)
        self.backend_selection = backend
        self.profile = profile
        self.instance = instance
        self.optimization_level = int(optimization_level)
        self.hardware_submission_enabled = bool(submit)
        self._service = service

    def _runtime_classes(self):
        try:
            from qiskit.transpiler import generate_preset_pass_manager
            from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise ImportError(
                "Install `gauge-ite[ibm]` to use IBM Quantum hardware."
            ) from exc
        return QiskitRuntimeService, SamplerV2, generate_preset_pass_manager

    @property
    def service(self):
        if self._service is None:
            service_class, _, _ = self._runtime_classes()
            kwargs = {}
            if self.profile:
                kwargs["name"] = self.profile
            if self.instance:
                kwargs["instance"] = self.instance
            self._service = service_class(**kwargs)
        return self._service

    def select_backend(self, bundle: CircuitBundle):
        selection = self.backend_selection
        if not isinstance(selection, str):
            backend = selection
        elif selection.lower() == "auto":
            filters = None
            if bundle.protocol.architecture == "measure_reset_trajectory":
                filters = supports_measure_reset
            backend = self.service.least_busy(
                min_num_qubits=bundle.circuit.num_qubits,
                operational=True,
                simulator=False,
                filters=filters,
            )
        else:
            kwargs = {"instance": self.instance} if self.instance else {}
            backend = self.service.backend(selection, **kwargs)

        if getattr(backend, "num_qubits", 0) < bundle.circuit.num_qubits:
            raise ValueError(
                f"Backend {backend_name(backend)!r} has {backend.num_qubits} qubits; "
                f"the circuit requires {bundle.circuit.num_qubits}."
            )
        status_method = getattr(backend, "status", None)
        status = status_method() if callable(status_method) else None
        if status is not None and not getattr(status, "operational", True):
            raise ValueError(f"Backend {backend_name(backend)!r} is not operational.")
        if (
            bundle.protocol.architecture == "measure_reset_trajectory"
            and not supports_measure_reset(backend)
        ):
            raise ValueError(
                f"Backend {backend_name(backend)!r} does not advertise the "
                "measure/reset support required by measure_reset_trajectory."
            )
        return backend

    def submit(self, bundle: CircuitBundle, basis: str = "ancilla") -> ExecutionJob:
        if not self.hardware_submission_enabled:
            raise PermissionError(
                "IBM hardware submission is locked. Set submit=True in Python "
                "or add --submit on the command line after reviewing the job."
            )
        _, sampler_class, pass_manager_factory = self._runtime_classes()
        backend = self.select_backend(bundle)
        measured = add_final_measurements(bundle, basis=basis)
        pass_manager = pass_manager_factory(
            backend=backend,
            optimization_level=self.optimization_level,
        )
        compiled = pass_manager.run(measured)
        sampler = sampler_class(mode=backend)
        raw_job = sampler.run([compiled], shots=self.shots)
        name = backend_name(backend)
        identifier = raw_job_id(raw_job)

        def load(job) -> CircuitExecution:
            pub_result = job.result()[0]
            data = pub_result.data
            register_counts: dict[str, dict[str, int]] = {}
            for register in compiled.cregs:
                bit_array = getattr(data, register.name, None)
                if bit_array is not None:
                    register_counts[register.name] = _counts(bit_array)

            if hasattr(data, "measurement"):
                counts = _counts(data.measurement)
            elif basis in {"z", "x"} and hasattr(data, "system") and hasattr(data, "record"):
                counts = _joint_register_counts(data)
            elif hasattr(data, "record"):
                counts = _counts(data.record)
            elif len(register_counts) == 1:
                counts = next(iter(register_counts.values()))
            else:
                raise RuntimeError(
                    "Could not identify a count register in the IBM Sampler result."
                )

            return CircuitExecution(
                bundle=bundle,
                basis=basis,
                shots=sum(counts.values()),
                seed=None,
                measured_circuit=measured,
                transpiled_circuit=compiled,
                counts=counts,
                backend_name=name,
                provider=self.provider,
                job_id=identifier,
                register_counts=register_counts,
                metadata=dict(getattr(pub_result, "metadata", {}) or {}),
            )

        return ExecutionJob(
            provider=self.provider,
            backend=name,
            raw_job=raw_job,
            loader=load,
            job_id=identifier,
        )


__all__ = ["IBMRuntimeExecutor", "supports_measure_reset"]
