"""Optional Quantinuum implementation of the common executor interface."""

from __future__ import annotations

from ..circuits import CircuitBundle, add_final_measurements
from ..execution import CircuitExecution
from .base import BackendExecutor, ExecutionJob, backend_name


class QuantinuumExecutor(BackendExecutor):
    """Execute with a caller-configured pytket Quantinuum backend."""

    provider = "quantinuum"

    def __init__(self, backend, *, shots: int = 4096, optimization_level: int = 2):
        if shots <= 0:
            raise ValueError("shots must be positive.")
        self.backend = backend
        self.shots = int(shots)
        self.optimization_level = int(optimization_level)

    def submit(self, bundle: CircuitBundle, basis: str = "ancilla") -> ExecutionJob:
        try:
            from pytket.extensions.qiskit import qiskit_to_tk
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError(
                "Install `gauge-ite[quantinuum]` to use this executor."
            ) from exc
        measured = add_final_measurements(bundle, basis=basis)
        tk_circuit = qiskit_to_tk(measured)
        compiled = self.backend.get_compiled_circuit(
            tk_circuit,
            optimisation_level=self.optimization_level,
        )
        handle = self.backend.process_circuit(compiled, n_shots=self.shots)
        name = backend_name(self.backend)

        def load(raw_handle) -> CircuitExecution:
            result = self.backend.get_result(raw_handle)
            counts = {str(key): int(value) for key, value in result.get_counts().items()}
            return CircuitExecution(
                bundle=bundle,
                basis=basis,
                shots=sum(counts.values()),
                seed=None,
                measured_circuit=measured,
                transpiled_circuit=None,
                counts=counts,
                backend_name=name,
                provider=self.provider,
                job_id=str(raw_handle),
                metadata={"compiled_type": type(compiled).__name__},
            )

        return ExecutionJob(
            provider=self.provider,
            backend=name,
            raw_job=handle,
            loader=load,
            job_id=str(handle),
        )


__all__ = ["QuantinuumExecutor"]
