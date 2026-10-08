"""Local Aer implementation of the common executor interface."""

from __future__ import annotations

import warnings

from qiskit import transpile
from scipy.sparse import SparseEfficiencyWarning

from ..circuits import CircuitBundle, add_final_measurements
from ..execution import CircuitExecution
from .base import BackendExecutor, ExecutionJob, backend_name, raw_job_id


class AerExecutor(BackendExecutor):
    """Execute circuits locally, using matrix-product-state Aer by default."""

    provider = "aer"

    def __init__(
        self,
        *,
        shots: int = 4096,
        seed: int | None = 2026,
        backend=None,
        optimization_level: int = 0,
        method: str = "matrix_product_state",
    ) -> None:
        if shots <= 0:
            raise ValueError("shots must be positive.")
        if optimization_level not in range(4):
            raise ValueError("optimization_level must be 0, 1, 2, or 3.")
        if backend is None or isinstance(backend, str):
            try:
                from qiskit_aer import Aer, AerSimulator
            except ImportError as exc:  # pragma: no cover - environment dependent
                raise ImportError(
                    "Aer is a separate package. Install this repository with "
                    "`python -m pip install -e .` to install qiskit-aer."
                ) from exc
            if backend is None or backend.lower() == "auto":
                backend = AerSimulator(method=method)
            else:
                backend = Aer.get_backend(backend)
        self.backend = backend
        self.shots = int(shots)
        self.seed = seed
        self.optimization_level = int(optimization_level)

    def submit(self, bundle: CircuitBundle, basis: str = "ancilla") -> ExecutionJob:
        measured = add_final_measurements(bundle, basis=basis)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)
            compiled = transpile(
                measured,
                backend=self.backend,
                optimization_level=self.optimization_level,
            )
            raw_job = self.backend.run(
                compiled,
                shots=self.shots,
                seed_simulator=self.seed,
            )

        name = backend_name(self.backend)
        identifier = raw_job_id(raw_job)

        def load(job) -> CircuitExecution:
            result = job.result()
            counts = {
                str(key): int(value)
                for key, value in result.get_counts(0).items()
            }
            metadata = {}
            if getattr(result, "results", None):
                metadata = dict(getattr(result.results[0], "metadata", {}) or {})
            return CircuitExecution(
                bundle=bundle,
                basis=basis,
                shots=self.shots,
                seed=self.seed,
                measured_circuit=measured,
                transpiled_circuit=compiled,
                counts=counts,
                backend_name=name,
                provider=self.provider,
                job_id=identifier,
                register_counts={"measurement": counts},
                metadata=metadata,
            )

        return ExecutionJob(
            provider=self.provider,
            backend=name,
            raw_job=raw_job,
            loader=load,
            job_id=identifier,
        )


__all__ = ["AerExecutor"]
