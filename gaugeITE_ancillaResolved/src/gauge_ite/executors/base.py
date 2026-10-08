"""Common executor and asynchronous-job interfaces."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable
from typing import Any

from ..circuits import CircuitBundle
from ..execution import CircuitExecution


def backend_name(backend: Any) -> str:
    """Return a stable readable name for a simulator or hardware backend."""

    name = getattr(backend, "name", None)
    if callable(name):
        name = name()
    return str(name or backend.__class__.__name__)


def raw_job_id(job: Any) -> str | None:
    """Read a provider job identifier without assuming a specific SDK."""

    value = getattr(job, "job_id", None)
    if callable(value):
        value = value()
    return None if value is None else str(value)


class ExecutionJob:
    """Provider-neutral handle for a submitted circuit.

    Submission and result collection are separate, which matters for hardware
    queues.  Calling :meth:`result` waits for completion and converts provider
    output into the same :class:`CircuitExecution` used by Aer.
    """

    def __init__(
        self,
        *,
        provider: str,
        backend: str,
        raw_job: Any,
        loader: Callable[[Any], CircuitExecution],
        job_id: str | None = None,
    ) -> None:
        self.provider = provider
        self.backend_name = backend
        self.raw_job = raw_job
        self.job_id = job_id if job_id is not None else raw_job_id(raw_job)
        self._loader = loader
        self._execution: CircuitExecution | None = None

    def status(self) -> str:
        value = getattr(self.raw_job, "status", None)
        if callable(value):
            value = value()
        name = getattr(value, "name", value)
        return str(name or "UNKNOWN")

    def result(self) -> CircuitExecution:
        if self._execution is None:
            self._execution = self._loader(self.raw_job)
        return self._execution


class BackendExecutor(ABC):
    """Minimal interface implemented by all execution providers."""

    provider: str

    @abstractmethod
    def submit(self, bundle: CircuitBundle, basis: str = "ancilla") -> ExecutionJob:
        """Submit one measured circuit and return its job handle."""

    def run(self, bundle: CircuitBundle, basis: str = "ancilla") -> CircuitExecution:
        return self.submit(bundle, basis=basis).result()


__all__ = ["BackendExecutor", "ExecutionJob", "backend_name", "raw_job_id"]
