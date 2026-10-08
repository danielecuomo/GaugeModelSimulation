"""Create provider executors from one small configuration object."""

from __future__ import annotations

from ..config import ExecutionConfig
from .aer import AerExecutor
from .base import BackendExecutor


def create_executor(config: ExecutionConfig | None = None, **overrides) -> BackendExecutor:
    """Return the requested executor; Aer is the default."""

    if config is None:
        config = ExecutionConfig(**overrides)
    elif overrides:
        raise ValueError("Pass either an ExecutionConfig or keyword overrides, not both.")

    if config.provider == "aer":
        backend = (
            None
            if isinstance(config.backend, str) and config.backend.lower() == "auto"
            else config.backend
        )
        return AerExecutor(
            shots=config.shots,
            seed=config.seed,
            backend=backend,
            optimization_level=config.resolved_optimization_level,
        )
    if config.provider == "ibm":
        from .ibm import IBMRuntimeExecutor

        return IBMRuntimeExecutor(
            shots=config.shots,
            backend=config.backend,
            profile=config.profile,
            instance=config.instance,
            optimization_level=config.resolved_optimization_level,
            submit=config.submit,
        )
    if config.provider == "quantinuum":
        if isinstance(config.backend, str):
            raise ValueError(
                "Quantinuum requires a configured backend object in Python."
            )
        from .quantinuum import QuantinuumExecutor

        return QuantinuumExecutor(
            config.backend,
            shots=config.shots,
            optimization_level=config.resolved_optimization_level,
        )
    raise AssertionError(f"Unhandled provider: {config.provider}")


__all__ = ["create_executor"]
