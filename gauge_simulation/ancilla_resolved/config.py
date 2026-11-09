"""Immutable configuration objects for reproducible experiments."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .models import TFIMSpec
from .protocols import ITEProtocol


PROVIDERS = {"aer", "ibm", "quantinuum"}
OBSERVABLES = {"class_probabilities", "trajectory_records", "energy"}


@dataclass(frozen=True)
class ExecutionConfig:
    """Describe where and how circuits are executed.

    Credentials are intentionally absent.  IBM credentials are loaded from a
    saved local profile, never embedded in a script, notebook, result, or
    object representation.
    """

    provider: str = "aer"
    backend: str | Any = "auto"
    shots: int = 4096
    seed: int | None = 2026
    optimization_level: int | None = None
    profile: str | None = None
    instance: str | None = None
    submit: bool = False

    def __post_init__(self) -> None:
        provider = str(self.provider).lower()
        if provider not in PROVIDERS:
            raise ValueError(f"provider must be one of {sorted(PROVIDERS)}.")
        if self.shots <= 0:
            raise ValueError("shots must be positive.")
        if self.optimization_level is not None and self.optimization_level not in range(4):
            raise ValueError("optimization_level must be 0, 1, 2, or 3.")
        object.__setattr__(self, "provider", provider)

    @property
    def resolved_optimization_level(self) -> int:
        if self.optimization_level is not None:
            return self.optimization_level
        return 3 if self.provider == "ibm" else 0


@dataclass(frozen=True)
class ExperimentConfig:
    """Complete scientific and execution configuration for one experiment."""

    spec: TFIMSpec
    protocol: ITEProtocol
    observable: str = "class_probabilities"
    execution: ExecutionConfig = ExecutionConfig()

    def __post_init__(self) -> None:
        observable = str(self.observable).lower()
        if observable not in OBSERVABLES:
            raise ValueError(f"observable must be one of {sorted(OBSERVABLES)}.")
        object.__setattr__(self, "observable", observable)


__all__ = [
    "ExecutionConfig",
    "ExperimentConfig",
    "OBSERVABLES",
    "PROVIDERS",
]
