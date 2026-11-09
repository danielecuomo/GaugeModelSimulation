"""Reusable study objects composed from the public experiment API."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Sequence

import numpy as np
import pandas as pd

from .config import ExecutionConfig
from .experiments import ExperimentResult, GaugeITEExperiment
from .models import TFIMSpec
from .operational import operational_class_energy_sweep
from .protocols import ARCHITECTURES, ITEProtocol
from .thermodynamics import (
    class_prior_from_random_bonds,
    exact_class_thermodynamics,
    quadrature_convergence_report,
    quenched_free_energy_summary,
)


@dataclass(frozen=True)
class BetaSweepStudy:
    """Run one observable over a caller-selected beta grid."""

    spec: TFIMSpec
    protocol: ITEProtocol
    beta_values: Sequence[float]
    execution: ExecutionConfig = ExecutionConfig()
    observable: str = "class_probabilities"

    def run(self) -> list[ExperimentResult]:
        return [
            GaugeITEExperiment(
                self.spec,
                replace(self.protocol, beta=float(beta)),
                execution=self.execution,
            ).run(self.observable)
            for beta in self.beta_values
        ]


@dataclass(frozen=True)
class ArchitectureComparisonStudy:
    """Run the same physical input through selected circuit architectures."""

    spec: TFIMSpec
    protocol: ITEProtocol
    execution: ExecutionConfig = ExecutionConfig()
    architectures: Sequence[str] = tuple(sorted(ARCHITECTURES))
    observable: str = "class_probabilities"

    def run(self) -> dict[str, ExperimentResult]:
        return {
            architecture: GaugeITEExperiment(
                self.spec,
                replace(self.protocol, architecture=architecture),
                execution=self.execution,
            ).run(self.observable)
            for architecture in self.architectures
        }


@dataclass(frozen=True)
class QuenchedFreeEnergyStudy:
    """Reconstruct quenched free energy from class-resolved energy curves."""

    spec: TFIMSpec
    beta_values: Sequence[float]
    shots_per_basis: int = 4096
    integration_method: str = "simpson"

    def run(self) -> dict[str, object]:
        beta_values = np.asarray(self.beta_values, dtype=float)
        exact = exact_class_thermodynamics(self.spec, beta_values)
        operational = operational_class_energy_sweep(
            self.spec,
            beta_values,
            shots_per_basis=self.shots_per_basis,
        )
        prior = class_prior_from_random_bonds(self.spec, p_negative=0.5)
        result = quenched_free_energy_summary(
            exact,
            operational,
            prior,
            dimension=2**self.spec.num_system_qubits,
            method=self.integration_method,
        )
        result["exact_class_thermodynamics"] = exact
        result["operational_class_energies"] = operational
        return result

    def convergence(self, point_counts: Sequence[int] = (5, 9, 15, 29)) -> pd.DataFrame:
        beta_max = float(np.max(np.asarray(self.beta_values, dtype=float)))
        return quadrature_convergence_report(
            self.spec,
            beta_max,
            point_counts=point_counts,
            method=self.integration_method,
        )


__all__ = [
    "ArchitectureComparisonStudy",
    "BetaSweepStudy",
    "QuenchedFreeEnergyStudy",
]
