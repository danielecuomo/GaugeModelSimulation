"""Imaginary-time and gauge-class simulation tools.

The package keeps hardware-provider imports optional. Importing
``gauge_ite`` therefore exposes the mathematical and circuit-building API
without opening a provider session or configuring Matplotlib.
"""

from .branches import BranchRecord, record_from_bits, record_from_signs
from .config import ExecutionConfig, ExperimentConfig
from .execution import CircuitExecution, decode_counts, run_aer_experiment
from .experiments import ExperimentResult, GaugeITEExperiment
from .executors import AerExecutor, BackendExecutor, ExecutionJob, create_executor
from .measurements import conditioned_energy_from_counts
from .models import TFIMSpec, Term
from .operational import operational_class_energy_sweep
from .protocols import ITEProtocol
from .simulation import clean_ground_energy, clean_spectrum, clean_thermodynamics
from .studies import (
    ArchitectureComparisonStudy,
    BetaSweepStudy,
    QuenchedFreeEnergyStudy,
)
from .thermodynamics import (
    quadrature_convergence_report,
    quenched_free_energy_summary,
    reconstruct_class_log_partitions,
)
from .topology import RectangularLattice

__all__ = [
    "BranchRecord",
    "ArchitectureComparisonStudy",
    "BetaSweepStudy",
    "CircuitExecution",
    "ExecutionConfig",
    "ExecutionJob",
    "ExperimentConfig",
    "ExperimentResult",
    "GaugeITEExperiment",
    "ITEProtocol",
    "QuenchedFreeEnergyStudy",
    "RectangularLattice",
    "TFIMSpec",
    "Term",
    "AerExecutor",
    "BackendExecutor",
    "clean_ground_energy",
    "clean_spectrum",
    "clean_thermodynamics",
    "decode_counts",
    "conditioned_energy_from_counts",
    "create_executor",
    "operational_class_energy_sweep",
    "quadrature_convergence_report",
    "quenched_free_energy_summary",
    "record_from_bits",
    "record_from_signs",
    "reconstruct_class_log_partitions",
    "run_aer_experiment",
]

__version__ = "0.3.0"
