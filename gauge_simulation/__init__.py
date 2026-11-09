"""Gauge-invariant imaginary-time simulation.

The historic TFIM/XY API remains available at the package root so the
published notebooks keep working.  The explicit TFIM extension is available
as :mod:`gauge_simulation.ancilla_resolved`.
"""

from .hamiltonian import GaugeModel
from . import models
from . import plotting
from . import ancilla_resolved

__all__ = ["GaugeModel", "ancilla_resolved", "models", "plotting"]
