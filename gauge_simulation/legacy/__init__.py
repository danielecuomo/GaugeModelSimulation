"""Stable compatibility namespace for the original notebook API.

New notebook code can import from this namespace to make its use of the
historic interface explicit.  Existing imports from :mod:`gauge_simulation`
remain supported unchanged.
"""

from .. import models, plotting
from ..hamiltonian import GaugeModel
from ..models import build_hamiltonian

__all__ = ["GaugeModel", "build_hamiltonian", "models", "plotting"]
