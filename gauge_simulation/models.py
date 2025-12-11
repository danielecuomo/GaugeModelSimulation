# gauge_simulation/models.py

from __future__ import annotations

from .core import tfim_generalized, xy_generalized
from .hamiltonian import GaugeModel 

def build_hamiltonian(
    model: str,
    lattice_size: list[int],
    trotter_steps: int = 2,
) -> GaugeModel:
    """
    Initialize a GaugeModel object from the generalized TFIM / XY routines.

    Parameters
    ----------
    model : "tfim" or "xy"
    lattice_size  : A list containing [rows, cols].
        [rows, cols]; e.g. [1, 3] for a 3-site 1D chain, or [2, 2] for a 4-site 2D lattice.
        This must be a 2-element iterable defining the lattice dimensions.
    trotter_steps : int
        Number of Trotter steps used in the imaginary-time evolution (ITE) circuit.

    Returns
    -------
    GaugeModel
        Object containing methods for both exact diagonalization and 
        circuit-based simulation (ITE).
    
    Raises
    ------
    ValueError
        If the total number of qubits exceeds 20.
    """
    if not (isinstance(lattice_size, list) and len(lattice_size) == 2):
        raise ValueError("lattice_size must be a list of two integers [rows, cols].")
    
    rows, cols = lattice_size

    if not (isinstance(rows, int) and isinstance(cols, int) and rows > 0 and cols > 0):
        raise ValueError("Both rows and cols must be positive integers.")

    model_low = model.lower()
    if model_low == "tfim":
        # The core function checks for qubit limit and throws ValueError
        sim_data = tfim_generalized(lattice_size, trotter_steps)
    elif model_low == "xy":
        sim_data = xy_generalized(lattice_size, trotter_steps)
    else:
        raise ValueError("model must be 'tfim' or 'xy'.")

    return GaugeModel(
        model=model_low,
        lattice_size=lattice_size,
        trotter_steps=trotter_steps,
        sim_data=sim_data,
    )