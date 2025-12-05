# gauge_simulation/initialization.py

from __future__ import annotations

from .core import tfim_generalized, xy_generalized
from .hamiltonian import Hamiltonian


def init(
    model: str = "TFIM",
    size: list[int] | tuple[int, int] = (1, 3),
    trotter_steps: int = 1,
) -> Hamiltonian:
    """
    Initialize a Hamiltonian object from the generalized TFIM / XY routines.

    Parameters
    ----------
    model : "TFIM" or "XY"
    size  : [rows, cols]; e.g. [1, 3] for a 3-site 1D chain.
    trotter_steps : int
        Number of Trotter steps used in the circuit construction (does not
        affect the exact Hamiltonian spectrum, only the underlying circuit size).

    Returns
    -------
    Hamiltonian
        Object with .ground_state() and .thermal_average(beta).
    """
    if size is None or len(size) != 2:
        raise ValueError("size must be a 2-element iterable [rows, cols].")

    rows, cols = int(size[0]), int(size[1])
    num_qubits = rows * cols

    if rows == 1 or cols == 1:
        dimension = 1
    else:
        dimension = 2

    model_low = model.lower()
    if model_low == "tfim":
        sim = tfim_generalized(dimension, num_qubits, trotter_steps)
    elif model_low == "xy":
        sim = xy_generalized(dimension, num_qubits, trotter_steps)
    else:
        raise ValueError("model must be 'TFIM' or 'XY'.")

    return Hamiltonian(
        model=model_low,
        num_qubits=num_qubits,
        hamiltonian_terms_0=sim["hamiltonian_terms_0"],
        hamiltonian_terms_1=sim["hamiltonian_terms_1"],
    )
