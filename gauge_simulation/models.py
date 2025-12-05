from .hamiltonian import Hamiltonian
from .core import tfim_generalized, xy_generalized

def build_hamiltonian(model, lattice_shape, trotter_steps=2):
    rows, cols = lattice_shape
    num_qubits = rows * cols
    dimension = 1 if (rows == 1 or cols == 1) else 2

    if model == "tfim":
        sim = tfim_generalized(dimension, num_qubits, trotter_steps)
    else:
        sim = xy_generalized(dimension, num_qubits, trotter_steps)

    # sim["hamiltonian"] is a SparsePauliOp
    return Hamiltonian(
        model=model,
        dimension=dimension,
        num_qubits=num_qubits,
        trotter_steps=trotter_steps,
        H_op=sim["hamiltonian"]
    )
