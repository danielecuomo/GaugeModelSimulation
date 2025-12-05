from .hamiltonian import Hamiltonian
from .core import tfim_generalized, xy_generalized

def build_hamiltonian(model, size, trotter_steps=2):
    dimension = size[0]
    num_qubits = size[1]**size[0]
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
