# gauge_simulation/core.py

import numpy as np
from numpy import pi

from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.circuit import Parameter


def get_nearest_neighbor_interactions(dimension: int, num_qubits: int) -> list[tuple[int, int]]:
    """
    Generate nearest-neighbor interaction pairs for any 1D or 2D spin model.

    Args:
        dimension (int): 1 or 2 (dimension of the model).
        num_qubits (int): Total number of qubits. For 2D, must be a perfect square.

    Returns:
        list[tuple[int,int]]: list of (i,j) index pairs for nearest-neighbor couplings
    """
    if dimension == 1:
        # Linear 1D chain, extreme ends not connected
        interactions = [(i, i + 1) for i in range(num_qubits - 1)]

    elif dimension == 2:
        # Square 2D lattice
        L = int(np.sqrt(num_qubits))
        if L * L != num_qubits:
            raise ValueError("Number of qubits must form a perfect square lattice.")
        interactions = []
        for i in range(num_qubits):
            # Right neighbor
            if (i + 1) % L != 0:
                interactions.append((i, i + 1))
            # Down neighbor
            if i + L < num_qubits:
                interactions.append((i, i + L))

    else:
        raise ValueError("Only 1D and 2D xy models are supported.")

    return interactions


def kron_matrices(pauli_string: str, pauli_dict: dict) -> np.ndarray:
    """
    Compute the Kronecker product for a given Pauli string.

    Args:
        pauli_string (str): String representing the Pauli operators (e.g., 'ZZX').
        pauli_dict (dict): Dictionary mapping Pauli labels to their matrix representations.

    Returns:
        np.ndarray: The resulting matrix from the Kronecker product.
    """
    result = pauli_dict[pauli_string[0]]
    for char in pauli_string[1:]:
        result = np.kron(result, pauli_dict[char])
    return result


def thermal_average_H(beta: float, lambda_n: np.ndarray) -> float:
    """
    Calculate the thermal average of the Hamiltonian at a given beta.

    Args:
        beta (float): Inverse temperature.
        lambda_n (np.ndarray): Eigenvalues of the Hamiltonian.

    Returns:
        float: Thermal average.
    """
    w = np.exp(-beta * lambda_n)
    Z_f = np.sum(w)
    return float(np.sum(w * lambda_n) / Z_f)


def exp_func(beta: np.ndarray | float, a: float, b: float, c: float) -> np.ndarray | float:
    """
    Exponential decay function for curve fitting.

    Args:
        beta (np.ndarray | float): Input value(s).
        a (float): Amplitude.
        b (float): Decay rate.
        c (float): Offset.

    Returns:
        np.ndarray | float: Function value(s).
    """
    return a * np.exp(-beta * b) + c


def tfim_generalized(dimension: int, num_qubits: int, trotter_steps: int) -> dict:
    """
    Cuomo's gauge-invariant ITE for a transverse field Ising model in 1D or 2D.:contentReference[oaicite:1]{index=1}
    """
    zz_interaction_list = get_nearest_neighbor_interactions(dimension, num_qubits)

    if dimension == 2:
        L = int(np.sqrt(num_qubits))
        if L * L != num_qubits:
            raise ValueError("Number of qubits must form a perfect square lattice for 2D.")
        total_num_qubits = (5 * num_qubits) - 2 * L
    else:  # dimension == 1
        total_num_qubits = (4 * num_qubits) - 1

    if total_num_qubits > 20:
        raise ValueError(
            f"Total qubits required ({total_num_qubits}) exceeds 20. "
            "Simulation may require excessive memory. "
            "Consider reducing num_qubits or trotter_steps."
        )

    γ = Parameter("γ")
    η = Parameter("η")
    β = Parameter("β")
    X = SparsePauliOp("X")
    Z = SparsePauliOp("Z")
    ZZZ = PauliEvolutionGate(Z ^ Z ^ Z, time=(β / (2 * trotter_steps)) * γ)
    ZX = PauliEvolutionGate(Z ^ X, time=(β / (2 * trotter_steps)) * η)

    tfim_qc_sym = QuantumCircuit(total_num_qubits)

    # maximally mixed state for system qubits
    for ind in range(num_qubits):
        tfim_qc_sym.h(ind)
        tfim_qc_sym.cx(ind, total_num_qubits - 1 - ind)
    tfim_qc_sym.barrier()

    # auxiliary qubits in |+>
    aux_start = num_qubits
    aux_end = total_num_qubits - num_qubits
    for ind in range(aux_start, aux_end):
        tfim_qc_sym.h(ind)
    tfim_qc_sym.barrier()

    # Trotter steps
    for _ in range(trotter_steps):
        # ZZZ
        for ind, zz_interaction in enumerate(zz_interaction_list):
            tfim_qc_sym.append(ZZZ, [zz_interaction[0], zz_interaction[1], num_qubits + ind])
        tfim_qc_sym.barrier()

        # ZX
        for ind in range(num_qubits):
            tfim_qc_sym.append(ZX, [ind, num_qubits + len(zz_interaction_list) + ind])
        tfim_qc_sym.barrier()

    # rotate auxiliary qubits
    for ind in range(aux_start, aux_end):
        tfim_qc_sym.sx(ind)

    # build circuit Hamiltonian Pauli strings
    ZZZ_list = []
    ZX_list = []
    for ind, zz_interaction in enumerate(zz_interaction_list):
        h_term = list("I" * total_num_qubits)
        h_term[zz_interaction[0]] = "Z"
        h_term[zz_interaction[1]] = "Z"
        h_term[num_qubits + ind] = "Z"
        ZZZ_list.append("".join(h_term)[::-1])

    for ind in range(num_qubits):
        h_term = list("I" * total_num_qubits)
        h_term[ind] = "X"
        h_term[num_qubits + len(zz_interaction_list) + ind] = "Z"
        ZX_list.append("".join(h_term)[::-1])

    hamiltonian_terms = ZZZ_list + ZX_list
    hamiltonian_terms_coeff = [γ] * len(ZZZ_list) + [η] * len(ZX_list)

    tfim_hamiltonian = SparsePauliOp(hamiltonian_terms, coeffs=hamiltonian_terms_coeff)
    tfim_hamiltonian = tfim_hamiltonian.assign_parameters({γ: pi / 4, η: pi / 4})
    tfim_hamiltonian = SparsePauliOp(
        tfim_hamiltonian.paulis,
        coeffs=np.asarray(tfim_hamiltonian.coeffs, dtype=np.complex128),
    )

    return {
        "qc": tfim_qc_sym,
        "hamiltonian": tfim_hamiltonian,
        "hamiltonian_terms_0": ZZZ_list,
        "hamiltonian_terms_1": ZX_list,
        "total_num_qubits": total_num_qubits,
        "gamma_param": γ,
        "eta_param": η,
        "beta_param": β,
    }


def xy_generalized(dimension: int, num_qubits: int, trotter_steps: int) -> dict:
    """
    Cuomo's gauge-invariant ITE for an XY model in 1D or 2D.:contentReference[oaicite:2]{index=2}
    """
    interaction_list = get_nearest_neighbor_interactions(dimension, num_qubits)

    if dimension == 2:
        L = int(np.sqrt(num_qubits))
        if L * L != num_qubits:
            raise ValueError("Number of qubits must form a perfect square lattice for 2D.")
        total_num_qubits = (6 * num_qubits) - 4 * L
    else:  # dimension == 1
        total_num_qubits = (4 * num_qubits) - 2

    if total_num_qubits > 20:
        raise ValueError(
            f"Total qubits required ({total_num_qubits}) exceeds 20. "
            "Simulation may require excessive memory. "
            "Consider reducing num_qubits or trotter_steps."
        )

    γ = Parameter("γ")
    β = Parameter("β")
    X = SparsePauliOp("X")
    Z = SparsePauliOp("Z")
    ZZZ = PauliEvolutionGate(Z ^ Z ^ Z, time=(β / (2 * trotter_steps)) * γ)
    ZXX = PauliEvolutionGate(Z ^ X ^ X, time=(β / (2 * trotter_steps)) * γ)

    xy_qc_sym = QuantumCircuit(total_num_qubits)

    # maximally mixed state for system qubits
    for ind in range(num_qubits):
        xy_qc_sym.h(ind)
        xy_qc_sym.cx(ind, total_num_qubits - 1 - ind)
    xy_qc_sym.barrier()

    # auxiliary qubits in |+>
    aux_start = num_qubits
    aux_end = total_num_qubits - num_qubits
    for ind in range(aux_start, aux_end):
        xy_qc_sym.h(ind)
    xy_qc_sym.barrier()

    # Trotter steps
    for _ in range(trotter_steps):
        for ind, interaction in enumerate(interaction_list):
            xy_qc_sym.append(ZZZ, [interaction[0], interaction[1], num_qubits + ind])
            xy_qc_sym.append(ZXX, [interaction[0], interaction[1], num_qubits + ind + len(interaction_list)])
        xy_qc_sym.barrier()

    # rotate auxiliary qubits
    for ind in range(aux_start, aux_end):
        xy_qc_sym.sx(ind)

    # build circuit Hamiltonian Pauli strings
    ZZZ_list = []
    ZXX_list = []
    for ind, interaction in enumerate(interaction_list):
        h_term = list("I" * total_num_qubits)
        h_term[interaction[0]] = "Z"
        h_term[interaction[1]] = "Z"
        h_term[num_qubits + ind] = "Z"
        ZZZ_list.append("".join(h_term)[::-1])

    for ind, interaction in enumerate(interaction_list):
        h_term = list("I" * total_num_qubits)
        h_term[interaction[0]] = "X"
        h_term[interaction[1]] = "X"
        h_term[num_qubits + ind + len(interaction_list)] = "Z"
        ZXX_list.append("".join(h_term)[::-1])

    hamiltonian_terms = ZZZ_list + ZXX_list
    hamiltonian_terms_coeff = [γ] * len(hamiltonian_terms)

    xy_hamiltonian = SparsePauliOp(hamiltonian_terms, coeffs=hamiltonian_terms_coeff)
    xy_hamiltonian = xy_hamiltonian.assign_parameters({γ: pi / 4})
    xy_hamiltonian = SparsePauliOp(
        xy_hamiltonian.paulis,
        coeffs=np.asarray(xy_hamiltonian.coeffs, dtype=np.complex128),
    )

    return {
        "qc": xy_qc_sym,
