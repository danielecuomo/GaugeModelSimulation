# gauge_simulation/core.py

import numpy as np
from numpy import pi
import psutil
import math
from scipy.integrate import cumulative_trapezoid

from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.circuit import Parameter


def get_nearest_neighbor_interactions(lattice_size: list[int]) -> list[tuple[int, int]]:
    """
    Generate nearest-neighbor interactions pairs for 1D chains or 2D rectangular lattices.

    Args:
        lattice_size (list[int]): A list containing [rows, cols]. 
                                  Example: [1, 5] or [5,1] for 1D, [4, 5] for 2D.

    Returns:
        list[tuple[int,int]]: list of (i,j) index pairs for nearest-neighbor couplings.
    """
    # Validate input
    if not (isinstance(lattice_size, list) and len(lattice_size) == 2):
        raise ValueError("lattice_size must be a list of two integers [rows, cols].")
    
    rows, cols = lattice_size

    if not (isinstance(rows, int) and isinstance(cols, int) and rows > 0 and cols > 0):
        raise ValueError("Both rows and cols must be positive integers.")

    num_qubits = rows * cols
    interactions = []

    # Get dimension and calculate interactions
    # Check for 1D (linear chain)
    if rows == 1 or cols == 1:
        
        interactions = [(i, i + 1) for i in range(num_qubits - 1)]

    # Check for 2D (full rectangular lattice)
    else:
        for i in range(num_qubits):
            # right neighbor 
            if (i + 1) % cols != 0:
                interactions.append((i, i + 1))

            # down neighbor
            if i + cols < num_qubits:
                interactions.append((i, i + cols))

    return interactions

def get_system_qubit_limit(memory_factor: float = 0.75) -> int:
    
    """
    Get the maximum number of qubits the local system can reasonably simulate
    based on available memory (RAM)
    
    Args:
        memory_factor (float): Fraction of the available RAM to dedicate to the 
                               simulation (setting it to 75% by default)
                               
    Returns:
        int: The maximum safe number of qubits that can be run for the simulation
    
    """
    
    # Get the available memory in bytes
    available_mem_bytes = psutil.virtual_memory().available
    
    # Get fraction of memory to use
    safe_mem_bytes = available_mem_bytes * memory_factor
    
    # Convert available memory to qubits
    # For N qubits, there are 2^N amplitudes, each amp takes about 16 bytes
    # total_memory = 16 * 2 ** N
    max_safe_qubit = math.floor(math.log2(safe_mem_bytes / 16))
    
    return max_safe_qubit


def kron_matrices(pauli_string: str) -> np.ndarray:
    """
    Compute the Kronecker product for a given Pauli string.

    Args:
        pauli_string (str): String representing the Pauli operators (e.g., 'ZZX').

    Returns:
        np.ndarray: The resulting matrix from the Kronecker product.
    """
    I_g = np.array([[1, 0], [0, 1]], dtype=complex)
    X_g = np.array([[0, 1], [1, 0]], dtype=complex)
    Z_g = np.array([[1, 0], [0, -1]], dtype=complex)
    pauli_dict = {'I': I_g, 'X': X_g, 'Z': Z_g}

    result = pauli_dict[pauli_string[0]]
    for char in pauli_string[1:]:
        result = np.kron(result, pauli_dict[char])
    return result


def thermal_average_H_H2(beta: float, lambda_n: np.ndarray) -> tuple[float, float]:
    """
    Calculate the thermal average of the Hamiltonian, and its square at a given beta.

    Args:
        beta (float): Inverse temperature.
        lambda_n (np.ndarray): Eigenvalues of the Hamiltonian.

    Returns:
        float: Thermal average <H>.
    """
    w = np.exp(-beta * lambda_n)
    Z_f = np.sum(w)
    H = float(np.sum(w * lambda_n) / Z_f)
    H2 = float(np.sum(w * lambda_n**2) / Z_f)
    return H, H2


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


def compute_free_energy_and_entropy(beta: np.ndarray, thermal_H: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the 'constant-free' free energy and entropy from thermal averages.

    Definitions used:
        F_tilde(β) = ( 1 / β ) * ∫_{β0}^{β} <H>(β') dβ'
        S_tilde(β) = β * ( <H>(β) - F_tilde(β) )

    These omit the additive constant A/β that would normally appear in the
    true free energy. This is intentional since both simulated and exact
    results are compared without the constant.

    Parameters
    ----------
    beta : np.ndarray
        1D array of β values (inverse temperature).
        Must be strictly increasing and beta[0] > 0.

    H_expect : np.ndarray
        1D array of thermal expectation values <H>(β).
        Must be the same shape as beta.

    Returns
    -------
    F_tilde : np.ndarray
        The constant-free free energy evaluated at each β.

    S_tilde : np.ndarray
        The constant-free entropy evaluated at each β.

    Notes
    -----
    - No constant term is added.
    - βeta = 0 is not supported (division by zero).
    - βeta must be sorted to exclude beta = 0.
    """
    beta_tilde = beta[1:]
    thermal_H_tilde = thermal_H[1:]

    F_tilde = 1/beta_tilde * (cumulative_trapezoid(thermal_H_tilde, beta_tilde, initial=0.0))

    S_tilde = beta_tilde * (thermal_H_tilde - F_tilde)

    return F_tilde, S_tilde


def tfim_generalized(lattice_size: list[int], trotter_steps: int) -> dict:
    """
    Implements Cuomo's guage-invarinat ITE for a transverse field Ising model in 1D or 2D.

    Args:
        lattice_size (list[int]): A list containing [rows, cols]. 
                                  Example: [1, 5] or [5,1] for 1D, [4, 5] for 2D.
        trotter_steps (int): Total number of Trotter steps to use in the circuit.

    Returns:
        dict: Dictionary containing the quantum circuit, results, and other computed values.

    Raises:
        ValueError: If lattice_size is not a list of two integers,
                    or if the total number of qubits exceeds 20 (to avoid high memory usage).
    """
    if not (isinstance(lattice_size, list) and len(lattice_size) == 2):
        raise ValueError("lattice_size must be a list of two integers [rows, cols].")
    
    rows, cols = lattice_size

    if not (isinstance(rows, int) and isinstance(cols, int) and rows > 0 and cols > 0):
        raise ValueError("Both rows and cols must be positive integers.")
    
    num_qubits = rows * cols

    if rows == 1 or cols == 1:
        dimension = 1
    else:
        dimension = 2

    zz_interaction_list = get_nearest_neighbor_interactions(lattice_size)
    total_hamiltonian_terms = num_qubits + len(zz_interaction_list)
    total_num_qubits = 2*num_qubits + total_hamiltonian_terms

    # Check for simulation feasibility; good when running on a low-end to mid-range PC
    local_system_limit = get_system_qubit_limit()
    
    if total_num_qubits > local_system_limit:
        raise ValueError(
            f"Total qubits required ({total_num_qubits}) exceeds local system limit ({local_system_limit}). "
            "Simulation requires more memory than is currently available. "
            "Consider reducing lattice_size."
        )

    # Define parameters and PauliEvolutionGates
    γ = Parameter('γ')
    η = Parameter('η')
    β = Parameter('β')
    X = SparsePauliOp('X')
    Z = SparsePauliOp('Z')
    ZZZ = PauliEvolutionGate(Z ^ Z ^ Z, time=(β / (2 * trotter_steps)) * γ)
    ZX = PauliEvolutionGate(Z ^ X, time=(β / (2 * trotter_steps)) * η)

    # Create quantum circuit
    tfim_qc_sym = QuantumCircuit(total_num_qubits)

    # Create maximally mixed state for system qubits
    for ind in range(num_qubits):
        tfim_qc_sym.h(ind)
        tfim_qc_sym.cx(ind, total_num_qubits - 1 - ind)  
    tfim_qc_sym.barrier()

    # Prepare auxiliary qubits in |+> state
    aux_start = num_qubits
    aux_end = total_num_qubits - num_qubits
    for ind in range(aux_start, aux_end):
        tfim_qc_sym.h(ind)
    tfim_qc_sym.barrier()

    # Apply evolution gates in Trotter steps
    for _ in range(trotter_steps):
        # ZZZ operations
        for ind, zz_interaction in enumerate(zz_interaction_list):
            tfim_qc_sym.append(ZZZ, [zz_interaction[0], zz_interaction[1], num_qubits + ind])
        tfim_qc_sym.barrier()

        # ZX operations
        for ind in range(num_qubits):
            tfim_qc_sym.append(ZX, [ind, num_qubits + len(zz_interaction_list) + ind])
        tfim_qc_sym.barrier()

    # Rotate auxiliary qubits to the measurement basis
    for ind in range(aux_start, aux_end):
        tfim_qc_sym.sx(ind)

    # Build circuit Hamiltonians for SparsePauliOp
    ZZZ_list = []
    ZX_list = []
    for ind, zz_interaction in enumerate(zz_interaction_list):
        h_term = list('I' * total_num_qubits)
        h_term[zz_interaction[0]], h_term[zz_interaction[1]], h_term[num_qubits + ind] = 'Z', 'Z', 'Z'
        ZZZ_list.append(''.join(h_term)[::-1])

    for ind in range(num_qubits):
        h_term = list('I' * total_num_qubits)
        h_term[ind], h_term[num_qubits + len(zz_interaction_list) + ind] = 'X', 'Z'
        ZX_list.append(''.join(h_term)[::-1])

    hamiltonian_terms = ZZZ_list + ZX_list
    hamiltonian_terms_coeff = [γ] * len(ZZZ_list) + [η] * len(ZX_list)

    # Build and parameterize Hamiltonian
    tfim_hamiltonian = SparsePauliOp(hamiltonian_terms, coeffs=hamiltonian_terms_coeff)
    tfim_hamiltonian = tfim_hamiltonian.assign_parameters({γ: pi / 4, η: pi / 4}) # using a fixed-value pi/4 for the params 

    # qiskit needs complex-typed coefficients for tfim_hamiltonian
    # only 'SparsePauliOp' with complex-typed coefficients can be converted to 'SparseObservable'
    tfim_hamiltonian = SparsePauliOp(tfim_hamiltonian.paulis, coeffs=np.asarray(tfim_hamiltonian.coeffs, dtype=np.complex128))
    tfim_hamiltonian_sq = (tfim_hamiltonian @ tfim_hamiltonian).simplify() # gets the SparsePauliOp of the square

    tfim_output = {
        'qc': tfim_qc_sym, 'hamiltonian': tfim_hamiltonian, 'hamiltonian_sq': tfim_hamiltonian_sq,
        'hamiltonian_terms_0': ZZZ_list, 'hamiltonian_terms_1': ZX_list,
        'total_num_qubits': total_num_qubits,
        'gamma_param': γ,
        'eta_param': η,
        'beta_param': β
    }

    return tfim_output

def xy_generalized(lattice_size: list[int], trotter_steps: int) -> dict:
    """
    Implements Cuomo's guage-invarinat ITE for an XY model in 1D or 2D.

    Args:
        lattice_size (list[int]): A list containing [rows, cols]. 
                                  Example: [1, 5] or [5,1] for 1D, [4, 5] for 2D.
        trotter_steps (int): Total number of Trotter steps to use in the circuit.

    Returns:
        dict: Dictionary containing the quantum circuit, results, and other computed values.

    Raises:
        ValueError: If lattice_size is not a list of two integers,
                    or if the total number of qubits exceeds 20 (to avoid high memory usage).
    """

    if not (isinstance(lattice_size, list) and len(lattice_size) == 2):
        raise ValueError("lattice_size must be a list of two integers [rows, cols].")
    
    rows, cols = lattice_size

    if not (isinstance(rows, int) and isinstance(cols, int) and rows > 0 and cols > 0):
        raise ValueError("Both rows and cols must be positive integers.")
    
    num_qubits = rows * cols

    if rows == 1 or cols == 1:
        dimension = 1
    else:
        dimension = 2

    interaction_list = get_nearest_neighbor_interactions(lattice_size)
    total_hamiltonian_terms = 2 * len(interaction_list)
    total_num_qubits = 2*num_qubits + total_hamiltonian_terms

    # Check for simulation feasibility; good when running on a low-end to mid-range PC
    local_system_limit = get_system_qubit_limit()
    
    if total_num_qubits > local_system_limit:
        raise ValueError(
            f"Total qubits required ({total_num_qubits}) exceeds local system limit ({local_system_limit}). "
            "Simulation requires more memory than is currently available. "
            "Consider reducing lattice_size."
        )
    
    # Define parameters and PauliEvolutionGates
    γ = Parameter('γ') # we are taking all the coefficients to be equal and == pi/4
    β = Parameter('β')
    X = SparsePauliOp('X')
    Z = SparsePauliOp('Z')
    ZZZ = PauliEvolutionGate(Z ^ Z ^ Z, time=(β / (2 * trotter_steps)) * γ)
    ZXX = PauliEvolutionGate(Z ^ X ^ X, time=(β / (2 * trotter_steps)) * γ)

    # Create quantum circuit
    xy_qc_sym = QuantumCircuit(total_num_qubits)

    # Create maximally mixed state for system qubits
    for ind in range(num_qubits):
        xy_qc_sym.h(ind)
        xy_qc_sym.cx(ind, total_num_qubits - 1 - ind)  
    xy_qc_sym.barrier()

    # Prepare auxiliary qubits in |+> state
    aux_start = num_qubits
    aux_end = total_num_qubits - num_qubits
    for ind in range(aux_start, aux_end):
        xy_qc_sym.h(ind)
    xy_qc_sym.barrier()

    # Apply evolution gates in Trotter steps
    for _ in range(trotter_steps):
        for ind, interaction in enumerate(interaction_list):
            xy_qc_sym.append(ZZZ, [interaction[0], interaction[1], num_qubits + ind])
            xy_qc_sym.append(ZXX, [interaction[0], interaction[1], num_qubits + ind + len(interaction_list)])

        xy_qc_sym.barrier()

    # Rotate auxiliary qubits to the measurement basis
    for ind in range(aux_start, aux_end):
        xy_qc_sym.sx(ind)

    # Build circuit Hamiltonians for SparsePauliOp
    ZZZ_list = []
    ZXX_list = []
    for ind, interaction in enumerate(interaction_list):
        h_term = list('I' * total_num_qubits)
        h_term[interaction[0]], h_term[interaction[1]], h_term[num_qubits + ind] = 'Z', 'Z', 'Z'
        ZZZ_list.append(''.join(h_term)[::-1])

    for ind, interaction in enumerate(interaction_list):
        h_term = list('I' * total_num_qubits)
        h_term[interaction[0]], h_term[interaction[1]], h_term[num_qubits + ind + len(interaction_list)]  = 'X', 'X', 'Z'
        ZZZ_list.append(''.join(h_term)[::-1])

    hamiltonian_terms = ZZZ_list + ZXX_list
    hamiltonian_terms_coeff = [γ] * len(hamiltonian_terms)

    # Build and parameterize Hamiltonian
    xy_hamiltonian = SparsePauliOp(hamiltonian_terms, coeffs=hamiltonian_terms_coeff)
    xy_hamiltonian = xy_hamiltonian.assign_parameters({γ: pi / 4})

    # qiskit needs complex-typed coefficients for xy_hamiltonian
    # only 'SparsePauliOp' with complex-typed coefficients can be converted to 'SparseObservable'
    xy_hamiltonian = SparsePauliOp(xy_hamiltonian.paulis, coeffs=np.asarray(xy_hamiltonian.coeffs, dtype=np.complex128))
    xy_hamiltonian_sq = (xy_hamiltonian @ xy_hamiltonian).simplify() # gets the SparsePauliOp of the square


    xy_output = {
        'qc': xy_qc_sym, 'hamiltonian': xy_hamiltonian, 'hamiltonian_sq': xy_hamiltonian_sq,
        'hamiltonian_terms_0': ZZZ_list, 'hamiltonian_terms_1': ZXX_list,
        'total_num_qubits': total_num_qubits,
        'gamma_param': γ,
        'beta_param': β
    }


    return xy_output