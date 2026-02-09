# gauge_simulation/hamiltonian.py

from __future__ import annotations

import numpy as np
import math
from numpy import pi
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
# EstimatorV2 is used for the simulation runs
from qiskit_ibm_runtime import EstimatorV2 as Estimator 
from qiskit.transpiler import generate_preset_pass_manager
from qiskit.quantum_info import SparsePauliOp

from .core import kron_matrices, thermal_average_H_H2, compute_free_energy_and_entropy, exp_func


class GaugeModel:
    """
    Unified container for a Generalized TFIM or XY model, holding both the 
    exact Hamiltonian and the Qiskit-based simulation circuit.

    Attributes
    ----------
    model : str
        'tfim' or 'xy'.
    lattice_size : list[int]
        [rows, cols] defining the lattice geometry.
    trotter_steps : int
        Number of Trotter steps used in the circuit.
    num_qubits : int
        Number of system qubits (rows * cols).
    qc : QuantumCircuit
        The parameterized imaginary-time evolution (ITE) circuit.
    H_op : SparsePauliOp
        The Qiskit Hamiltonian observable for $\langle H \rangle$.
    H_sq_op : SparsePauliOp
        The Qiskit Hamiltonian squared observable for $\langle H^2 \rangle$.
    params : dict
        Qiskit Parameter objects (beta, gamma, eta) for assignment.
    """

    def __init__(
        self,
        model: str,
        lattice_size: list[int],
        trotter_steps: int,
        sim_data: dict,
    ):
        self.model = model.lower()
        self.lattice_size = lattice_size
        self.trotter_steps = trotter_steps
        self.num_qubits = lattice_size[0] * lattice_size[1]
        
        # Qiskit simulation data
        self.qc: QuantumCircuit = sim_data['qc']
        self.H_op: SparsePauliOp = sim_data['hamiltonian']
        self.H_sq_op: SparsePauliOp = sim_data['hamiltonian_sq']
        self.params = {
            'beta': sim_data['beta_param'],
            'gamma': sim_data['gamma_param'],
            'eta': sim_data.get('eta_param'), # Only for TFIM
        }
        
        # Hamiltonian terms (for exact diagonalization)
        self.hamiltonian_terms_0 = sim_data['hamiltonian_terms_0']
        self.hamiltonian_terms_1 = sim_data['hamiltonian_terms_1']
        
        self._lambda_n: np.ndarray | None = None

    # ---------- Internal Helpers ----------

    def _build_exact_matrix(self) -> np.ndarray:
        """
        Rebuilds the exact $2^{N} \times 2^{N}$ Hamiltonian matrix for diagonalization, 
        using the coefficients $\gamma = \pi/4$ and $\eta = \pi/4$.
        """
        gamma = pi / 4
        eta = pi / 4 if self.model == "tfim" else None

        # Truncate Pauli strings to only cover system qubits (not auxiliary)
        ZZ_list_exact = [term[::-1][: self.num_qubits] for term in self.hamiltonian_terms_0]
        list_1_exact = [term[::-1][: self.num_qubits] for term in self.hamiltonian_terms_1]

        if self.model == "tfim":
            X_list_exact = list_1_exact
            H_exact = (
                gamma * sum(kron_matrices(pauli_string) for pauli_string in ZZ_list_exact)
                + eta * sum(kron_matrices(pauli_string) for pauli_string in X_list_exact)
            )
        else: # xy
            XX_list_exact = list_1_exact
            H_exact = (
                gamma * sum(kron_matrices(pauli_string) for pauli_string in ZZ_list_exact)
                + gamma * sum(kron_matrices(pauli_string) for pauli_string in XX_list_exact)
            )

        return H_exact

    def _compute_spectrum(self):
        """Computes and caches the eigenvalues (spectrum) of the exact Hamiltonian."""
        if self._lambda_n is None:
            H_exact = self._build_exact_matrix()
            # lambda_n is the array of eigenvalues sorted in ascending order
            self._lambda_n, _ = np.linalg.eigh(H_exact)

    # ---------- Exact Calculation API ----------

    def exact_ground_state(self) -> float:
        """Ground-state energy $\lambda_0$ (E_min)."""
        self._compute_spectrum()
        return float(self._lambda_n[0])

    def get_exact_thermal_quantities(self, beta_vals: np.ndarray) -> dict:
        """
        Calculates the exact thermal average $\langle H \rangle$, variance, 
        free energy, and entropy for a range of inverse temperatures $\beta$.
        
        Parameters
        ----------
        beta_vals : np.ndarray
            Array of inverse temperatures ($\beta$).
        
        Returns
        -------
        dict
            Contains keys: 'betas', 'thermal_avgs', 'variance', 'free_energy', 'entropy', 
            'free_energy_betas', 'entropy_betas'.
        """
        self._compute_spectrum()
        
        # Calculate $\langle H \rangle$ and $\langle H^2 \rangle$ for all beta values
        H_H2_avgs = [thermal_average_H_H2(beta, self._lambda_n) for beta in beta_vals]
        
        exact_thermal_H_avgs = np.array([H for H, H2 in H_H2_avgs])
        exact_H_sq_avgs = np.array([H2 for H, H2 in H_H2_avgs])
        
        # Calculate Variance $\text{Var}(H) = \langle H^2 \rangle - \langle H \rangle^2$
        exact_variance = exact_H_sq_avgs - (exact_thermal_H_avgs**2)
        
        # Calculate Free Energy $F(\beta)$ and Entropy $S(\beta)$
        # Need at least two points to calculate the integral for Free Energy/Entropy
        if len(beta_vals) > 1:
            exact_free_energy, exact_entropy = compute_free_energy_and_entropy(beta_vals, exact_thermal_H_avgs)
        else:
            exact_free_energy = np.array([])
            exact_entropy = np.array([])

        # Get curve fitting for $\langle H \rangle$
        if len(beta_vals) > 0:
            fit_exact_params, _ = curve_fit(exp_func, beta_vals, exact_thermal_H_avgs, p0=[1, 1, 0])
            fitted_exact_avgs = exp_func(beta_vals, *fit_exact_params)
        else:
            fitted_exact_avgs = np.array([])
            
        return {
            'betas': beta_vals,
            'thermal_avgs': exact_thermal_H_avgs,
            'variance': exact_variance,
            'free_energy': exact_free_energy,
            'entropy': exact_entropy,
            'fitted_exact_avgs': fitted_exact_avgs,
            'free_energy_betas': beta_vals[1:], # $\beta$ values for $F(\beta)$ and $S(\beta)$ (starts at second point !=0)
            'entropy_betas': beta_vals[1:],
        }


    # ---------- Simulation API ----------

    def get_approximate_thermal_quantities(self, beta_vals: np.ndarray, backend=AerSimulator(method='statevector')) -> dict:
        """
        Runs the imaginary-time evolution (ITE) simulation using Qiskit's Estimator 
        for a range of inverse temperatures $\beta$.
        
        Parameters
        ----------
        beta_vals : np.ndarray
            Array of inverse temperatures ($\beta$).
        backend : AerSimulator | Qiskit backend to run the simulation on.
            
        Returns
        -------
        dict
            Contains keys: 'betas', 'thermal_avgs', 'variance', 'free_energy', 'entropy', 
            'fitted_thermal_avgs', 'free_energy_betas', 'entropy_betas'.
        """
        
        # 1. Transpile circuits for all beta values
        pm = generate_preset_pass_manager(backend=backend, optimization_level=3)
        transpile_result = []
        
        beta_param = self.params['beta']
        fixed_params = {
            self.params['gamma']: pi / 4,
        }
        if self.model == "tfim":
            fixed_params[self.params['eta']] = pi / 4

        for val in beta_vals:
            # Assign parameters: $\gamma = \pi/4$, $\eta = \pi/4$ (if TFIM), $\beta = -\text{val}$
            psi = self.qc.assign_parameters({**fixed_params, beta_param: -val})
            isa_psi = pm.run(psi)
            
            # Apply layout to observables
            isa_hamiltonian = self.H_op.apply_layout(isa_psi.layout)
            isa_hamiltonian_sq = self.H_sq_op.apply_layout(isa_psi.layout)
            
            transpile_result.append({
                'beta_val': val, 
                'isa_psi': isa_psi, 
                'isa_hamiltonian': isa_hamiltonian, 
                'isa_hamiltonian_sq': isa_hamiltonian_sq
            })

        # 2. Estimate expectation values
        estimator = Estimator(mode=backend, options={'default_shots': 10000})

        estimated_store = []
        for item in transpile_result:
            job = estimator.run([(item['isa_psi'], item['isa_hamiltonian'])]).result()[0]
            job_sq = estimator.run([(item['isa_psi'], item['isa_hamiltonian_sq'])]).result()[0]
            expectation_val = job.data.evs
            expectation_val_sq = job_sq.data.evs
            expectation_val_std = job.data.stds
            estimated_variance = expectation_val_sq - (expectation_val**2)
            estimated_store.append((item['beta_val'], expectation_val, expectation_val_sq, estimated_variance, expectation_val_std))

        final_betas, estimated_thermal_avgs, estimated_thermal_avgs_sq, estimated_variance, expectation_val_std = zip(*estimated_store)
        final_betas, estimated_thermal_avgs, expectation_val_std = np.array(final_betas), np.array(estimated_thermal_avgs), np.array(expectation_val_std)



        # Calculate Free Energy $F(\beta)$ and Entropy $S(\beta)$     
        if len(final_betas) > 1:
            sim_free_energy, sim_entropy = compute_free_energy_and_entropy(final_betas, estimated_thermal_avgs)
        else:
            sim_free_energy = np.array([])
            sim_entropy = np.array([])
        
        # Get curve fitting for $\langle H \rangle$
        if len(final_betas) > 1:
            fit_sim_params, _ = curve_fit(exp_func, final_betas, estimated_thermal_avgs, p0=[1, 1, 0])
            fitted_estimated_avgs = exp_func(final_betas, *fit_sim_params)
        else:
            fitted_estimated_avgs = np.array([])

        return {
            'betas': final_betas,
            'thermal_avgs': estimated_thermal_avgs,
            'variance': estimated_variance,
            'free_energy': sim_free_energy,
            'entropy': sim_entropy,
            'fitted_estimated_avgs': fitted_estimated_avgs,
            'free_energy_betas': final_betas[1:],
            'entropy_betas': final_betas[1:],
            'quantum_circuit': self.qc,
        }
    
    # ---------- Get Approximate Average ----------
    def approximate_thermal_average(self, beta: float) -> float:
        """
        Computes the simulated thermal average <H> for a given beta.
        If beta is infinity, it automatically uses the optimal simulation 
        beta value (1.4) to approximate the ground state.
        """

        effective_beta = 1.4 if math.isinf(beta) or beta > 1.4 else beta
        beta_array = np.array([effective_beta])
        results = self.get_approximate_thermal_quantities(beta_array)
        
        return float(results['thermal_avgs'][0])
    
    def approximate_ground_state(self) -> float:
        """
        Calls approximate_thermal_average with beta = infinity
        
        """

        return self.approximate_thermal_average(math.inf)    
        

    # ---------- Plot QuantumCircuit ----------

    def plot_circuit(self, filename: str = "gauge_circuit.png", **draw_kwargs) -> str:
        """
        Draws the parameterized imaginary-time evolution (ITE) quantum circuit.

        Parameters
        ----------
        filename : str, optional
            The file name to save the image to. Defaults to "gauge_circuit.png".
        **draw_kwargs :
            Keyword arguments passed directly to QuantumCircuit.draw('mpl', ...).
            Common arguments include: scale, fold, cregbundle.

        Returns
        -------
        str
            The filename of the saved circuit diagram.
        """
        # Define default drawing arguments
        default_kwargs = {
            'output': 'mpl',
            'scale': 0.5, 
            'fold': 1000,
        }
        
        # Merge defaults with user-provided arguments
        merged_kwargs = {**default_kwargs, **draw_kwargs}
        
        # Draw the circuit and save it
        self.qc.draw(**merged_kwargs)
        plt.show()
        
        return None