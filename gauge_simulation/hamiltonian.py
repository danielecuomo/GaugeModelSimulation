from __future__ import annotations

import numpy as np
from numpy import pi
from scipy.optimize import curve_fit

from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp

from .core import (
    kron_matrices,
    thermal_average_H_H2,
    compute_free_energy_and_entropy,
    exp_func,
)


class GaugeModel:
    """
    Unified container for a Generalized TFIM or XY model.

    Exact diagonalization and thermodynamic quantities
    are intrinsic model properties.

    Approximate / backend-based methods are handled by Engine.
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

        self.qc: QuantumCircuit = sim_data["qc"]
        self.H_op: SparsePauliOp = sim_data["hamiltonian"]
        self.H_sq_op: SparsePauliOp = sim_data["hamiltonian_sq"]

        self.params = {
            "beta": sim_data["beta_param"],
            "gamma": sim_data["gamma_param"],
            "eta": sim_data.get("eta_param"),
        }

        self.hamiltonian_terms_0 = sim_data["hamiltonian_terms_0"]
        self.hamiltonian_terms_1 = sim_data["hamiltonian_terms_1"]

        self._lambda_n: np.ndarray | None = None

        # thermal cache
        self._last_betas: np.ndarray | None = None
        self._last_thermal: np.ndarray | None = None
        self._last_variance: np.ndarray | None = None

    # ------------------------------------------------------------------
    # Exact Matrix & Spectrum
    # ------------------------------------------------------------------

    def _build_exact_matrix(self) -> np.ndarray:
        gamma = pi / 4
        eta = pi / 4 if self.model == "tfim" else None

        ZZ_terms = [
            term[::-1][: self.num_qubits]
            for term in self.hamiltonian_terms_0
        ]
        list_1 = [
            term[::-1][: self.num_qubits]
            for term in self.hamiltonian_terms_1
        ]

        if self.model == "tfim":
            X_terms = list_1
            H_exact = (
                gamma * sum(kron_matrices(p) for p in ZZ_terms)
                + eta * sum(kron_matrices(p) for p in X_terms)
            )
        else:
            XX_terms = list_1
            H_exact = (
                gamma * sum(kron_matrices(p) for p in ZZ_terms)
                + gamma * sum(kron_matrices(p) for p in XX_terms)
            )

        return H_exact

    def _compute_spectrum(self):
        if self._lambda_n is None:
            H_exact = self._build_exact_matrix()
            self._lambda_n, _ = np.linalg.eigh(H_exact)

    # ------------------------------------------------------------------
    # Exact Thermal Core (cached)
    # ------------------------------------------------------------------

    def _compute_thermal_core(self, beta_vals: np.ndarray):
        if (
            self._last_betas is not None
            and np.array_equal(beta_vals, self._last_betas)
        ):
            return

        self._compute_spectrum()

        H_H2 = [
            thermal_average_H_H2(beta, self._lambda_n)
            for beta in beta_vals
        ]

        thermal = np.array([h for h, _ in H_H2])
        h_sq = np.array([h2 for _, h2 in H_H2])
        variance = h_sq - thermal**2

        self._last_betas = beta_vals
        self._last_thermal = thermal
        self._last_variance = variance

    # ------------------------------------------------------------------
    # Exact Public API (single-purpose)
    # ------------------------------------------------------------------

    def exact_ground_state(self) -> float:
        self._compute_spectrum()
        return float(self._lambda_n[0])

    def exact_thermal_average(
        self,
        beta_vals: np.ndarray,
    ):
        self._compute_thermal_core(beta_vals)
        return self._last_thermal

    def exact_variance(
        self,
        beta_vals: np.ndarray,
    ):
        self._compute_thermal_core(beta_vals)
        return self._last_variance

    def exact_free_energy(
        self,
        beta_vals: np.ndarray,
    ):
        self._compute_thermal_core(beta_vals)

        if len(beta_vals) > 1:
            free_energy, _ = compute_free_energy_and_entropy(
                beta_vals,
                self._last_thermal,
            )
            return free_energy

        return np.array([])

    def exact_entropy(
        self,
        beta_vals: np.ndarray,
    ):
        self._compute_thermal_core(beta_vals)

        if len(beta_vals) > 1:
            _, entropy = compute_free_energy_and_entropy(
                beta_vals,
                self._last_thermal,
            )
            return entropy

        return np.array([])

    def exact_exponential_fit(
        self,
        beta_vals: np.ndarray,
    ):
        self._compute_thermal_core(beta_vals)

        if len(beta_vals) > 0:
            params, _ = curve_fit(
                exp_func,
                beta_vals,
                self._last_thermal,
                p0=[1, 1, 0],
            )
            return exp_func(beta_vals, *params)

        return np.array([])
