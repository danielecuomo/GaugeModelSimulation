# gauge_simulation/hamiltonian.py

from __future__ import annotations

import numpy as np
from numpy import pi

from .core import kron_matrices, thermal_average_H


class Hamiltonian:
    """
    Exact Hamiltonian built from the Pauli strings used in the circuits.

    Attributes
    ----------
    model : str
        'tfim' or 'xy'
    num_qubits : int
        Number of system qubits.
    hamiltonian_terms_0, hamiltonian_terms_1 : list[str]
        Pauli strings from the circuit construction (as in TFIM_XY_generalized.ipynb).
    """

    def __init__(
        self,
        model: str,
        num_qubits: int,
        hamiltonian_terms_0: list[str],
        hamiltonian_terms_1: list[str],
    ):
        self.model = model.lower()
        self.num_qubits = int(num_qubits)
        self.hamiltonian_terms_0 = list(hamiltonian_terms_0)
        self.hamiltonian_terms_1 = list(hamiltonian_terms_1)

        self._lambda_n: np.ndarray | None = None

    # ---------- internal helpers ----------

    def _build_exact_matrix(self) -> np.ndarray:
        """
        Rebuild H_exact using the same recipe as the notebook.
        """

        gamma = pi / 4
        eta = pi / 4 if self.model == "tfim" else None

        I_g = np.array([[1, 0], [0, 1]], dtype=complex)
        X_g = np.array([[0, 1], [1, 0]], dtype=complex)
        Z_g = np.array([[1, 0], [0, -1]], dtype=complex)
        pauli_dict = {"I": I_g, "X": X_g, "Z": Z_g}

        ZZ_list_exact = [term[::-1][: self.num_qubits] for term in self.hamiltonian_terms_0]
        list_1_exact = [term[::-1][: self.num_qubits] for term in self.hamiltonian_terms_1]

        if self.model == "tfim":
            X_list_exact = list_1_exact
            H_exact = (
                gamma * sum(kron_matrices(ps, pauli_dict) for ps in ZZ_list_exact)
                + eta * sum(kron_matrices(ps, pauli_dict) for ps in X_list_exact)
            )
        elif self.model == "xy":
            XX_list_exact = list_1_exact
            H_exact = (
                gamma * sum(kron_matrices(ps, pauli_dict) for ps in ZZ_list_exact)
                + gamma * sum(kron_matrices(ps, pauli_dict) for ps in XX_list_exact)
            )
        else:
            raise ValueError(f"Unsupported model {self.model!r}.")

        return H_exact

    def _compute_spectrum(self):
        if self._lambda_n is None:
            H_exact = self._build_exact_matrix()
            self._lambda_n, _ = np.linalg.eigh(H_exact)

    # ---------- public API ----------

    def ground_state(self) -> float:
        """
        Ground-state energy λ₀.
        """
        self._compute_spectrum()
        return float(self._lambda_n[0])

    def thermal_average(self, beta: float) -> float:
        """
        Thermal average ⟨H⟩_β.
        """
        self._compute_spectrum()
        if np.isinf(beta):
            return self.ground_state()
        return thermal_average_H(beta, self._lambda_n)
