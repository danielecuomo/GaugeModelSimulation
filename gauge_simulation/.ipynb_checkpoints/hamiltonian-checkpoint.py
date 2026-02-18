from __future__ import annotations

import numpy as np
import math
import warnings
from numpy import pi
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit.transpiler import generate_preset_pass_manager
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
    Deprecated computational interface — use Engine instead.
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

    # ------------------------------------------------------------------
    # Internal Helpers
    # ------------------------------------------------------------------

    def _build_exact_matrix(self) -> np.ndarray:
        gamma = pi / 4
        eta = pi / 4 if self.model == "tfim" else None

        ZZ_list_exact = [
            term[::-1][: self.num_qubits]
            for term in self.hamiltonian_terms_0
        ]
        list_1_exact = [
            term[::-1][: self.num_qubits]
            for term in self.hamiltonian_terms_1
        ]

        if self.model == "tfim":
            X_list_exact = list_1_exact
            H_exact = (
                gamma * sum(kron_matrices(p) for p in ZZ_list_exact)
                + eta * sum(kron_matrices(p) for p in X_list_exact)
            )
        else:
            XX_list_exact = list_1_exact
            H_exact = (
                gamma * sum(kron_matrices(p) for p in ZZ_list_exact)
                + gamma * sum(kron_matrices(p) for p in XX_list_exact)
            )

        return H_exact

    def _compute_spectrum(self):
        if self._lambda_n is None:
            H_exact = self._build_exact_matrix()
            self._lambda_n, _ = np.linalg.eigh(H_exact)

    # ------------------------------------------------------------------
    # Exact API (DEPRECATED)
    # ------------------------------------------------------------------

    def exact_ground_state(self) -> float:
        warnings.warn(
            "GaugeModel.exact_ground_state() is deprecated. "
            "Use Engine.exact_ground_state(model) instead.",
            DeprecationWarning,
            stacklevel=2,
        )

        self._compute_spectrum()
        return float(self._lambda_n[0])

    def get_exact_thermal_quantities(
        self,
        beta_vals: np.ndarray,
    ) -> dict:
        warnings.warn(
            "GaugeModel.get_exact_thermal_quantities() is deprecated. "
            "Use Engine.get_exact_thermal_quantities(model, beta_vals) instead.",
            DeprecationWarning,
            stacklevel=2,
        )

        self._compute_spectrum()

        H_H2_avgs = [
            thermal_average_H_H2(beta, self._lambda_n)
            for beta in beta_vals
        ]

        exact_thermal_H_avgs = np.array([H for H, _ in H_H2_avgs])
        exact_H_sq_avgs = np.array([H2 for _, H2 in H_H2_avgs])
        exact_variance = exact_H_sq_avgs - (exact_thermal_H_avgs ** 2)

        if len(beta_vals) > 1:
            exact_free_energy, exact_entropy = (
                compute_free_energy_and_entropy(
                    beta_vals,
                    exact_thermal_H_avgs,
                )
            )
        else:
            exact_free_energy = np.array([])
            exact_entropy = np.array([])

        if len(beta_vals) > 0:
            fit_exact_params, _ = curve_fit(
                exp_func,
                beta_vals,
                exact_thermal_H_avgs,
                p0=[1, 1, 0],
            )
            fitted_exact_avgs = exp_func(
                beta_vals,
                *fit_exact_params,
            )
        else:
            fitted_exact_avgs = np.array([])

        return {
            "betas": beta_vals,
            "thermal_avgs": exact_thermal_H_avgs,
            "variance": exact_variance,
            "free_energy": exact_free_energy,
            "entropy": exact_entropy,
            "fitted_exact_avgs": fitted_exact_avgs,
            "free_energy_betas": beta_vals[1:],
            "entropy_betas": beta_vals[1:],
        }

    # ------------------------------------------------------------------
    # Simulation API (DEPRECATED)
    # ------------------------------------------------------------------

    def get_approximate_thermal_quantities(
        self,
        beta_vals: np.ndarray,
        backend=AerSimulator(method="statevector"),
    ) -> dict:
        warnings.warn(
            "GaugeModel.get_approximate_thermal_quantities() is deprecated. "
            "Use Engine.get_approximate_thermal_quantities(model, beta_vals) instead.",
            DeprecationWarning,
            stacklevel=2,
        )

        pm = generate_preset_pass_manager(
            backend=backend,
            optimization_level=3,
        )

        beta_param = self.params["beta"]
        fixed_params = {self.params["gamma"]: pi / 4}

        if self.model == "tfim":
            fixed_params[self.params["eta"]] = pi / 4

        transpile_result = []

        for val in beta_vals:
            psi = self.qc.assign_parameters(
                {**fixed_params, beta_param: -val}
            )
            isa_psi = pm.run(psi)

            isa_hamiltonian = self.H_op.apply_layout(
                isa_psi.layout
            )
            isa_hamiltonian_sq = self.H_sq_op.apply_layout(
                isa_psi.layout
            )

            transpile_result.append(
                (
                    val,
                    isa_psi,
                    isa_hamiltonian,
                    isa_hamiltonian_sq,
                )
            )

        estimator = Estimator(
            mode=backend,
            options={"default_shots": 10000},
        )

        estimated_store = []

        for val, isa_psi, isa_H, isa_H_sq in transpile_result:
            job = estimator.run([(isa_psi, isa_H)]).result()[0]
            job_sq = estimator.run([(isa_psi, isa_H_sq)]).result()[0]

            ev = job.data.evs
            ev_sq = job_sq.data.evs
            var = ev_sq - (ev ** 2)

            estimated_store.append((val, ev, var))

        final_betas, estimated_thermal_avgs, estimated_variance = zip(
            *estimated_store
        )

        final_betas = np.array(final_betas)
        estimated_thermal_avgs = np.array(estimated_thermal_avgs)
        estimated_variance = np.array(estimated_variance)

        if len(final_betas) > 1:
            sim_free_energy, sim_entropy = (
                compute_free_energy_and_entropy(
                    final_betas,
                    estimated_thermal_avgs,
                )
            )
        else:
            sim_free_energy = np.array([])
            sim_entropy = np.array([])

        if len(final_betas) > 1:
            fit_sim_params, _ = curve_fit(
                exp_func,
                final_betas,
                estimated_thermal_avgs,
                p0=[1, 1, 0],
            )
            fitted_estimated_avgs = exp_func(
                final_betas,
                *fit_sim_params,
            )
        else:
            fitted_estimated_avgs = np.array([])

        return {
            "betas": final_betas,
            "thermal_avgs": estimated_thermal_avgs,
            "variance": estimated_variance,
            "free_energy": sim_free_energy,
            "entropy": sim_entropy,
            "fitted_estimated_avgs": fitted_estimated_avgs,
            "free_energy_betas": final_betas[1:],
            "entropy_betas": final_betas[1:],
        }

    # ------------------------------------------------------------------
    # Convenience Wrappers (DEPRECATED)
    # ------------------------------------------------------------------

    def approximate_thermal_average(self, beta: float) -> float:
        warnings.warn(
            "GaugeModel.approximate_thermal_average() is deprecated. "
            "Use Engine.approximate_thermal_average(model, beta) instead.",
            DeprecationWarning,
            stacklevel=2,
        )

        effective_beta = (
            1.4
            if math.isinf(beta) or beta > 1.4
            else beta
        )

        results = self.get_approximate_thermal_quantities(
            np.array([effective_beta])
        )

        return float(results["thermal_avgs"][0])

    def approximate_ground_state(self) -> float:
        warnings.warn(
            "GaugeModel.approximate_ground_state() is deprecated. "
            "Use Engine.approximate_ground_state(model) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.approximate_thermal_average(math.inf)
