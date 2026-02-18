from __future__ import annotations

import numpy as np
import math
from numpy import pi
from scipy.optimize import curve_fit

from .hamiltonian import GaugeModel
from .backend import Backend


class Engine:
    """
    Backend-driven approximation engine.

    Exact diagonalization and exact thermodynamics
    are handled by GaugeModel.

    This class only manages circuit execution and
    approximate thermodynamic quantities.
    """

    def __init__(self, backend: Backend):
        self.backend = backend

        # simulation cache
        self._last_betas = None
        self._last_thermal = None
        self._last_variance = None

    # ------------------------------------------------------------------
    # Simulation Thermal Core (cached)
    # ------------------------------------------------------------------

    def _compute_sim_thermal_core(
        self,
        H: GaugeModel,
        beta_vals: np.ndarray,
    ):
        if (
            self._last_betas is not None
            and np.array_equal(beta_vals, self._last_betas)
        ):
            return

        pm = self.backend.pass_manager()
        estimator = self.backend.estimator()

        beta_param = H.params["beta"]
        fixed_params = {H.params["gamma"]: pi / 4}

        if H.model == "tfim":
            fixed_params[H.params["eta"]] = pi / 4

        results = []

        for val in beta_vals:
            psi = H.qc.assign_parameters(
                {**fixed_params, beta_param: -val}
            )
            isa_psi = pm.run(psi)

            isa_H = H.H_op.apply_layout(isa_psi.layout)
            isa_H_sq = H.H_sq_op.apply_layout(isa_psi.layout)

            ev = estimator.run([(isa_psi, isa_H)]).result()[0].data.evs
            ev_sq = estimator.run([(isa_psi, isa_H_sq)]).result()[0].data.evs

            variance = ev_sq - ev**2
            results.append((val, ev, variance))

        betas, thermal, variance = zip(*results)

        self._last_betas = np.array(betas)
        self._last_thermal = np.array(thermal)
        self._last_variance = np.array(variance)

    # ------------------------------------------------------------------
    # Public Approximate API (single-purpose)
    # ------------------------------------------------------------------

    def approximate_thermal_average(
        self,
        H: GaugeModel,
        beta_vals: np.ndarray,
    ):
        self._compute_sim_thermal_core(H, beta_vals)
        return self._last_thermal

    def approximate_variance(
        self,
        H: GaugeModel,
        beta_vals: np.ndarray,
    ):
        self._compute_sim_thermal_core(H, beta_vals)
        return self._last_variance

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    def approximate_ground_state(
        self,
        H: GaugeModel,
        beta_cap: float = 1.4,
    ) -> float:
        """
        Approximate ground state via large-beta limit.
        """
        effective_beta = (
            beta_cap if not math.isinf(beta_cap) else 1.4
        )

        value = self.approximate_thermal_average(
            H,
            np.array([effective_beta]),
        )

        return float(value[0])
