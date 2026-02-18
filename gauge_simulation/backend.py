from __future__ import annotations

from abc import ABC, abstractmethod

from qiskit_aer import AerSimulator
from qiskit.transpiler import generate_preset_pass_manager
from qiskit_ibm_runtime import EstimatorV2 as Estimator


# ============================================================
# Abstract Backend
# ============================================================

class Backend(ABC):
    """
    Abstract execution backend.

    By default, instantiates a local AerSimulator.

    Subclasses (e.g. IonQBackend) override backend-specific behavior.
    """

    def __init__(self, shots: int = 10000):
        self.shots = shots
        self._backend = self._create_backend()

    # --------------------------------------------------------

    @abstractmethod
    def _create_backend(self):
        """Create and return the provider-specific backend object."""
        pass

    # --------------------------------------------------------

    def pass_manager(self):
        """Return preset pass manager for this backend."""
        return generate_preset_pass_manager(
            backend=self._backend,
            optimization_level=3,
        )

    # --------------------------------------------------------

    def estimator(self):
        """
        Return Estimator-like interface.

        Default implementation works for AerSimulator.
        Subclasses may override.
        """
        return Estimator(
            mode=self._backend,
            options={"default_shots": self.shots},
        )

    # --------------------------------------------------------

    @property
    def backend(self):
        """Return underlying provider backend."""
        return self._backend


# ============================================================
# Local Backend (Default)
# ============================================================


class LocalBackend(Backend):

    def _create_backend(self):
        return AerSimulator(method="statevector")

    def estimator(self):
        return Estimator(
            mode=self._backend,
            options={"default_shots": self.shots},
        )


# ============================================================
# IonQ Backend
# ============================================================

class IonQBackend(Backend):
    """
    IonQ backend via qiskit_ionq provider.

    Requires:
        pip install qiskit-ionq
    """

    def __init__(self, token: str, backend_name: str, shots: int = 10000):
        self.token = token
        self.backend_name = backend_name
        super().__init__(shots=shots)

    # --------------------------------------------------------

    def _create_backend(self):
        from qiskit_ionq import IonQProvider

        provider = IonQProvider(token=self.token)
        return provider.get_backend(self.backend_name)

    # --------------------------------------------------------

    def estimator(self):
        """
        IonQ does not support EstimatorV2 natively.

        You must implement expectation value extraction
        from sampled measurement results.

        For now, raise explicit error.
        """
        raise NotImplementedError(
            "IonQBackend does not support EstimatorV2. "
            "You must implement a custom estimator adapter."
        )
