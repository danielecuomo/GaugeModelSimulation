"""Optional pytket/Quantinuum adapter.

The adapter deliberately accepts a configured backend object so credentials
remain outside the package and outside saved notebooks.
"""

from __future__ import annotations

from ..circuits import CircuitBundle, add_final_measurements


def run_quantinuum(
    bundle: CircuitBundle,
    backend,
    basis: str = "ancilla",
    shots: int = 4096,
    optimisation_level: int = 2,
) -> dict:
    try:
        from pytket.extensions.qiskit import qiskit_to_tk
    except ImportError as exc:
        raise ImportError(
            "Install pytket, pytket-qiskit, and pytket-quantinuum to use "
            "the Quantinuum adapter."
        ) from exc

    measured = add_final_measurements(bundle, basis=basis)
    tk_circuit = qiskit_to_tk(measured)
    compiled = backend.get_compiled_circuit(
        tk_circuit, optimisation_level=optimisation_level
    )
    handle = backend.process_circuit(compiled, n_shots=shots)
    result = backend.get_result(handle)
    return {
        "handle": handle,
        "counts": dict(result.get_counts()),
        "basis": basis,
        "shots": shots,
        "compiled_circuit": compiled,
    }

