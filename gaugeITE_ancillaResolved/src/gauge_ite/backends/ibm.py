"""Backward-compatible IBM adapter.

New code should use :class:`gauge_ite.executors.IBMRuntimeExecutor`, which
adds account profiles, smart backend selection, submission safety, and a
provider-neutral result object.
"""

from __future__ import annotations

from ..circuits import CircuitBundle


def run_ibm_sampler(
    bundle: CircuitBundle,
    backend,
    basis: str = "ancilla",
    shots: int = 4096,
    optimization_level: int = 1,
) -> dict:
    """Compile and submit one measured circuit through IBM SamplerV2."""

    from ..executors.ibm import IBMRuntimeExecutor

    execution = IBMRuntimeExecutor(
        backend=backend,
        shots=shots,
        optimization_level=optimization_level,
        submit=True,
    ).run(bundle, basis=basis)
    return {
        "job_id": execution.job_id,
        "counts": execution.counts,
        "register_counts": execution.register_counts,
        "basis": basis,
        "shots": shots,
        "isa_circuit": execution.transpiled_circuit,
        "metadata": execution.metadata,
    }
