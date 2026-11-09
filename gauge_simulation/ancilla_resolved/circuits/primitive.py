"""The original ancilla-Pauli circuit primitive."""

from __future__ import annotations

import numpy as np
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp

from ..models import Term


def ite_angle(
    beta: float,
    steps: int,
    signed_coefficient: float,
    mapping: str,
) -> float:
    """Return the Pauli-evolution angle used by the original circuit.

    The minus sign is the numeric form of binding the legacy symbolic beta
    parameter to ``-beta``.
    """

    value = -float(beta) * float(signed_coefficient) / (2 * int(steps))
    if mapping == "linear":
        return value
    if mapping == "exact_kraus":
        return float(np.arctan(np.tanh(value)))
    raise ValueError("mapping must be 'linear' or 'exact_kraus'.")


def evolution_gate_for_term(term: Term, beta: float, steps: int, mapping: str):
    angle = ite_angle(beta, steps, term.signed_coefficient, mapping)
    if term.kind == "ZZ":
        operator = SparsePauliOp("ZZZ")
    elif term.kind == "X":
        # This is the same Z ^ X ordering used by the legacy package. When
        # appended as [system, ancilla], it implements X_system Z_ancilla.
        operator = SparsePauliOp("ZX")
    else:
        raise ValueError(f"Unsupported term kind: {term.kind!r}.")
    return PauliEvolutionGate(operator, time=angle)


def append_term_interaction(
    circuit,
    term: Term,
    ancilla: int,
    beta: float,
    steps: int,
    mapping: str,
) -> None:
    gate = evolution_gate_for_term(term, beta, steps, mapping)
    if term.kind == "ZZ":
        circuit.append(gate, [term.sites[0], term.sites[1], ancilla])
    else:
        circuit.append(gate, [term.sites[0], ancilla])

