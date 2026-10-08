"""Physical TFIM model specifications and exact system Hamiltonians."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Sequence

import numpy as np

from .topology import RectangularLattice


def _coefficient_tuple(value: float | Sequence[float], length: int, name: str) -> tuple[float, ...]:
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        return (float(array),) * length
    if array.shape != (length,):
        raise ValueError(f"{name} must be a scalar or a length-{length} sequence.")
    return tuple(float(item) for item in array)


@dataclass(frozen=True)
class Term:
    index: int
    kind: str
    sites: tuple[int, ...]
    coefficient: float
    compiled_sign: int = 1

    @property
    def signed_coefficient(self) -> float:
        return self.coefficient * self.compiled_sign


@dataclass(frozen=True)
class TFIMSpec:
    lattice: RectangularLattice
    gamma: tuple[float, ...]
    eta: tuple[float, ...]
    compiled_bond_signs: tuple[int, ...]
    terms: tuple[Term, ...] = field(init=False)

    def __post_init__(self) -> None:
        if len(self.gamma) != self.lattice.num_edges:
            raise ValueError("gamma must contain one value per lattice edge.")
        if len(self.eta) != self.lattice.num_vertices:
            raise ValueError("eta must contain one value per lattice vertex.")
        if len(self.compiled_bond_signs) != self.lattice.num_edges:
            raise ValueError("compiled_bond_signs must contain one sign per edge.")
        if any(sign not in (-1, 1) for sign in self.compiled_bond_signs):
            raise ValueError("compiled_bond_signs must contain only +1 and -1.")

        terms: list[Term] = []
        for index, (edge, coefficient, sign) in enumerate(
            zip(self.lattice.edges, self.gamma, self.compiled_bond_signs)
        ):
            terms.append(Term(index, "ZZ", edge, coefficient, sign))
        offset = self.lattice.num_edges
        for vertex, coefficient in enumerate(self.eta):
            terms.append(Term(offset + vertex, "X", (vertex,), coefficient, 1))
        object.__setattr__(self, "terms", tuple(terms))

    @classmethod
    def rectangular(
        cls,
        rows: int,
        cols: int,
        gamma: float | Sequence[float] = np.pi / 4,
        eta: float | Sequence[float] = np.pi / 4,
        compiled_bond_signs: Sequence[int] | None = None,
        boundary: str = "open",
    ) -> "TFIMSpec":
        lattice = RectangularLattice(rows, cols, boundary=boundary)
        compiled = (
            (1,) * lattice.num_edges
            if compiled_bond_signs is None
            else tuple(int(sign) for sign in compiled_bond_signs)
        )
        return cls(
            lattice=lattice,
            gamma=_coefficient_tuple(gamma, lattice.num_edges, "gamma"),
            eta=_coefficient_tuple(eta, lattice.num_vertices, "eta"),
            compiled_bond_signs=compiled,
        )

    @property
    def num_system_qubits(self) -> int:
        return self.lattice.num_vertices

    @property
    def num_terms(self) -> int:
        return len(self.terms)

    @property
    def num_bond_terms(self) -> int:
        return self.lattice.num_edges

    def branch_hamiltonian(
        self,
        bond_signs: Sequence[int],
        field_signs: Sequence[int],
    ) -> np.ndarray:
        bonds = tuple(int(sign) for sign in bond_signs)
        fields = tuple(int(sign) for sign in field_signs)
        if len(bonds) != self.num_bond_terms:
            raise ValueError("bond_signs must contain one sign per bond term.")
        if len(fields) != self.num_system_qubits:
            raise ValueError("field_signs must contain one sign per field term.")
        if any(sign not in (-1, 1) for sign in (*bonds, *fields)):
            raise ValueError("Branch signs must contain only +1 and -1.")

        dimension = 2**self.num_system_qubits
        hamiltonian = np.zeros((dimension, dimension), dtype=complex)
        basis = np.arange(dimension, dtype=np.int64)

        for coefficient, sign, (u, v) in zip(
            self.gamma, bonds, self.lattice.edges
        ):
            z_u = 1 - 2 * ((basis >> u) & 1)
            z_v = 1 - 2 * ((basis >> v) & 1)
            hamiltonian[basis, basis] += coefficient * sign * z_u * z_v

        for vertex, (coefficient, sign) in enumerate(zip(self.eta, fields)):
            flipped = basis ^ (1 << vertex)
            hamiltonian[flipped, basis] += coefficient * sign
        return hamiltonian

    def clean_hamiltonian(self) -> np.ndarray:
        return self.branch_hamiltonian(
            self.compiled_bond_signs,
            (1,) * self.num_system_qubits,
        )


@lru_cache(maxsize=None)
def pauli_term_matrix(kind: str, sites: tuple[int, ...], num_qubits: int) -> np.ndarray:
    dimension = 2**num_qubits
    matrix = np.zeros((dimension, dimension), dtype=complex)
    basis = np.arange(dimension, dtype=np.int64)
    if kind == "ZZ":
        u, v = sites
        values = (1 - 2 * ((basis >> u) & 1)) * (1 - 2 * ((basis >> v) & 1))
        matrix[basis, basis] = values
    elif kind == "X":
        (vertex,) = sites
        matrix[basis ^ (1 << vertex), basis] = 1.0
    else:
        raise ValueError(f"Unsupported term kind: {kind!r}.")
    return matrix


def local_pauli_unitary(
    x_signs: Sequence[int],
    z_signs: Sequence[int],
) -> np.ndarray:
    """Build a system unitary from local X and Z conjugations.

    A sign of ``-1`` means apply that Pauli. Global phases are irrelevant for
    the gauge-covariance comparisons performed by the package.
    """

    x_signs = tuple(int(sign) for sign in x_signs)
    z_signs = tuple(int(sign) for sign in z_signs)
    if len(x_signs) != len(z_signs):
        raise ValueError("x_signs and z_signs must have the same length.")
    identity = np.eye(2, dtype=complex)
    x_matrix = np.array([[0, 1], [1, 0]], dtype=complex)
    z_matrix = np.array([[1, 0], [0, -1]], dtype=complex)
    local = []
    for x_sign, z_sign in zip(x_signs, z_signs):
        operator = identity
        if x_sign == -1:
            operator = x_matrix @ operator
        if z_sign == -1:
            operator = z_matrix @ operator
        local.append(operator)
    result = np.array([[1.0 + 0.0j]])
    for operator in reversed(local):
        result = np.kron(result, operator)
    return result

