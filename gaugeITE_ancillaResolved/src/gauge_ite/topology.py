"""Lattice topology and :math:`Z_2` cycle-flux classification."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import prod
from typing import Iterable, Sequence


Edge = tuple[int, int]
Cycle = tuple[int, ...]


def _edge(u: int, v: int) -> Edge:
    return (u, v) if u < v else (v, u)


@dataclass(frozen=True)
class RectangularLattice:
    """Open or periodic rectangular nearest-neighbour lattice.

    Open lattices use elementary plaquettes as the independent cycle basis.
    Periodic lattices use all but one plaquette plus two non-contractible
    Wilson loops, matching the classification in the manuscript appendix.
    """

    rows: int
    cols: int
    boundary: str = "open"
    vertices: tuple[int, ...] = field(init=False)
    edges: tuple[Edge, ...] = field(init=False)
    plaquettes: tuple[Cycle, ...] = field(init=False)
    cycle_basis: tuple[Cycle, ...] = field(init=False)
    cycle_labels: tuple[str, ...] = field(init=False)

    def __post_init__(self) -> None:
        if self.rows <= 0 or self.cols <= 0:
            raise ValueError("rows and cols must be positive integers.")
        boundary = self.boundary.lower()
        if boundary not in {"open", "periodic"}:
            raise ValueError("boundary must be 'open' or 'periodic'.")
        if boundary == "periodic" and (self.rows < 3 or self.cols < 3):
            raise ValueError(
                "Periodic lattices require at least 3 rows and 3 columns so "
                "the undirected interaction graph has no duplicate edges."
            )
        object.__setattr__(self, "boundary", boundary)

        vertices = tuple(range(self.rows * self.cols))
        edges = self._build_edges()
        edge_index = {edge: index for index, edge in enumerate(edges)}
        plaquettes = self._build_plaquettes(edge_index)

        if boundary == "open":
            cycle_basis = plaquettes
            cycle_labels = tuple(f"p{index}" for index in range(len(plaquettes)))
        else:
            horizontal = tuple(
                edge_index[_edge(col, (col + 1) % self.cols)]
                for col in range(self.cols)
            )
            vertical = tuple(
                edge_index[
                    _edge(row * self.cols, ((row + 1) % self.rows) * self.cols)
                ]
                for row in range(self.rows)
            )
            cycle_basis = (*plaquettes[:-1], horizontal, vertical)
            cycle_labels = (
                *(f"p{index}" for index in range(len(plaquettes) - 1)),
                "W_x",
                "W_y",
            )

        expected_cycles = len(edges) - len(vertices) + 1
        if len(cycle_basis) != expected_cycles:
            raise RuntimeError(
                "Internal cycle-basis construction failed: expected "
                f"{expected_cycles}, constructed {len(cycle_basis)}."
            )

        object.__setattr__(self, "vertices", vertices)
        object.__setattr__(self, "edges", edges)
        object.__setattr__(self, "plaquettes", plaquettes)
        object.__setattr__(self, "cycle_basis", tuple(cycle_basis))
        object.__setattr__(self, "cycle_labels", tuple(cycle_labels))

    def _build_edges(self) -> tuple[Edge, ...]:
        if self.boundary == "open" and (self.rows == 1 or self.cols == 1):
            return tuple((index, index + 1) for index in range(self.rows * self.cols - 1))

        edges: list[Edge] = []
        seen: set[Edge] = set()
        for vertex in range(self.rows * self.cols):
            row, col = divmod(vertex, self.cols)
            neighbours: list[int] = []
            if col + 1 < self.cols:
                neighbours.append(vertex + 1)
            elif self.boundary == "periodic":
                neighbours.append(row * self.cols)
            if row + 1 < self.rows:
                neighbours.append(vertex + self.cols)
            elif self.boundary == "periodic":
                neighbours.append(col)

            for neighbour in neighbours:
                edge = _edge(vertex, neighbour)
                if edge not in seen:
                    seen.add(edge)
                    edges.append(edge)
        return tuple(edges)

    def _build_plaquettes(self, edge_index: dict[Edge, int]) -> tuple[Cycle, ...]:
        row_range = range(self.rows if self.boundary == "periodic" else self.rows - 1)
        col_range = range(self.cols if self.boundary == "periodic" else self.cols - 1)
        plaquettes = []
        for row in row_range:
            for col in col_range:
                top_left = row * self.cols + col
                top_right = row * self.cols + (col + 1) % self.cols
                bottom_left = ((row + 1) % self.rows) * self.cols + col
                bottom_right = (
                    ((row + 1) % self.rows) * self.cols
                    + (col + 1) % self.cols
                )
                boundary_edges = (
                    _edge(top_left, top_right),
                    _edge(top_right, bottom_right),
                    _edge(bottom_left, bottom_right),
                    _edge(top_left, bottom_left),
                )
                plaquettes.append(tuple(edge_index[edge] for edge in boundary_edges))
        return tuple(plaquettes)

    @property
    def num_vertices(self) -> int:
        return len(self.vertices)

    @property
    def num_edges(self) -> int:
        return len(self.edges)

    @property
    def cycle_rank(self) -> int:
        return len(self.cycle_basis)

    @property
    def num_classes(self) -> int:
        return 2**self.cycle_rank

    def _validate_signs(self, bond_signs: Sequence[int]) -> tuple[int, ...]:
        signs = tuple(int(sign) for sign in bond_signs)
        if len(signs) != self.num_edges:
            raise ValueError(f"Expected {self.num_edges} bond signs, received {len(signs)}.")
        if any(sign not in (-1, 1) for sign in signs):
            raise ValueError("Bond signs must contain only +1 and -1.")
        return signs

    def fluxes_for_cycles(
        self, bond_signs: Sequence[int], cycles: Iterable[Cycle]
    ) -> tuple[int, ...]:
        signs = self._validate_signs(bond_signs)
        return tuple(prod(signs[index] for index in cycle) for cycle in cycles)

    def cycle_fluxes(self, bond_signs: Sequence[int]) -> tuple[int, ...]:
        return self.fluxes_for_cycles(bond_signs, self.cycle_basis)

    def plaquette_fluxes(self, bond_signs: Sequence[int]) -> tuple[int, ...]:
        return self.fluxes_for_cycles(bond_signs, self.plaquettes)

    def class_id(self, bond_signs: Sequence[int]) -> int:
        return sum(
            (flux == -1) << index
            for index, flux in enumerate(self.cycle_fluxes(bond_signs))
        )

    def format_fluxes(self, fluxes: Sequence[int]) -> str:
        fluxes = tuple(int(value) for value in fluxes)
        if not fluxes:
            return "no cycles"
        if self.boundary == "open" and len(fluxes) == (self.rows - 1) * (self.cols - 1):
            symbols = ["+" if value == 1 else "-" for value in fluxes]
            width = max(self.cols - 1, 1)
            return "/".join(
                "".join(symbols[start : start + width])
                for start in range(0, len(symbols), width)
            )
        return " ".join(
            f"{label}={'+' if value == 1 else '-'}"
            for label, value in zip(self.cycle_labels, fluxes)
        )

    def find_vertex_gauge(
        self,
        source_bonds: Sequence[int],
        target_bonds: Sequence[int],
    ) -> tuple[int, ...]:
        """Find vertex signs ``g`` with target = g_u source g_v.

        Raises ``ValueError`` when the two configurations lie in different
        cycle-flux classes.
        """

        source = self._validate_signs(source_bonds)
        target = self._validate_signs(target_bonds)
        if self.cycle_fluxes(source) != self.cycle_fluxes(target):
            raise ValueError("Bond configurations are not gauge equivalent.")

        adjacency: dict[int, list[tuple[int, int]]] = {
            vertex: [] for vertex in self.vertices
        }
        relative = tuple(t * s for s, t in zip(source, target))
        for edge_index, (u, v) in enumerate(self.edges):
            adjacency[u].append((v, edge_index))
            adjacency[v].append((u, edge_index))

        gauges: list[int | None] = [None] * self.num_vertices
        gauges[0] = 1
        stack = [0]
        while stack:
            u = stack.pop()
            for v, edge_index in adjacency[u]:
                candidate = int(gauges[u]) * relative[edge_index]
                if gauges[v] is None:
                    gauges[v] = candidate
                    stack.append(v)
                elif gauges[v] != candidate:
                    raise ValueError("Inconsistent gauge relation on a graph cycle.")
        return tuple(int(value) for value in gauges)

