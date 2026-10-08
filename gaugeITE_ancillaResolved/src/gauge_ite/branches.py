"""Measured sign records and gauge-class labels."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .models import TFIMSpec


def bits_to_signs(bits: Sequence[int]) -> tuple[int, ...]:
    bits = tuple(int(bit) for bit in bits)
    if any(bit not in (0, 1) for bit in bits):
        raise ValueError("bits must contain only 0 and 1.")
    return tuple(1 if bit == 0 else -1 for bit in bits)


def signs_to_bits(signs: Sequence[int]) -> tuple[int, ...]:
    signs = tuple(int(sign) for sign in signs)
    if any(sign not in (-1, 1) for sign in signs):
        raise ValueError("signs must contain only +1 and -1.")
    return tuple(0 if sign == 1 else 1 for sign in signs)


@dataclass(frozen=True)
class BranchRecord:
    measured_signs_by_step: tuple[tuple[int, ...], ...]
    effective_bond_signs_by_step: tuple[tuple[int, ...], ...]
    field_signs_by_step: tuple[tuple[int, ...], ...]
    cycle_fluxes_by_step: tuple[tuple[int, ...], ...]
    plaquette_fluxes_by_step: tuple[tuple[int, ...], ...]
    class_ids_by_step: tuple[int, ...]
    raw_bits: tuple[int, ...]

    @property
    def final_measured_signs(self) -> tuple[int, ...]:
        return self.measured_signs_by_step[-1]

    @property
    def bond_signs(self) -> tuple[int, ...]:
        return self.effective_bond_signs_by_step[-1]

    @property
    def field_signs(self) -> tuple[int, ...]:
        return self.field_signs_by_step[-1]

    @property
    def cycle_fluxes(self) -> tuple[int, ...]:
        return self.cycle_fluxes_by_step[-1]

    @property
    def plaquette_fluxes(self) -> tuple[int, ...]:
        return self.plaquette_fluxes_by_step[-1]

    @property
    def class_id(self) -> int:
        return self.class_ids_by_step[-1]

    @property
    def persistent_signs(self) -> bool:
        return all(
            signs == self.measured_signs_by_step[0]
            for signs in self.measured_signs_by_step[1:]
        )

    @property
    def ancilla_bitstring(self) -> str:
        return "".join(str(bit) for bit in self.raw_bits)


def record_from_signs(
    spec: TFIMSpec,
    measured_signs_by_step: Sequence[Sequence[int]],
) -> BranchRecord:
    steps = tuple(tuple(int(sign) for sign in row) for row in measured_signs_by_step)
    if not steps:
        raise ValueError("At least one sign step is required.")
    if any(len(row) != spec.num_terms for row in steps):
        raise ValueError(f"Each sign step must contain {spec.num_terms} entries.")
    if any(sign not in (-1, 1) for row in steps for sign in row):
        raise ValueError("Measured signs must contain only +1 and -1.")

    effective_bonds = []
    fields = []
    cycle_fluxes = []
    plaquette_fluxes = []
    class_ids = []
    for signs in steps:
        measured_bonds = signs[: spec.num_bond_terms]
        branch_bonds = tuple(
            measured * compiled
            for measured, compiled in zip(measured_bonds, spec.compiled_bond_signs)
        )
        branch_fields = signs[spec.num_bond_terms :]
        effective_bonds.append(branch_bonds)
        fields.append(branch_fields)
        cycle_fluxes.append(spec.lattice.cycle_fluxes(branch_bonds))
        plaquette_fluxes.append(spec.lattice.plaquette_fluxes(branch_bonds))
        class_ids.append(spec.lattice.class_id(branch_bonds))

    raw_bits = tuple(0 if sign == 1 else 1 for row in steps for sign in row)
    return BranchRecord(
        measured_signs_by_step=steps,
        effective_bond_signs_by_step=tuple(effective_bonds),
        field_signs_by_step=tuple(fields),
        cycle_fluxes_by_step=tuple(cycle_fluxes),
        plaquette_fluxes_by_step=tuple(plaquette_fluxes),
        class_ids_by_step=tuple(class_ids),
        raw_bits=raw_bits,
    )


def record_from_bits(
    spec: TFIMSpec,
    bits: Sequence[int],
    steps: int = 1,
    architecture: str = "coherent_reuse",
) -> BranchRecord:
    bits = tuple(int(bit) for bit in bits)
    rows = 1 if architecture == "coherent_reuse" else int(steps)
    expected = rows * spec.num_terms
    if len(bits) != expected:
        raise ValueError(f"Expected {expected} record bits, received {len(bits)}.")
    signs = bits_to_signs(bits)
    sign_rows = tuple(
        signs[start : start + spec.num_terms]
        for start in range(0, len(signs), spec.num_terms)
    )
    return record_from_signs(spec, sign_rows)

