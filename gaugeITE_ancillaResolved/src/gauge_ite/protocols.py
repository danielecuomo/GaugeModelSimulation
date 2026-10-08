"""Protocol choices with explicit multi-step interpretation."""

from __future__ import annotations

from dataclasses import dataclass


ARCHITECTURES = {
    "coherent_reuse",
    "fresh_trajectory",
    "measure_reset_trajectory",
}
ANGLE_MAPPINGS = {"linear", "exact_kraus"}
INPUT_MODES = {"purification", "basis", "zero"}


@dataclass(frozen=True)
class ITEProtocol:
    beta: float
    steps: int = 1
    angle_mapping: str = "exact_kraus"
    architecture: str = "coherent_reuse"
    input_mode: str = "purification"
    basis_bits: tuple[int, ...] | None = None

    def __post_init__(self) -> None:
        angle_mapping = self.angle_mapping.lower()
        architecture = self.architecture.lower()
        input_mode = self.input_mode.lower()
        if self.beta < 0:
            raise ValueError("beta must be non-negative.")
        if self.steps <= 0:
            raise ValueError("steps must be a positive integer.")
        if angle_mapping not in ANGLE_MAPPINGS:
            raise ValueError(f"angle_mapping must be one of {sorted(ANGLE_MAPPINGS)}.")
        if architecture not in ARCHITECTURES:
            raise ValueError(f"architecture must be one of {sorted(ARCHITECTURES)}.")
        if input_mode not in INPUT_MODES:
            raise ValueError(f"input_mode must be one of {sorted(INPUT_MODES)}.")
        if input_mode == "basis" and self.basis_bits is None:
            raise ValueError("basis_bits are required when input_mode='basis'.")
        if self.basis_bits is not None and any(bit not in (0, 1) for bit in self.basis_bits):
            raise ValueError("basis_bits must contain only 0 and 1.")
        object.__setattr__(self, "angle_mapping", angle_mapping)
        object.__setattr__(self, "architecture", architecture)
        object.__setattr__(self, "input_mode", input_mode)

    @property
    def interpretation(self) -> str:
        if self.architecture == "coherent_reuse" and self.steps == 1:
            return "one-step sign-configured product formula"
        if self.architecture == "coherent_reuse":
            return "operational final-record branch; static H_s is not assumed"
        return "measured time-ordered sign trajectory"

