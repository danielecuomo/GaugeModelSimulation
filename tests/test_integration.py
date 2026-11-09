"""Compatibility checks for the integrated legacy and ancilla-resolved APIs."""

from __future__ import annotations

import json
from pathlib import Path
import unittest

import numpy as np
from qiskit.quantum_info import Operator

from gauge_simulation import models
from gauge_simulation.ancilla_resolved import (
    ITEProtocol,
    TFIMSpec,
    clean_ground_energy,
)
from gauge_simulation.ancilla_resolved.circuits import build_tfim_circuit


class IntegrationTests(unittest.TestCase):
    """Exercise the public paths used by the legacy and new notebooks."""

    def test_uniform_tfim_exact_energy_is_compatible(self) -> None:
        legacy = models.build_hamiltonian("tfim", [1, 2], trotter_steps=1)
        extension = TFIMSpec.rectangular(rows=1, cols=2)
        self.assertAlmostEqual(
            legacy.exact_ground_state(), clean_ground_energy(extension)
        )

    def test_linear_coherent_circuit_matches_the_legacy_primitive(self) -> None:
        beta = 0.75
        legacy = models.build_hamiltonian("tfim", [1, 2], trotter_steps=1)
        legacy_circuit = legacy.qc.assign_parameters(
            {
                legacy.params["gamma"]: np.pi / 4,
                legacy.params["eta"]: np.pi / 4,
                legacy.params["beta"]: -beta,
            }
        )
        extension_circuit = build_tfim_circuit(
            TFIMSpec.rectangular(rows=1, cols=2),
            ITEProtocol(
                beta=beta,
                steps=1,
                angle_mapping="linear",
                architecture="coherent_reuse",
                input_mode="purification",
            ),
        ).circuit

        legacy_operator = Operator(legacy_circuit).data
        extension_operator = Operator(extension_circuit).data
        self.assertTrue(np.allclose(legacy_operator, extension_operator))

    def test_ancilla_resolved_notebooks_are_valid_and_compilable(self) -> None:
        notebook_dir = Path("notebooks/ancilla_resolved")
        notebooks = sorted(notebook_dir.glob("*.ipynb"))
        self.assertEqual(len(notebooks), 7)
        for path in notebooks:
            notebook = json.loads(path.read_text(encoding="utf-8"))
            def source_text(cell: dict[str, object]) -> str:
                source = cell["source"]
                return "".join(source) if isinstance(source, list) else str(source)

            code = "\n".join(
                source_text(cell)
                for cell in notebook["cells"]
                if cell["cell_type"] == "code"
            )
            compile(code, str(path), "exec")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
