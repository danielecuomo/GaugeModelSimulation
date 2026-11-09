# Quantum simulation in imaginary time for gauge-invariant models

[![DOI](https://img.shields.io/badge/DOI-10.1088%2F1402--4896%2Fae4dce-blue)](https://doi.org/10.1088/1402-4896/ae4dce)
[![Qiskit](https://img.shields.io/badge/Qiskit-Mentorship-6929C4)](https://github.com/qiskit-advocate/qamp-2025/issues/5)
![Python](https://img.shields.io/badge/Python-3.12%2B-3776AB?logo=python&logoColor=white)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
[![Qiskit Ecosystem](https://qisk.it/e-0747c243)](https://qisk.it/e)


This repository contains Jupyter Notebooks for
simulating gauge-invariant models in imaginary time and computing thermal
averages.  It has two integrated paths:

- **Legacy TFIM/XY:** the original object-oriented API and published
  notebooks.  These remain stable for reproducibility.
- **Ancilla-resolved TFIM:** an explicit extension which retains every
  ancilla record, classifies gauge fluxes, and reconstructs
  class-conditioned thermodynamic quantities.

<p align="center">
  <img src="https://github.com/user-attachments/assets/f941d69b-26e8-4fc3-8688-a06b39ce46ed" width="350" alt="Description">
  <img src="https://github.com/user-attachments/assets/e4a8e1ce-cea8-42ba-a1a2-1c562cac03f4" width="350" alt="Description">
</p>

### Ground–State Evaluation

The package allows you to obtain both the **approximate** ground–state energy
(from the circuit-based imaginary–time method) and the **exact** ground–state
energy (from diagonalization of the stored Hamiltonian operator). 
Check `gauge_simulation_run.py` for other accessible quantities.

```python
from gauge_simulation import models

H = models.build_hamiltonian(model = "tfim", lattice_size=[1,3])

exact_groundstate = H.exact_ground_state()

approx_groundstate = H.approximate_ground_state()

print(f"Exact:  {exact_groundstate:.4f}")
print(f"Approx: {approx_groundstate:.4f}")
```

```text
Exact:  -2.7441
Approx: -2.5255
```

The bond parameters are currently set to π/4.

## Installation

Install the required dependencies and the package itself:

```bash
pip install -r requirements.txt
pip install -e .
```

The second command installs the package in editable mode, allowing the example notebooks to import `gauge_simulation` directly.


## Legacy notebooks

The original notebooks remain in [`notebooks/`](notebooks), unchanged:

- [`Notebook_Implementation_v1_1.ipynb`](notebooks/Notebook_Implementation_v1_1.ipynb)
  is the object-oriented TFIM/XY workflow;
- [`Notebook_Theory_v1.ipynb`](notebooks/Notebook_Theory_v1.ipynb) documents
  the theoretical foundations; and
- the TFIM and XY notebooks retain the original model-specific derivations.

Imports such as `from gauge_simulation import models, plotting` remain
supported.  New code may make this choice explicit with
`gauge_simulation.legacy`.

## Ancilla-resolved TFIM extension

The successor implementation is an explicit subpackage:
`gauge_simulation.ancilla_resolved`.  It reuses the original TFIM
ancilla–Pauli primitive and adds configurable circuit architectures,
gauge-flux classification, class-conditioned observables, quenched
thermodynamic reconstruction, and modular simulator or hardware execution.

```python
from gauge_simulation.ancilla_resolved import (
    ExecutionConfig,
    GaugeITEExperiment,
    ITEProtocol,
    TFIMSpec,
)

spec = TFIMSpec.rectangular(rows=2, cols=2)
protocol = ITEProtocol(
    beta=0.75,
    angle_mapping="exact_kraus",
    architecture="coherent_reuse",
    input_mode="purification",
)
execution = ExecutionConfig(provider="aer", shots=4096, seed=2026)
result = GaugeITEExperiment(spec, protocol, execution=execution).run(
    "class_probabilities"
)
print(result.classes)
```

The ancilla-resolved notebooks are in
[`notebooks/ancilla_resolved/`](notebooks/ancilla_resolved).  They begin with
[`00_from_legacy_to_ancilla_resolved.ipynb`](notebooks/ancilla_resolved/00_from_legacy_to_ancilla_resolved.ipynb),
which compares the two APIs, then cover energy estimation, circuit
architectures, record classification, class thermodynamics, quenched
reconstruction, and backend selection.

The corresponding executable scripts remain in [`examples/`](examples).
The command-line entry point is also available after installation:

```bash
gauge-ite run --rows 2 --cols 2 --beta 0.75 --observable energy
```

For the complete ancilla-resolved API, thermodynamics, and hardware workflow,
see [the extension guide](docs/ancilla_resolved.md).
