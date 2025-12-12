# Quantum simulation in imaginary time for gauge-invariant models

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17060897.svg)](https://doi.org/10.5281/zenodo.17060897)
![Mathematica](https://img.shields.io/badge/Mathematica-14.2-red?logo=wolfram)
![Python](https://img.shields.io/badge/Python-3.12%2B-3776AB?logo=python&logoColor=white)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Status](https://img.shields.io/badge/status-In%20Progress-yellow)


This repository contains Wolfram Language code and Jupyter Notebooks for simulating gauge-invariant models in imaginary time, computing thermal averages.

<p align="center">
  <img src="https://github.com/user-attachments/assets/e4a8e1ce-cea8-42ba-a1a2-1c562cac03f4"
       width="400"
       />
</p>

### Ground–State Evaluation

The package allows you to obtain both the **approximate** ground–state energy
(from the circuit-based imaginary–time method) and the **exact** ground–state
energy (from diagonalization of the stored Hamiltonian operator).

```python
from gauge_simulation import models

H = models.build_hamiltonian(model = "tfim", size = [2,2])

ext = H.exact_ground_state()

apx = H.run_simulation(1.4)

print(f"Ground State Energy:  {ext:.4f}")
print(f"Approximated:  {apx:.4f}")
```

```text
Ground State Energy:  -4.1047
Approximated:  -3.7006 
```

The bond parameters are currently set to π/4.
