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

H = models.build_hamiltonian("tfim", [1,3])

apx = H.approx_ground_state()     # circuit-based approximation
ext = H.exact_ground_state()      # exact diagonalization

print("Approximated ground:", apx)
print("Exact ground:",        ext)
```

```text
Approximated ground: -2.522029838042481
Exact ground:        -2.744149144505002
```

The bond parameters are currently set to π/4.
