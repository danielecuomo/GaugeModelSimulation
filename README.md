# Gauge Model Simulation

![Wolfram Language](https://img.shields.io/badge/language-Wolfram%20Language-orange?logo=wolfram)
![Mathematica](https://img.shields.io/badge/Mathematica-14.2-red?logo=wolfram)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Status](https://img.shields.io/badge/status-In%20Progress-yellow)

This repository contains Wolfram Language code and notebooks for simulating **lattice gauge models** in imaginary time.  
The project provides tools for analyzing observables, correlation functions, and dynamics of gauge-invariant systems.

---

## 🔹 Repository structure
- `Notebooks/` → interactive `.nb` notebooks for development and exploration.  
- `Code/` → `.wl` code exports for GitHub-friendly reading.  
- `Figures/` → generated plots and diagrams.  
- `Data/` → optional datasets.  

---

## 🔹 Example (Mathematica)

```mathematica
(* Example Hamiltonian evolution *)
L = 6;  (* lattice size *)
g = 1.0; J = 1.0;

ψ0 = RandomVector[NormalDistribution[0, 1], {2^L}];
τ = 5;

H = (* define Hamiltonian here *);
Uτ = MatrixExp[-τ H];
ψτ = Normalize[Uτ.ψ0];
Eτ = Conjugate[ψτ].(H.ψτ) // Chop
