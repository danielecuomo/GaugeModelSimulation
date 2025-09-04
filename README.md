# Imaginary-Time Simulations – TFIM & XY Model

This repository contains Wolfram Language code and notebooks for studying **imaginary-time dynamics** in the **Transverse-Field Ising Model (TFIM)** and the **XY model**.  
The scripts compute eigenvalues, correlation functions, and physical observables from lattice Hamiltonians.

---

## 🔹 Repository structure
- `Notebooks/` → original `.nb` notebooks (interactive exploration and development).  
- `Code/` → `.wl` code exports, directly readable on GitHub.  
- `Data/` → optional saved numerical data.  
- `Figures/` → plots and generated images.  

---

## 🔹 Example Mathematica code

```mathematica
(* Parameters *)
L = 6;  (* chain length *)
g = 1.0;  (* transverse field *)
J = 1.0;  (* coupling *)

(* Pauli matrices *)
σx = PauliMatrix[1];
σz = PauliMatrix[3];

(* TFIM Hamiltonian *)
Hising[L_, g_, J_] := 
  -J Sum[
      KroneckerProduct[
        IdentityMatrix[2, SparseArray]^(i - 1),
        σz,
        σz,
        IdentityMatrix[2, SparseArray]^(L - i - 1)
      ]
    , {i, 1, L - 1}] 
  -g Sum[
      KroneckerProduct[
        IdentityMatrix[2, SparseArray]^(i - 1),
        σx,
        IdentityMatrix[2, SparseArray]^(L - i)
      ]
    , {i, 1, L}];

(* Imaginary-time evolution *)
ψ0 = RandomVector[NormalDistribution[0, 1], {2^L}];
τ = 5;
Uτ = MatrixExp[-τ Hising[L, g, J]];
ψτ = Normalize[Uτ.ψ0];

(* Expectation value of the energy *)
Eτ = Conjugate[ψτ].(Hising[L, g, J].ψτ) // Chop
