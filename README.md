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
import numpy as np
from gauge_simulation import models
from gauge_simulation import plotting

# User input and changes: 
model = "tfim" # or "xy"
lattice_size = [2, 2]
beta_vals = np.linspace(0.0, 1.4, 11)
plot_circuit = True # or False

# Build the Hamiltonian Model (set model type, and lattice_size)
H = models.build_hamiltonian(model, lattice_size) # H is an instance of GaugeModel


# Get Exact Ground State
exact_groundstate = H.exact_ground_state()

# Get Exact Thermodynamic Quantities
exact_results = H.get_exact_thermal_quantities(beta_vals)

# Run the full simulation across the range of $\beta$ values
sim_results = H.run_simulation(beta_vals)

# Get Approximated Ground State
approx_groundstate = sim_results['thermal_avgs'][-1]

print(f"Exact Ground State Energy (E_min):      {exact_groundstate:.4f}")
print(f"Exact <H>:                              {exact_results['thermal_avgs'][-1]:.4f}")
print(f"Simulated <H>:                          {sim_results['thermal_avgs'][-1]:.4f} \n"
      f"Simulation variance:                        {sim_results['variance'][-1]:.4f}"
      )

# Circuit plotting
if plot_circuit:
    #circuit = H.plot_circuit()
    H.plot_circuit()


# Plotting all four thermodynamic quantities in a 2x2 grid
# quantities_to_plot = ['thermal_energy', 'free_energy', 'entropy', 'variance']
quantities_to_plot = ['thermal_energy', 'free_energy', 'entropy', 'variance']

plot_filename = plotting.plot_simulation_results(
    exact_results=exact_results, 
    sim_results=sim_results,
    quantities=quantities_to_plot,
    min_eigenval=exact_groundstate,
)

```

```text
Exact Ground State Energy (E_min):      -4.1047
Exact <H>:                              -3.8079
Simulated <H>:                          -3.7006 
Simulation variance:                        0.7475
```

The bond parameters are currently set to π/4.
