import numpy as np
from gauge_simulation import models
from gauge_simulation import plotting


# User input and changes: 
model = "tfim" # or "xy"
lattice_size = [2, 2]
beta_vals = np.linspace(0.0, 1.4, 11)
plot_circuit = True #True or False

# Build the Hamiltonian Model (set model type, and lattice_size)
H = models.build_hamiltonian(model, lattice_size) # H is an instance of GaugeModel


# Get Exact Ground State
exact_groundstate = H.exact_ground_state()

# Get Approximate Ground State
approx_groundstate = H.approximate_ground_state()

# Get Exact Thermodynamic Quantities
exact_results = H.get_exact_thermal_quantities(beta_vals)

# Get Approximate (simulated) Thermodynamic Quantities
sim_results = H.get_approximate_thermal_quantities(beta_vals)

print(f"Exact Ground State Energy:      {exact_groundstate:.4f}")
print(f"Approx. Ground State Energy:    {approx_groundstate:.4f}")
print(f"Exact <H>:                      {exact_results['thermal_avgs'][-1]:.4f}")
print(f"Approx. <H>:                    {sim_results['thermal_avgs'][-1]:.4f} \n"
      f"Result variance:                 {sim_results['variance'][-1]:.4f}"
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
