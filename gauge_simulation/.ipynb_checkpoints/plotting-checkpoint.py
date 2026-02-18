# gauge_simulation/plotting.py

import matplotlib.pyplot as plt
import numpy as np

# A dictionary to map user-friendly names to internal data keys and plot details
PLOT_DATA_MAP = {
    'thermal_energy': {
        'key': 'thermal_avgs', 
        'y_label': r'Thermal average ($\langle H \rangle$)', 
        'x_label': r'$\beta$',
        'title': 'Thermal Energy',
        'exact_color': 'skyblue', 
        'sim_color': 'lightgreen',
        'sim_plot_type': 'plot_scatter', # Scatter + dashed fit
        'x_key': 'betas'
    },
    'free_energy': {
        'key': 'free_energy', 
        'y_label': r'Free Energy ($\tilde{F}(\beta)$)', 
        'x_label': r'$\beta$',
        'title': 'Constant-Free Energy',
        'exact_color': 'skyblue', 
        'sim_color': 'lightgreen',
        'sim_plot_type': 'plot_scatter_line', # Scatter + solid line
        'x_key': 'free_energy_betas'
    },
    'entropy': {
        'key': 'entropy', 
        'y_label': r'Entropy ($\tilde{S}(\beta)$)', 
        'x_label': r'$\beta$',
        'title': 'Constant-Free Entropy',
        'exact_color': 'skyblue', 
        'sim_color': 'lightgreen',
        'sim_plot_type': 'plot_scatter_line',
        'x_key': 'entropy_betas'
    },
    'variance': {
        'key': 'variance', 
        'y_label': r'Variance ($\text{Var}(H)$)', 
        'x_label': r'$\beta$',
        'title': 'Hamiltonian Variance',
        'exact_color': 'skyblue', 
        'sim_color': 'lightgreen',
        'sim_plot_type': 'plot_scatter_line',
        'x_key': 'betas'
    }
}

def plot_simulation_results(
    exact_results: dict, 
    sim_results: dict, 
    quantities: list[str] | None = None,
    min_eigenval: float | None = None
) -> str:
    """
    Generates a figure with flexible subplots for specified thermal quantities.

    Parameters
    ----------
    exact_results : dict
        Dictionary of exact quantities from GaugeModel.get_exact_thermal_quantities.
    sim_results : dict
        Dictionary of simulation quantities from GaugeModel.run_simulation.
    quantities : list[str] | None
        List of quantities to plot. Options: 
        'thermal_energy', 'free_energy', 'entropy', 'variance'.
        If None, plots only 'thermal_energy'.
    min_eigenval : float | None
        The exact ground state energy for the thermal energy plot.
        
    Returns
    -------
    str
        The filename of the saved plot image.
    """
    valid_quantities = list(PLOT_DATA_MAP.keys())

    # Make sure at least thermal_energy is plot is plotting is called but no quantities indicated
    if quantities is None or not quantities:
        quantities = ['thermal_energy']
    
    # Remove duplicates 
    quantities = list(dict.fromkeys(quantities))

    invalid_quantities = []
    for quantity in quantities:
        if quantity not in valid_quantities:
            invalid_quantities.append(quantity)
            quantities.remove(quantity)

    if invalid_quantities:
        raise ValueError(
            f"Invalid quantity/quantities: {invalid_quantities}. \n"
            f"Valid options are: {valid_quantities}"
        )
    
    num_plots = len(quantities)
    
    # Determine the grid layout based on the number of plots
    if num_plots == 1:
        layout = (1, 1)
    elif num_plots == 2:
        layout = (1, 2)
    else:
        # 3 or 4 plots, use a 2x2 grid
        layout = (2, 2)
    rows, cols = layout
    figsize = (4*cols, 3*rows)
    fig = plt.figure(figsize=figsize)

    for i, q_name in enumerate(quantities):

        plot_def = PLOT_DATA_MAP[q_name]
        data_key = plot_def['key']
        x_key = plot_def['x_key']

        # Get the subplot index (1-based)
        ax = fig.add_subplot(rows, cols, i + 1)
        
        # --- Plot Exact Results ---
        exact_x = exact_results[x_key]
        exact_y = exact_results[data_key]

        if len(exact_y) > 0:
            # Scatter plot for exact data
            ax.scatter(exact_x, exact_y, color=plot_def['exact_color'], 
                    #    label=f'Exact $\\langle H \\rangle$' if q_name=='thermal_energy' else f'Exact $\\tilde{{ {q_name[0].upper()} }}$')
                        label = "exact" )
            
            # Plot line for non-thermal energy quantities
            if plot_def['sim_plot_type'] == 'plot_scatter_line':
                ax.plot(exact_x, exact_y, color=plot_def['exact_color'], linestyle='-')

            # For thermal energy, plot the fitted line (dashed)
            if q_name == 'thermal_energy':
                ax.plot(exact_x, exact_results['fitted_exact_avgs'], color=plot_def['sim_color'], linestyle='--')
            
        # --- Plot Simulation Results ---
        sim_x = sim_results[x_key]
        sim_y = sim_results[data_key]
        
        if len(sim_y) > 0:
            # Scatter plot for simulation data
            ax.scatter(sim_x, sim_y, color=plot_def['sim_color'], 
                    #    label=f'Sim $\\langle H \\rangle$' if q_name=='thermal_energy' else f'Sim $\\tilde{{ {q_name[0].upper()} }}$')
                        label = "sim")
            
            # Plot line for non-thermal energy quantities
            if plot_def['sim_plot_type'] == 'plot_scatter_line':
                ax.plot(sim_x, sim_y, color=plot_def['sim_color'], linestyle='-')
            
            # For thermal energy, plot the fitted line (dashed)
            if q_name == 'thermal_energy':
                ax.plot(sim_x, sim_results['fitted_estimated_avgs'], color=plot_def['sim_color'], linestyle='--')


        # --- Add Ground State Line and Axes Labels for Thermal Energy ---
        if q_name == 'thermal_energy' and min_eigenval is not None:
            ax.axhline(y=min_eigenval, linestyle=':', linewidth=2, color='gray', label='Ground State')
            # Adjust y-limit
            current_ylim = ax.get_ylim()
            if min_eigenval < current_ylim[0]:
                ax.set_ylim(min_eigenval * 1.1, current_ylim[1]) # Ensure min_eigenval is visible
            
        ax.set_xlabel(plot_def['x_label'])
        ax.set_ylabel(plot_def['y_label'])
        # ax.set_title(plot_def['title'])
        ax.grid(True, alpha=0.5)
        ax.legend()

    plt.tight_layout()
    plt.show()
    
    return 'gauge_sim_plot'