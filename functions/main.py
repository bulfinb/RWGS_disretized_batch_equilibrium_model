import numpy as np
from functions import build_arrays
from functions import build_mu_O_functions
from functions.discretized_reduction import compute_reduction_step
from functions.discretized_oxidation import compute_oxidation_step

def simulate_cycle(material="CeO2", first_cycle=True, delta_x_0=np.zeros(100),
                   T=1073, n_CO2=1.0, n_H2=1.01, n_oxide=20,
                   x_H2O_0=0.005, x_CO2_0=0.998,
                   oxide_mesh=100, gas_mesh=100):
    """
    Simulate a single chemical looping cycle (reduction + oxidation) in a 1D reactor.

    Args:
        material (str): Material name (options, "CeO2", "LSF", "CeO2_simple", "CeZr05", "CeZr15").
        first_cycle (bool): True if this is the first cycle (used to set starting codition).
        delta_x_0 (ndarray): Initial delta profile if not first cycle.
        T (float): Operating temperature in Kelvin.
        n_CO2 (float): Moles of CO2 used in oxidation.
        n_H2 (float): Moles of H2 used in reduction.
        n_oxide (float): Moles of solid oxide in the reactor.
        x_H2O_0 (float): Initial mole fraction of H2O in feed gas.
        x_CO2_0 (float): Initial mole fraction of CO2 in feed gas.
        oxide_mesh (int): Number of discretization points in solid.
        gas_mesh (int): Number of discretization points in gas (time-like).

    Returns:
        tuple:
            - delta_t_x_red (ndarray): Delta solution matrix after reduction.
            - x_H2O_t_x_red (ndarray): H2O solution matrix after reduction.
            - delta_t_x_ox (ndarray): Delta solution matrix after oxidation.
            - x_CO2_t_x_ox (ndarray): CO2 solution matrix after oxidation.
    """
    # Mass balance scaling factors
    mbf_red = (n_oxide / n_H2) * (gas_mesh / oxide_mesh)
    mbf_ox = (n_oxide / n_CO2) * (gas_mesh / oxide_mesh)

    # Set delta_max and delta range
    delta_max, delta_range = build_arrays.set_delta_max_and_delta_range(material)

    # Generate mu_O(delta) function for the oxide and its inverse
    mu_O_delta_func, mu_O_delta_func_inv = build_mu_O_functions.set_mu_O_delta_functions(material, T, delta_range)

    # Generate X (gas composition) range
    X_range = build_arrays.X_array()

    # Generate mu_O(X) functions for gas streams
    mu_O_CO2_func = build_mu_O_functions.set_mu_O_CO2_function(T, X_range)
    mu_O_H2O_func = build_mu_O_functions.set_mu_O_H2O_function(T, X_range)

    # Determine delta bounds from gas stream equilibria
    max_oxidising_mu_O = mu_O_CO2_func(x_CO2_0)
    delta_min = mu_O_delta_func_inv(max_oxidising_mu_O)
    min_reducing_mu_O = mu_O_H2O_func(x_H2O_0)
    delta_max_reduction = mu_O_delta_func_inv(min_reducing_mu_O)

    if delta_max < delta_max_reduction:
        print("Warning: Gas thermodynamics outside material range (H2 too reducing).")
    else:
        delta_max = delta_max_reduction

    # Initialize simulation arrays
    delta_t_x_red = np.zeros((oxide_mesh, gas_mesh)).T
    x_H2O_t_x_red = np.zeros((oxide_mesh, gas_mesh)).T
    delta_t_x_ox = np.zeros((oxide_mesh, gas_mesh)).T
    x_CO2_t_x_ox = np.zeros((oxide_mesh, gas_mesh)).T

    # Initial condition for reduction
    if first_cycle:
        delta_t_x_red[0] = delta_min
    else:
        delta_t_x_red[0] = delta_x_0

    # Create differential change in oxygen mass balance arrays for batch equilibrium solver
    d_delta_red, d_X_red = build_arrays.mass_balance_arrays(
        delta_max=delta_max, n_gas=n_H2,
        n_oxide=n_oxide, gas_mesh=gas_mesh, oxide_mesh=oxide_mesh
    )
    d_delta_ox, d_X_ox = build_arrays.mass_balance_arrays(
        delta_max=delta_max, n_gas=n_CO2,
        n_oxide=n_oxide, gas_mesh=gas_mesh, oxide_mesh=oxide_mesh
    )

    # Run reduction simulation
    delta_t_x_red, x_H2O_t_x_red = compute_reduction_step(
        mu_O_delta_func, mu_O_H2O_func, x_H2O_0,
        delta_t_x_red, x_H2O_t_x_red, d_delta_red, d_X_red,
        mbf_red, gas_mesh=gas_mesh, oxide_mesh=oxide_mesh
    )

    # Use reduction output as initial condition for oxidation
    delta_t_x_ox[0] = delta_t_x_red[-1]

    # Run oxidation simulation
    delta_t_x_ox, x_CO2_t_x_ox = compute_oxidation_step(
        mu_O_delta_func, mu_O_CO2_func, x_CO2_0, delta_min,
        delta_t_x_ox, x_CO2_t_x_ox, d_delta_ox, d_X_ox,
        mbf_ox, gas_mesh=gas_mesh, oxide_mesh=oxide_mesh
    )

    return delta_t_x_red, x_H2O_t_x_red, delta_t_x_ox, x_CO2_t_x_ox


def cycle_until_balanced(max_cycles = 20, O_balance_tolerance = 0.001, material="CeO2",
                   T=1073, n_CO2=1.0, n_H2=1.01, n_oxide=20,
                   x_H2O_0=0.005, x_CO2_0=0.998,
                   oxide_mesh=100, gas_mesh=100):
    """
    Run chemical looping cycle until a mass balance is reached for the gas phase streams.

    Args:
        material (str): Material name (options, "CeO2", "LSF", "CeO2_simple", "CeZr05", "CeZr15").
        delta_x_0 (ndarray): Initial delta profile if not first cycle.
        T (float): Operating temperature in Kelvin.
        n_CO2 (float): Moles of CO2 used in oxidation.
        n_H2 (float): Moles of H2 used in reduction.
        n_oxide (float): Moles of solid oxide in the reactor.
        x_H2O_0 (float): Initial mole fraction of H2O in feed gas.
        x_CO2_0 (float): Initial mole fraction of CO2 in feed gas.
        oxide_mesh (int): Number of discretization points in solid.
        gas_mesh (int): Number of discretization points in gas (time-like).

    Returns:
        tuple:
            - delta_t_x_red (ndarray): Delta solution matrix after reduction.
            - x_H2O_t_x_red (ndarray): H2O solution matrix after reduction.
            - delta_t_x_ox (ndarray): Delta solution matrix after oxidation.
            - x_CO2_t_x_ox (ndarray): CO2 solution matrix after oxidation.
    """
    # Run the first cycle with some excess H2 as we start fully oxidised
    delta_t_x_red, x_H2O_t_x_red, delta_t_x_ox, x_CO2_t_x_ox = simulate_cycle(material=material,
                                                                              first_cycle=True,
                                                                              T=T, n_CO2 = n_CO2, n_oxide=n_oxide,
                                                                              n_H2=n_H2,
                                                                              x_CO2_0=x_CO2_0, x_H2O_0=x_H2O_0,
                                                                              oxide_mesh=oxide_mesh, gas_mesh=gas_mesh
                                                                              )

    # Set Starting condition for the next cycle
    delta_x_0 = delta_t_x_ox[-1]
    # Set O_bal_gas low enough that at least one more cycle is run.
    O_bal_gas = 0.9
    cycles = 1
    while (abs(1.0 - O_bal_gas) > O_balance_tolerance and cycles < max_cycles):
        delta_t_x_red, x_H2O_t_x_red, delta_t_x_ox, x_CO2_t_x_ox = simulate_cycle(material=material,
                                                                                  first_cycle=False,
                                                                                  delta_x_0 = delta_x_0,
                                                                                  T=T, n_CO2 = n_CO2, n_oxide=n_oxide,
                                                                                  n_H2=n_H2,
                                                                                  x_CO2_0=x_CO2_0, x_H2O_0=x_H2O_0,
                                                                                  oxide_mesh=oxide_mesh, gas_mesh=gas_mesh
                                                                                 )
        # Set Starting condition for the next cycle
        delta_x_0 = delta_t_x_ox[-1]
        # Calculate the mass balance
        n_H2O = (x_H2O_t_x_red.T[-1].mean() - x_H2O_0) * n_H2
        n_CO = (x_CO2_0 - x_CO2_t_x_ox.T[0].mean()) * n_CO2
        O_bal_gas = n_CO / n_H2O
        cycles += 1
    return delta_t_x_red, x_H2O_t_x_red, delta_t_x_ox, x_CO2_t_x_ox, cycles
