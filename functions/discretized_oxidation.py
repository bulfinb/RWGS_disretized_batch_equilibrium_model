import numpy as np
from functions import batch_equilibrium_solver

# Tolerances used to determine if the reaction should be terminated
delta_tolerance = 0.0001  # Stop if delta is close enough to min
x_gas_tolerance = 0.001   # Stop if gas is almost fully converted

def oxidation_t0_x0_bc(mu_O_delta_func, mu_O_CO2_func, x_CO2_0,
                       delta_min_oxidation, delta_t0_x0,
                       d_delta, d_X, mbf_ox):
    """
    Compute oxidation equilibrium at the initial grid node (t=0, x=0).

    Args:
        mu_O_delta_func (callable): Function for mu_O in the solid.
        mu_O_CO2_func (callable): Function for mu_O in the CO2 gas stream.
        x_CO2_0 (float): Initial mole fraction of CO2.
        delta_min_oxidation (float): Minimum allowed delta for oxidation.
        delta_t0_x0 (float): Initial delta at t=0, x=0.
        d_delta (ndarray): Array of delta differentials.
        d_X (ndarray): Array of gas composition differentials.
        mbf_ox (float): Mass balance factor for oxidation.

    Returns:
        tuple: Updated delta and CO2 mole fraction at (t=1, x=0).
    """
    if delta_t0_x0 < delta_min_oxidation - delta_tolerance:
        return delta_t0_x0, x_CO2_0

    d_d = batch_equilibrium_solver.oxidation(mu_O_delta_func, mu_O_CO2_func,
                                             delta_t0_x0, x_CO2_0, d_delta, d_X)
    delta_t1_x0 = delta_t0_x0 - d_d[0]
    x_CO2_t1_x0 = x_CO2_0 - d_d[0] * mbf_ox
    return delta_t1_x0, x_CO2_t1_x0

def oxidation_x0_bc(mu_O_delta_func, mu_O_CO2_func, x_CO2_0, delta_min_oxidation,
                    delta_t_x0, x_CO2_t_x0, d_delta, d_X, mbf_ox, gas_mesh=100):
    """
    Compute oxidation progression along time (t) at x = 0 element.

    Args:
        mu_O_delta_func (callable): Function for mu_O in the solid.
        mu_O_CO2_func (callable): Function for mu_O in the CO2 gas stream.
        x_CO2_0 (float): Initial mole fraction of CO2.
        delta_min_oxidation (float): Minimum delta to permit oxidation.
        delta_t_x0 (ndarray): Array of delta values at x=0 over time.
        x_CO2_t_x0 (ndarray): CO2 mole fraction array at x=0 over time.
        d_delta (ndarray): Delta differentials.
        d_X (ndarray): Gas composition differentials.
        mbf_ox (float): Mass balance factor for oxidation.
        gas_mesh (int): Number of gas-phase elements (like time steps).

    Returns:
        tuple: Updated arrays for delta and CO2 at x=0.
    """
    for t_step in range(1, gas_mesh):
        # Stop if oxide element is fully oxidized and fill the solution forward in time
        if delta_t_x0[t_step - 1] < delta_min_oxidation + delta_tolerance:
            delta_t_x0[t_step - 1:] = delta_t_x0[t_step - 1]
            x_CO2_t_x0[t_step - 1:] = x_CO2_0
            break

        d_d = batch_equilibrium_solver.oxidation(mu_O_delta_func, mu_O_CO2_func,
                                                 delta_t_x0[t_step - 1], x_CO2_0, d_delta, d_X)
        delta_t_x0[t_step] = delta_t_x0[t_step - 1] - d_d[0]
        x_CO2_t_x0[t_step] = x_CO2_0 - d_d[0] * mbf_ox
    return delta_t_x0, x_CO2_t_x0

def oxidation_t0_bc(mu_O_delta_func, mu_O_CO2_func, delta_min_oxidation,
                    delta_t0_x, x_CO2_t0_x, d_delta, d_X, mbf_ox, oxide_mesh=100):
    """
    Compute oxidation progression along space (x) at t = 0.

    Args:
        mu_O_delta_func (callable): Function for mu_O in the solid.
        mu_O_CO2_func (callable): Function for mu_O in the CO2 gas stream.
        delta_min_oxidation (float): Minimum delta to allow oxidation.
        delta_t0_x (ndarray): Initial delta values at t=0 along x.
        x_CO2_t0_x (ndarray): CO2 mole fractions at t=0 along x.
        d_delta (ndarray): Delta differential array.
        d_X (ndarray): Gas differential array.
        mbf_ox (float): Mass balance factor.
        oxide_mesh (int): Number of solid phase spatial steps.

    Returns:
        tuple: Updated delta and CO2 profiles at t = 0.
    """
    for x_step in range(1, oxide_mesh):
        # Skip over fully oxidized regions
        if delta_t0_x[x_step] < delta_min_oxidation - delta_tolerance:
            x_CO2_t0_x[x_step] = x_CO2_t0_x[x_step - 1]
            continue

        d_d = batch_equilibrium_solver.oxidation(mu_O_delta_func, mu_O_CO2_func,
                                                 delta_t0_x[x_step], x_CO2_t0_x[x_step - 1],
                                                 d_delta, d_X)
        delta_t0_x[x_step] -= d_d[0]
        x_CO2_t0_x[x_step] = x_CO2_t0_x[x_step - 1] - d_d[0] * mbf_ox
    return delta_t0_x, x_CO2_t0_x

def oxidation_x_t(mu_O_delta_func, mu_O_CO2_func, delta_min_oxidation,
                  delta_t_x, x_CO2_t_x, d_delta, d_X, mbf_ox,
                  gas_mesh=100, oxide_mesh=100):
    """
    Compute oxidation across the full t-x grid.

    Args:
        mu_O_delta_func (callable): mu_O in solid.
        mu_O_CO2_func (callable): mu_O in gas.
        delta_min_oxidation (float): Minimum delta threshold.
        delta_t_x (ndarray): Delta grid.
        x_CO2_t_x (ndarray): CO2 mole fraction grid.
        d_delta (ndarray): Delta differential array.
        d_X (ndarray): Gas fraction differential array.
        mbf_ox (float): Mass balance factor.
        gas_mesh (int): Number of gas-phase elements (like time steps).
        oxide_mesh (int): Number of oxide elements (spatial x steps).

    Returns:
        tuple: Updated delta and gas composition grids.
    """
    for x_step in range(1, oxide_mesh):
        for t_step in range(1, gas_mesh):
            # Stop if oxide element is fully oxidized and fill the solution forward in time
            if delta_t_x[t_step - 1, x_step] < delta_min_oxidation + delta_tolerance:
                delta_t_x[t_step:, x_step] = delta_t_x[t_step - 1, x_step]
                x_CO2_t_x[t_step:, x_step] = x_CO2_t_x[t_step, x_step - 1]
                break

            d_d = batch_equilibrium_solver.oxidation(mu_O_delta_func, mu_O_CO2_func,
                                                     delta_t_x[t_step - 1, x_step],
                                                     x_CO2_t_x[t_step, x_step - 1],
                                                     d_delta, d_X)
            delta_t_x[t_step, x_step] = delta_t_x[t_step - 1, x_step] - d_d[0]
            x_CO2_t_x[t_step, x_step] = x_CO2_t_x[t_step, x_step - 1] - d_d[0] * mbf_ox
    return delta_t_x, x_CO2_t_x

def compute_oxidation_step(mu_O_delta_func, mu_O_CO2_func, x_CO2_0, delta_min_oxidation,
                           delta_t_x, x_CO2_t_x, d_delta, d_X, mbf_ox,
                           gas_mesh=100, oxide_mesh=100):
    """
    Perform full oxidation step across the reactor grid.

    This includes setting boundary conditions and propagating oxidation across the full t-x mesh.

    Args:
        mu_O_delta_func (callable): mu_O(delta) for solid.
        mu_O_CO2_func (callable): mu_O(X) for CO2 gas.
        x_CO2_0 (float): Initial CO2 mole fraction.
        delta_min_oxidation (float): Minimum delta threshold.
        delta_t_x (ndarray): Delta grid.
        x_CO2_t_x (ndarray): Gas mole fraction grid.
        d_delta (ndarray): Delta differential array.
        d_X (ndarray): Gas differential array.
        mbf_ox (float): Mass balance factor.
        gas_mesh (int): Number of gas elements.
        oxide_mesh (int): Number of solid elements.

    Returns:
        tuple: Updated delta and CO2 fraction grids after oxidation.
    """
    # Reverse solid profile to simulate counter-current oxidation
    delta_t_x[0] = np.flip(delta_t_x[0])

    # Apply boundary conditions and solve over the mesh
    delta_t_x[0, 0], x_CO2_t_x[0, 0] = oxidation_t0_x0_bc(
        mu_O_delta_func, mu_O_CO2_func, x_CO2_0,
        delta_min_oxidation, delta_t_x[0, 0], d_delta, d_X, mbf_ox
    )
    delta_t_x[:, 0], x_CO2_t_x[:, 0] = oxidation_x0_bc(
        mu_O_delta_func, mu_O_CO2_func, x_CO2_0, delta_min_oxidation,
        delta_t_x[:, 0], x_CO2_t_x[:, 0], d_delta, d_X, mbf_ox, gas_mesh=gas_mesh
    )
    delta_t_x[0], x_CO2_t_x[0] = oxidation_t0_bc(
        mu_O_delta_func, mu_O_CO2_func, delta_min_oxidation,
        delta_t_x[0], x_CO2_t_x[0], d_delta, d_X, mbf_ox, oxide_mesh=oxide_mesh
    )
    delta_t_x, x_CO2_t_x = oxidation_x_t(
        mu_O_delta_func, mu_O_CO2_func, delta_min_oxidation,
        delta_t_x, x_CO2_t_x, d_delta, d_X, mbf_ox,
        gas_mesh=gas_mesh, oxide_mesh=oxide_mesh
    )

    # Flip back to original orientation
    return np.flip(delta_t_x, axis=1), np.flip(x_CO2_t_x, axis=1)

