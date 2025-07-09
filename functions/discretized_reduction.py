from functions import batch_equilibrium_solver

# Tolerances for stopping criteria
delta_tolerance = 0.0001  # Stop if delta approaches maximum allowable value
x_gas_tolerance = 0.001   # Stop if gas composition reaches near complete conversion

def reduction_t0_x0_bc(mu_O_delta_func, mu_O_H2O_func, x_H2O_0, delta_t0_x0,
                        d_delta, d_X, mbf_red):
    """
    Compute reduction at the initial grid node (t=0, x=0).

    Args:
        mu_O_delta_func (callable): Function for mu_O in the solid.
        mu_O_H2O_func (callable): Function for mu_O in the H2O gas stream.
        x_H2O_0 (float): Initial H2O mole fraction.
        delta_t0_x0 (float): Initial delta at t=0, x=0.
        d_delta (ndarray): Delta differential array.
        d_X (ndarray): Gas differential array.
        mbf_red (float): Mass balance factor for reduction.

    Returns:
        tuple: Updated delta and H2O mole fraction at (t=1, x=0).
    """
    d_d = batch_equilibrium_solver.reduction(mu_O_delta_func, mu_O_H2O_func,
                                             delta_t0_x0, x_H2O_0, d_delta, d_X)
    delta_t1_x0 = delta_t0_x0 + d_d[0]
    x_H2O_t1_x0 = x_H2O_0 + d_d[0] * mbf_red
    return delta_t1_x0, x_H2O_t1_x0

def reduction_t0_bc(mu_O_delta_func, mu_O_H2O_func, delta_t0_x, x_H2O_t0_x,
                    d_delta, d_X, mbf_red, oxide_mesh=100):
    """
    Compute reduction progression along space (x) at t = 0.

    Args:
        mu_O_delta_func (callable): Function for mu_O in solid.
        mu_O_H2O_func (callable): Function for mu_O in gas.
        delta_t0_x (ndarray): Delta values at t=0 along x.
        x_H2O_t0_x (ndarray): H2O mole fraction at t=0 along x.
        d_delta (ndarray): Delta differentials.
        d_X (ndarray): Gas composition differentials.
        mbf_red (float): Mass balance factor for reduction.
        oxide_mesh (int): Number of oxide elements (spatial x steps).

    Returns:
        tuple: Updated delta and H2O mole fraction arrays at t = 0.
    """
    for x_step in range(1, oxide_mesh):
        # Stop if gas is nearly fully oxidized
        if x_H2O_t0_x[x_step - 1] > 1 - x_gas_tolerance:
            x_H2O_t0_x[x_step - 1:] = x_H2O_t0_x[x_step - 1]
            break
        d_d = batch_equilibrium_solver.reduction(mu_O_delta_func, mu_O_H2O_func,
                                                 delta_t0_x[x_step], x_H2O_t0_x[x_step - 1],
                                                 d_delta, d_X)
        if len(d_d) < 1:
            x_H2O_t0_x[x_step] = x_H2O_t0_x[x_step - 1]
            continue
        delta_t0_x[x_step] += d_d[0]
        x_H2O_t0_x[x_step] = x_H2O_t0_x[x_step - 1] + d_d[0] * mbf_red
    return delta_t0_x, x_H2O_t0_x

def reduction_x0_bc(mu_O_delta_func, mu_O_H2O_func, x_H2O_0, delta_t_x0,
                    x_H2O_t_x0, d_delta, d_X, mbf_red, gas_mesh=100):
    """
    Compute reduction progression along time (t) at x = 0.

    Args:
        mu_O_delta_func (callable): Function for mu_O in the solid.
        mu_O_H2O_func (callable): Function for mu_O in the H2O gas stream.
        x_H2O_0 (float): Initial mole fraction of H2O.
        delta_t_x0 (ndarray): Delta values at x=0 over time.
        x_H2O_t_x0 (ndarray): H2O mole fraction values at x=0 over time.
        d_delta (ndarray): Delta differential array.
        d_X (ndarray): Gas composition differential array.
        mbf_red (float): Mass balance factor for reduction.
        gas_mesh (int): Number of gas-phase elements (like time steps).

    Returns:
        tuple: Updated delta and H2O mole fraction arrays at x = 0.
    """
    for t_step in range(1, gas_mesh):
        d_d = batch_equilibrium_solver.reduction(mu_O_delta_func, mu_O_H2O_func,
                                                 delta_t_x0[t_step - 1], x_H2O_0,
                                                 d_delta, d_X)
        delta_t_x0[t_step] = delta_t_x0[t_step - 1] + d_d[0]
        x_H2O_t_x0[t_step] = x_H2O_0 + d_d[0] * mbf_red
    return delta_t_x0, x_H2O_t_x0

def reduction_x_t(mu_O_delta_func, mu_O_H2O_func, delta_t_x, x_H2O_t_x,
                  d_delta, d_X, mbf_red, gas_mesh=100, oxide_mesh=100):
    """
    Compute reduction over full space-time reactor grid.

    Args:
        mu_O_delta_func (callable): Function for mu_O in the solid.
        mu_O_H2O_func (callable): Function for mu_O in the gas stream.
        delta_t_x (ndarray): Delta values grid.
        x_H2O_t_x (ndarray): H2O mole fraction grid.
        d_delta (ndarray): Delta differential array.
        d_X (ndarray): Gas composition differential array.
        mbf_red (float): Mass balance factor.
        gas_mesh (int): Number of gas-phase elements (like time steps).
        oxide_mesh (int): Number of oxide elements (spatial x steps).

    Returns:
        tuple: Updated delta and gas composition grids.
    """
    for t_step in range(1, gas_mesh):
        for x_step in range(1, oxide_mesh):
            # Stop if gas is nearly fully oxidized
            if x_H2O_t_x[t_step, x_step - 1] > 1 - x_gas_tolerance:
                x_H2O_t_x[t_step, x_step - 1:] = x_H2O_t_x[t_step, x_step - 1]
                delta_t_x[t_step, x_step - 1:] = delta_t_x[t_step - 1, x_step - 1:]
                break
            d_d = batch_equilibrium_solver.reduction(mu_O_delta_func, mu_O_H2O_func,
                                                     delta_t_x[t_step - 1, x_step],
                                                     x_H2O_t_x[t_step, x_step - 1],
                                                     d_delta, d_X)
            if len(d_d) < 1:
                delta_t_x[t_step, x_step] = delta_t_x[t_step - 1, x_step]
                x_H2O_t_x[t_step, x_step] = x_H2O_t_x[t_step, x_step - 1]
                continue
            delta_t_x[t_step, x_step] = delta_t_x[t_step - 1, x_step] + d_d[0]
            x_H2O_t_x[t_step, x_step] = x_H2O_t_x[t_step, x_step - 1] + d_d[0] * mbf_red
    return delta_t_x, x_H2O_t_x

def compute_reduction_step(mu_O_delta_func, mu_O_H2O_func, x_H2O_0,
                           delta_x_t, x_H2O_x_t, d_delta, d_X, mbf_red,
                           gas_mesh=100, oxide_mesh=100):
    """
    Execute full reduction step over the reactor grid.

    Args:
        mu_O_delta_func (callable): mu_O(delta) for solid.
        mu_O_H2O_func (callable): mu_O(X) for H2O gas.
        x_H2O_0 (float): Initial H2O mole fraction.
        delta_x_t (ndarray): Delta grid.
        x_H2O_x_t (ndarray): H2O mole fraction grid.
        d_delta (ndarray): Delta differential array.
        d_X (ndarray): Gas composition differential array.
        mbf_red (float): Mass balance factor for reduction.
        gas_mesh (int): Number of gas-phase elements (like time steps).
        oxide_mesh (int): Number of oxide elements (spatial x steps).

    Returns:
        tuple: Updated delta and H2O mole fraction grids.
    """
    delta_x_t[0, 0], x_H2O_x_t[0, 0] = reduction_t0_x0_bc(
        mu_O_delta_func, mu_O_H2O_func, x_H2O_0, delta_x_t[0, 0], d_delta, d_X, mbf_red
    )
    delta_x_t[0], x_H2O_x_t[0] = reduction_t0_bc(
        mu_O_delta_func, mu_O_H2O_func, delta_x_t[0], x_H2O_x_t[0],
        d_delta, d_X, mbf_red, oxide_mesh=oxide_mesh
    )
    delta_x_t[:, 0], x_H2O_x_t[:, 0] = reduction_x0_bc(
        mu_O_delta_func, mu_O_H2O_func, x_H2O_0, delta_x_t[:, 0],
        x_H2O_x_t[:, 0], d_delta, d_X, mbf_red, gas_mesh=gas_mesh
    )
    delta_x_t, x_H2O_x_t = reduction_x_t(
        mu_O_delta_func, mu_O_H2O_func, delta_x_t, x_H2O_x_t,
        d_delta, d_X, mbf_red, gas_mesh=gas_mesh, oxide_mesh=oxide_mesh
    )
    return delta_x_t, x_H2O_x_t
