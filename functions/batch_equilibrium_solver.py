import numpy as np

def reduction(mu_O_delta_func, mu_O_H2O_func, delta_i, x_H2O_i, d_delta, d_X):
    """
    Solve equilibrium condition for batch reactor reduction.

    Finds the intersection of solid and gas phase chemical potential
    curves under mass balance constraints.

    Args:
        mu_O_delta_func (callable): Solid-phase oxygen chemical potential as a function of delta.
        mu_O_H2O_func (callable): Gas-phase oxygen chemical potential as a function of H2O mole fraction.
        delta_i (float): Initial delta value.
        x_H2O_i (float): Initial H2O mole fraction.
        d_delta (ndarray): Array of delta value shifts to probe.
        d_X (ndarray): Corresponding array of gas composition shifts.

    Returns:
        ndarray: Array of delta shifts that satisfy the equilibrium condition.

    Raises:
        ValueError: If no intersection of chemical potential curves is found.
    """
    diff_mu_O = mu_O_delta_func(delta_i + d_delta) - mu_O_H2O_func(x_H2O_i + d_X)
    idx = np.argwhere(np.diff(np.sign(diff_mu_O))).flatten()
    d_d = d_delta[idx] - diff_mu_O[idx] / (
        (diff_mu_O[idx + 1] - diff_mu_O[idx]) / (d_delta[idx + 1] - d_delta[idx])
    )
    return d_d

def oxidation(mu_O_delta_func, mu_O_CO2_func, delta_i, x_CO2_i, d_delta, d_X):
    """
    Solve equilibrium condition for batch reactor oxidation.

    Finds the intersection of chemical potential curves between
    the reduced solid and oxidizing gas stream.

    Args:
        mu_O_delta_func (callable): Solid-phase oxygen chemical potential as a function of delta.
        mu_O_CO2_func (callable): Gas-phase oxygen chemical potential as a function of CO2 mole fraction.
        delta_i (float): Initial delta value.
        x_CO2_i (float): Initial CO2 mole fraction.
        d_delta (ndarray): Array of delta value shifts to probe.
        d_X (ndarray): Corresponding array of gas composition shifts.

    Returns:
        ndarray: Array of delta shifts that satisfy the equilibrium condition.

    Raises:
        ValueError: If no equilibrium point (intersection of chemical potential curves) is found.
    """
    diff_mu_O = -mu_O_delta_func(delta_i - d_delta) + mu_O_CO2_func(x_CO2_i - d_X)
    idx = np.argwhere(np.diff(np.sign(diff_mu_O))).flatten()

    if len(idx) == 0:
        raise ValueError("No equilibrium point found: no intersection in chemical potential curves.")

    d_d = d_delta[idx] - diff_mu_O[idx] / (
        (diff_mu_O[idx + 1] - diff_mu_O[idx]) / (d_delta[idx + 1] - d_delta[idx])
    )
    return d_d

