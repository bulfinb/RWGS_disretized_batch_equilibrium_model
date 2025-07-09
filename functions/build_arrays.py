import numpy as np

def X_array(n_points=2000):
    """
    Generate a non-linear gas composition array spanning 0 to 1.

    The array is dense near 0 and 1 to better capture edge behavior
    during gas phase equilibrium.

    Args:
        n_points (int): Total number of array points.

    Returns:
        ndarray: Gas composition extent array.
    """
    X_half = np.logspace(-6, np.log10(0.5), int(n_points / 2))
    X_half_2 = 1.0 - X_half
    X = np.concatenate((X_half, np.flip(X_half_2)))
    return X


def delta_array(delta_max=0.5, n_points=2000):
    """
    Generate a delta array with finer resolution near zero.

    Combines log spacing and linear spacing to efficiently sample
    the delta space for materials with sharp equilibrium behavior.

    Args:
        delta_max (float): Maximum delta value.
        n_points (int): Total number of points in the array.

    Returns:
        ndarray: Discretized delta values.
    """
    delta_log = np.logspace(-6.0, np.log10(delta_max / 1.5), int(n_points / 2))
    delta_lin = np.linspace(delta_log[-1], delta_max, int(n_points / 2))
    delta_range = np.concatenate((delta_log, delta_lin))
    return delta_range


def set_delta_max_and_delta_range(material, n_points=2000):
    """
    Assign delta_max and build corresponding delta array based on material.

    Args:
        material (str): Material name (e.g., 'LSF', 'CeO2').
        n_points (int): Number of points for the delta array.

    Returns:
        tuple:
            delta_max (float): Maximum delta value for material.
            delta_range (ndarray): Delta values from 0 to delta_max.

    Raises:
        ValueError: If material name is not recognized.
    """
    if material == "LSF":
        delta_max = 0.5
    elif material in ["CeO2", "CeO2_simple"]:
        delta_max = 0.335
    elif material == "CeZr05":
        delta_max = 0.405
    elif material == "CeZr15":
        delta_max = 0.424
    elif material == "CeZr20":
        delta_max = 0.399
    else:
        raise ValueError("Material not recognized. Must be 'LSF', 'CeO2', 'CeO2_simple', 'CeZr05','CeZr15' or 'CeZr20'.")

    delta_range = delta_array(delta_max=delta_max, n_points=n_points)
    return delta_max, delta_range


def mass_balance_arrays(
    delta_max=0.5, n_gas=1.0, n_oxide=20.0,
    gas_mesh=100, oxide_mesh=100, n_points=400
):
    """
    Construct d_delta and d_X arrays for probing equilibrium conditions.

    The arrays represent potential shifts in delta and gas composition
    constrained by stoichiometry and mesh ratios.

    Args:
        delta_max (float): Maximum possible delta shift.
        n_gas (float): Moles of gas phase.
        n_oxide (float): Moles of oxide phase.
        gas_mesh (int): Number of gas discretization steps.
        oxide_mesh (int): Number of solid discretization steps.
        n_points (int): Number of points in d_delta and d_X arrays.

    Returns:
        tuple:
            d_delta_array (ndarray): Array of delta shift values.
            d_X_array (ndarray): Corresponding gas composition shifts.
    """
    d_delta_max_oxide = delta_max
    d_delta_max_gas = 1 / ((n_oxide / n_gas) * (gas_mesh / oxide_mesh))
    d_delta_max = min(d_delta_max_oxide, d_delta_max_gas)

    a = np.logspace(-7, np.log10(d_delta_max / 1.5), int(n_points / 2))
    b = np.linspace(a[-1], d_delta_max, int(n_points / 2))
    d_delta_array = np.concatenate(([0.0], a, b))

    d_X_array = d_delta_array * (n_oxide / n_gas) * (gas_mesh / oxide_mesh)

    return d_delta_array, d_X_array
