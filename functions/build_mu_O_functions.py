import numpy as np
from scipy import interpolate
from functions.gas_thermodynamics import mu_O2_H2O, mu_O2_CO2
from functions.oxide_thermodynamics import (
    Ce_mu_O, Ce_simple_mu_O, CeZr05_mu_O, CeZr15_mu_O, CeZr20_mu_O, LSF_mu_O
)

def set_mu_O_delta_functions(material, T, delta_range):
    """
    Generate interpolation functions for μ_O(delta) and delta(μ_O).

    Args:
        material (str): Material name ("LSF", "CeO2", etc.).
        T (float): Temperature in Kelvin.
        delta_range (ndarray): Delta values for interpolation.

    Returns:
        tuple:
            mu_O_delta_func (interp1d): Interpolated μ_O(delta).
            mu_O_delta_func_inv (interp1d): Interpolated inverse delta(μ_O).

    Raises:
        ValueError: If the material is not supported.
    """
    if material == "LSF":
        mu_O_delta = LSF_mu_O(delta_range, T)
    elif material == "CeO2":
        mu_O_delta = Ce_mu_O(delta_range, T)
    elif material == "CeO2_simple":
        mu_O_delta = Ce_simple_mu_O(delta_range, T)
    elif material == "CeZr05":
        mu_O_delta = CeZr05_mu_O(delta_range, T)
    elif material == "CeZr15":
        mu_O_delta = CeZr15_mu_O(delta_range, T)
    elif material == "CeZr20":
        delta_range = np.clip(delta_range, 0.001, 0.399)
        mu_O_delta = CeZr20_mu_O(delta_range, T)
    else:
        raise ValueError("Material not recognized. Must be 'LSF', 'CeO2', 'CeO2_simple', 'CeZr05', 'CeZr15', 'CeZr20'.")

    mu_O_delta_func = interpolate.interp1d(
        delta_range, mu_O_delta,
        fill_value=(mu_O_delta.max(), mu_O_delta.min()),
        bounds_error=False
    )
    mu_O_delta_func_inv = interpolate.interp1d(mu_O_delta, delta_range)
    return mu_O_delta_func, mu_O_delta_func_inv

def set_mu_O_H2O_function(T, X_range):
    """
    Generate interpolation function for μ_O(X) in the H2/H2O gas stream.

    Args:
        T (float): Temperature in Kelvin.
        X_range (ndarray): Array of gas composition extents.

    Returns:
        interp1d: Interpolated μ_O(X) for the H2O gas phase.
    """
    mu_O_H2O = 0.5 * mu_O2_H2O(T, X_range)
    mu_O_H2O_func = interpolate.interp1d(
        X_range, mu_O_H2O,
        fill_value=(mu_O_H2O.min(), mu_O_H2O.max()),
        bounds_error=False
    )
    return mu_O_H2O_func

def set_mu_O_CO2_function(T, X_range):
    """
    Generate interpolation function for μ_O(X) in the CO/CO2 gas stream.

    Args:
        T (float): Temperature in Kelvin.
        X_range (ndarray): Array of gas composition extents.

    Returns:
        interp1d: Interpolated μ_O(X) for the CO2 gas phase.
    """
    mu_O_CO2 = 0.5 * mu_O2_CO2(T, X_range)
    mu_O_CO2_func = interpolate.interp1d(
        X_range, mu_O_CO2,
        fill_value=(mu_O_CO2.min(), mu_O_CO2.max()),
        bounds_error=False
    )
    return mu_O_CO2_func


