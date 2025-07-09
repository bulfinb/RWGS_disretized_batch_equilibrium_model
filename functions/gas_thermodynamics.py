import numpy as np

# === GAS PHASE THERMODYNAMICS ===
# Polynomial fits derived from Cantera thermodynamic data (573–1273 K), see RWGS_thermo

dG_WS = np.poly1d([2.70345921e-06, -1.13933543e-02, -4.06918554e+01,  2.41896010e+05])  # H2 + 0.5O2 = H2O
dG_CS = np.poly1d([2.32796741e-07,  1.06605461e-03, -8.99506394e+01,  2.84170742e+05])  # CO + 0.5O2 = CO2
dH_WS = np.poly1d([-7.63562376e-07, -1.35614677e-03,  1.13909881e+01,  2.38586737e+05]) # H2 + 0.5O2 = H2O
dH_CS = np.poly1d([1.81197383e-06, -7.59809076e-03,  6.05499810e+00,  2.82355370e+05]) # CO + 0.5O2 = CO2



mu_O2_std = np.poly1d([6.06299815e-06, -3.57632155e-02, -1.90402910e+02, -8.84474990e+02])  # µ_O2 at 1 bar

def pO2_H2O(T, X):
    """
    Calculate effective pO2 from H2/H2O gas mixture.

    Args:
        T (float): Temperature in Kelvin.
        X (float or ndarray): H2 to H2O conversion extent (X = p_H2O / (p_H2 + p_H2O)).

    Returns:
        float or ndarray: Partial pressure of O2 in bar.
    """
    p_O2 = 1 * np.exp(-2 * dG_WS(T) / (8.3145 * T)) * (X / (1 - X)) ** 2
    return p_O2

def pO2_CO2(T, X):
    """
    Calculate effective pO2 from CO/CO2 gas mixture.

    Args:
        T (float): Temperature in Kelvin.
        X (float or ndarray): Mole fraction of CO2 (X = p_CO2 / (p_CO + p_CO2)).

    Returns:
        float or ndarray: Partial pressure of O2 in bar.
    """
    p_O2 = 1 * np.exp(-2 * dG_CS(T) / (8.3145 * T)) * (X / (1 - X)) ** 2
    return p_O2

def mu_O2_H2O(T, X):
    """
    Calculate oxygen chemical potential from H2/H2O system.

    Args:
        T (float): Temperature in Kelvin.
        X (float or ndarray): Mole fraction of H2O.

    Returns:
        float or ndarray: μ_O in kJ/mol.
    """
    mu_O2 = (8.3145 * T * np.log(pO2_H2O(T, X)) + mu_O2_std(T)) / 1000
    return mu_O2

def mu_O2_CO2(T, X):
    """
    Calculate oxygen chemical potential from CO/CO2 system.

    Args:
        T (float): Temperature in Kelvin.
        X (float or ndarray): Mole fraction of CO2.

    Returns:
        float or ndarray: μ_O in kJ/mol.
    """
    mu_O2 = (8.3145 * T * np.log(pO2_CO2(T, X)) + mu_O2_std(T)) / 1000
    return mu_O2


def Cp_H2O(T):
    """Molar specific heat capacity of H₂O (steam) [J/(mol·K)], T in K"""
    return 0.018 * (1745.964 + 0.1851146 * T + 6.194487e-4 * T**2 - 3.026785e-7 * T**3 + 4.19053e-11 * T**4)


def Cp_H2(T):
    """Molar specific heat capacity of H₂ [J/(mol·K)], T in K"""
    return 29.9897 - 3.818e-3 * T + 5.335e-6 * T**2 - 1.1609e-9 * T**3


def Cp_CO2(T):
    """Molar specific heat capacity of CO₂ [J/(mol·K)], T in K"""
    return 20.989 + 6.7679e-2 * T - 4.9598e-5 * T**2 + 1.7795e-8 * T**3 - 2.4949e-12 * T**4


def Cp_CO(T):
    """Molar specific heat capacity of CO [J/(mol·K)], T in K"""
    return 28.983 - 3.4117e-3 * T + 1.4801e-5 * T**2 - 8.9586e-9 * T**3 + 1.6545e-12 * T**4
