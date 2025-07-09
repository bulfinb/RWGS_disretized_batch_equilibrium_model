import numpy as np
from functions.oxide_thermodynamics import (
    Ce_Dh, CeZr05_Dh, CeZr15_Dh, CeZr20_Dh, LSF_Dh,
    Cm_CeO2, Cm_LSF
)
from functions.gas_thermodynamics import (
    dH_WS, dH_CS, Cp_H2, Cp_H2O, Cp_CO2, Cp_CO
)


def calculate_energy_balance(T, delta_t_x_red, material="CeO2",
                              F_CO2_in=1, X_CO2=0.95):
    """
    Calculate the heat absorbed during reduction and released during oxidation.

    Args:
        T (float): Temperature in Kelvin.
        delta_t_x_red (ndarray): Oxygen non-stoichiometry grid after reduction.
        material (str): Oxide material identifier.
        F_CO2_in (float): Inlet molar flow rate of CO2 [mol/s]. Default is 1.
        X_CO2 (float): CO2 conversion. Default is 0.95.

    Returns:
        tuple:
            Q_red_absorbed (float): Heat absorbed during reduction [kW].
            Q_ox_released (float): Heat released during oxidation [kW].
    """
    if material == "LSF":
        dH_red = LSF_Dh(delta_t_x_red, T)
    elif material in ["CeO2", "CeO2_simple"]:
        dH_red = Ce_Dh(delta_t_x_red)
    elif material == "CeZr05":
        dH_red = CeZr05_Dh(delta_t_x_red)
    elif material == "CeZr15":
        dH_red = CeZr15_Dh(delta_t_x_red)
    elif material == "CeZr20":
        dH_red = CeZr20_Dh(delta_t_x_red)
    else:
        raise ValueError(f"Unknown material: {material}")

    dH_red_mean = np.mean(dH_red)

    # Endothermic heat flow during reduction (heat input)
    Q_red_absorbed = F_CO2_in * (dH_red_mean - dH_WS(T)) * X_CO2 / 1000.0  # kW

    # Exothermic heat flow during oxidation (heat release)
    Q_ox_released = F_CO2_in * (dH_CS(T) - dH_red_mean) * X_CO2 / 1000.0  # kW

    return Q_red_absorbed, Q_ox_released


def oxidation_dT_quasiadiabatic(T, Q_ox_released, t_ox, delta_t_x_red,
                                 material="CeO2", F_CO2_in=1,
                                 X_CO2=0.95, n_oxide=10.0):
    """
    Estimate the temperature rise during quasi-adiabatic oxidation.

    Args:
        T (float): Temperature in Kelvin.
        Q_ox_released (float): Heat released during oxidation [kW].
        t_ox (float): Duration of oxidation step [s].
        delta_t_x_red (ndarray): Oxygen non-stoichiometry after reduction.
        material (str): Oxide material.
        F_CO2_in (float): CO2 molar flow rate [mol/s].
        X_CO2 (float): CO2 conversion.
        n_oxide (float): Moles of oxide in the reactor.

    Returns:
        dT (float): Estimated temperature change [K].
    """
    # Select appropriate molar heat capacity
    if material == "LSF":
        Cm_oxide = Cm_LSF(T, np.mean(delta_t_x_red))
    elif material in ["CeO2", "CeO2_simple", "CeZr05", "CeZr15", "CeZr20"]:
        Cm_oxide = Cm_CeO2(T) * (3 - np.mean(delta_t_x_red))
    else:
        raise ValueError(f"Unknown material: {material}")

    # Heat capacities of gas and oxide
    C_gas_ox = t_ox * F_CO2_in * (X_CO2 * Cp_CO(T) + (1 - X_CO2) * Cp_CO2(T))
    C_oxide_red = n_oxide * Cm_oxide

    # Total energy released during oxidation (J)
    Q_ox_total = Q_ox_released * t_ox * 1000  # convert kW·s to J

    # Adiabatic temperature change
    dT = -Q_ox_total / (C_gas_ox + C_oxide_red)

    return dT
