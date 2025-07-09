"""
This script compares several oxide materials by computing:
- Redox conversions
- Heat duties (reduction and oxidation)
- Volumetric flow volumes to reactor volume rations ratios
"""
import numpy as np
import math as mt
import time
import matplotlib.pyplot as plt

from functions.mass_balance import calculate_mass_balance
from functions.energy_balance import calculate_energy_balance, oxidation_dT_quasiadiabatic
from functions.main import cycle_until_balanced

# === Material and bed properties ===
materials = ["CeO2", "LSF", "CeZr05", "CeZr15", "CeZr20"]

# Molar densities [mol/m^3] estimated from lattice parameters
molar_density_materials = np.asarray([
    41948,
    27982.0,
    41948 * (5.38 / 5.397) ** 3,
    41948 * (5.376 / 5.397) ** 3,
    41948 * (5.357 / 5.397) ** 3
])

bed_void_fraction = 0.37
granule_void_fraction = 0.5

# === Process conditions ===
T = 1023  # Temperature [K]
p = 15e5  # Pressure [Pa]
n_CO2 = 1.0
n_H2 = 2.0
n_oxide_array = np.asarray([78, 27.3, 16.3, 10.4, 9.9])
x_H2O_0 = 0.005
x_CO2_0 = 0.998

# Reactor and gas volumes
V_reactor = n_oxide_array / (molar_density_materials * (1 - bed_void_fraction) * granule_void_fraction)
V_free = V_reactor * bed_void_fraction
V_CO2 = n_CO2 * 8.314 * T / p
V_ratio = V_CO2 / V_free

# Initialize output arrays
heat_reduction_array = np.zeros(len(materials))
heat_oxidation_array = np.zeros(len(materials))
conversion_CO2_array = np.zeros(len(materials))

# === Simulation loop ===
for i, material in enumerate(materials):
    delta_red, x_H2O_red, delta_ox, x_CO2_ox, cycles = cycle_until_balanced(
        max_cycles=6,
        material=material,
        T=T,
        n_CO2=n_CO2,
        n_H2=n_H2,
        n_oxide=n_oxide_array[i],
        x_CO2_0=x_CO2_0,
        x_H2O_0=x_H2O_0
    )

    X_CO2, X_H2, _, _ = calculate_mass_balance(
        delta_red, x_H2O_red, delta_ox, x_CO2_ox,
        n_CO2=n_CO2,
        n_H2=n_H2,
        n_oxide=n_oxide_array[i],
        x_CO2_0=x_CO2_0,
        x_H2O_0=x_H2O_0
    )
    conversion_CO2_array[i] = X_CO2

    Q_red_abs, Q_ox_rel = calculate_energy_balance(
        T, delta_red, material=material, X_CO2=X_CO2
    )
    heat_reduction_array[i] = Q_red_abs
    heat_oxidation_array[i] = Q_ox_rel

    print(f"T = {T} K | V_ratio = {V_ratio[i]:.3f}  | n_oxide = {n_oxide_array[i]:.1f} | "
          f"X_CO2 = {X_CO2:.4f} | Q_red = {Q_red_abs:.3f} kJ/mol | Q_ox = {Q_ox_rel:.3f} kJ/mol")
    print(n_oxide_array[i], ' & ', Q_red_abs, ' & ', Q_ox_rel, ' & ', V_ratio[i])


