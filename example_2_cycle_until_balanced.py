"""
This script runs cycles until the bed reaches steady state from cycle to cycle,
i.e. 0.999 < n_CO/n_H2O < 1.001
"""

from functions.mass_balance import calculate_mass_balance
from functions.energy_balance import calculate_energy_balance, oxidation_dT_quasiadiabatic
from functions.main import cycle_until_balanced

# === Set material from options - "CeO2", "CeZr05", "CeZr15", "CeZr20", "LSF"
material = "CeZr15"

# === Define process conditions
T = 1023  # Temperature in Kelvin
n_CO2 = 1.0  # Moles of CO2 per oxidation cycle
n_H2 = 1.5  # Moles of H2 per reduction cycle
n_oxide = 15  # Moles of oxide material
x_CO2_0 = 0.998  # Initial CO2 mole fraction
x_H2O_0 = 0.005  # Initial H2O mole fraction


# === Simulate cycles===
delta_red, x_H2O_red, delta_ox, x_CO2_ox, cycles = cycle_until_balanced(
    max_cycles=50,
    material=material,
    T=T,
    n_CO2=n_CO2,
    n_H2=n_H2,
    n_oxide=n_oxide,
    x_CO2_0=x_CO2_0,
    x_H2O_0=x_H2O_0
)

# === Post processing ===
X_CO2, X_H2, O_bal_gas, O_bal_oxide= calculate_mass_balance(
    delta_red, x_H2O_red, delta_ox, x_CO2_ox,
    n_CO2=n_CO2,
    n_H2=n_H2,
    n_oxide=n_oxide,
    x_CO2_0=x_CO2_0,
    x_H2O_0=x_H2O_0
)


Q_red_abs, Q_ox_rel = calculate_energy_balance(
    T, delta_red, material=material, X_CO2=X_CO2
)

print(f" # of cycles till balanced= {cycles} | O balance oxide = {O_bal_oxide:.6f} | CO/H2O = {O_bal_gas:.5f} ")
print(f"X_CO2 = {X_CO2:.4f} | X_H2 = {X_H2:.4f} | Q_red = {Q_red_abs:.3f} kJ/mol | Q_ox = {Q_ox_rel:.3f} kJ/mol")