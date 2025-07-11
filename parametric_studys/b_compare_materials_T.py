from pathlib import Path
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add project root to sys.path
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

from functions.mass_balance import calculate_mass_balance
from functions.main import cycle_until_balanced
from RWGS_thermo.X_CO2_ideal import X_CO2_catalytic

"""
This script compares CO2 conversion across several oxide materials as a function of temperature.
It includes comparison against idealized counterflow and cofeed reactors.
"""

# === Materials and Labels
materials = ["CeO2", "LSF", "CeZr05", "CeZr15", "CeZr20"]
labels = [
    "CeO$_2$",
    "LSF",
    "CeZr05",
    "CeZr15",
    "CeZr20"
]

# === Process Conditions
n_CO2 = 1.0  # Moles CO2 per cycle
n_H2 = 1.5  # Moles H2 per cycle
n_oxide = 15  # Moles oxide
x_CO2_0 = 0.998  # Initial CO2 concentration
x_H2O_0 = 0.005  # Initial H2O concentration
T_range = np.linspace(923, 1173, 30)  # Temperature sweep

# === Run simulations for each material
results = []
for material in materials:
    X_CO2_results = []
    X_CO2 = 0
    for T in T_range:
        if X_CO2 < 0.99:
            delta_t_x_red, x_H2O_t_x_red, delta_t_x_ox, x_CO2_t_x_ox, cycles = cycle_until_balanced(
                material=material,
                T=T,
                n_CO2=n_CO2,
                n_H2=n_H2,
                n_oxide=n_oxide,
                x_CO2_0=x_CO2_0,
                x_H2O_0=x_H2O_0
            )

            X_CO2, X_H2, O_bal_gas, O_bal_oxide = calculate_mass_balance(
                delta_t_x_red, x_H2O_t_x_red,
                delta_t_x_ox, x_CO2_t_x_ox,
                n_CO2=n_CO2, n_H2=n_H2, n_oxide=n_oxide,
                x_CO2_0=x_CO2_0, x_H2O_0=x_H2O_0
            )
            print(f"{material}, T = {T:.1f}, cycles = {cycles}, X_CO2 = {X_CO2:.4f}")

        X_CO2_results.append(X_CO2)
    results.append(X_CO2_results)

# === Ideal Coffed Limit
catalytic_max = [X_CO2_catalytic(T, n_H2) for T in T_range]

os.makedirs("publication_plots", exist_ok=True)

# === Plot Results
fig, ax = plt.subplots(figsize=(3.8, 3.4), facecolor='white')
ax.set_xlabel('Temperature [K]')
ax.set_ylabel('CO$_2$ Conversion $X_{\mathrm{CO_2}}$ [-]')

for i, material in enumerate(materials):
    ax.plot(T_range, results[i], ls='-', lw=1.2, label=labels[i])

# Plot Ideal Limits
#ax.plot(T_range, membrane_max, ls='-', lw=1.0, color='black', label='Ideal Counterflow')
ax.plot(T_range, catalytic_max, ls='--', lw=1.0, color='black', label='Ideal Cofeed')
ax.set_xlim(923,1173)
ax.set_ylim(0.1,1.01)
ax.text(0.7, 0.1, "$n_\mathrm{CO_2} = $" + str(n_CO2) + '\n' + "$n_\mathrm{H_2} = $" + str(n_H2) + '\n' + "$n_\mathrm{oxide} = $" + str(n_oxide), color ="black",
        transform = ax.transAxes)

# Save and show plot
filename = Path("publication_plots") / "materials_comparison_vs_T.png"
plt.savefig(filename, dpi=300, bbox_inches='tight')
plt.show()
