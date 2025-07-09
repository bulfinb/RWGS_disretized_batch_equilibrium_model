import numpy as np
import matplotlib.pyplot as plt

from functions.mass_balance import calculate_mass_balance
from functions.main import cycle_until_balanced
from RWGS_thermo.X_CO2_ideal import X_CO2_membrane, X_CO2_catalytic

np.set_printoptions(threshold=10)

"""
This script evaluates CO2 conversion as a function of oxide-to-gas ratio (n_oxide / n_CO2)
for a range of oxide materials in a packed-bed redox reactor at constant temperature.
It also compares results to ideal counterflow and cofeed reactor limits.
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
T = 1023  # Temperature [K]
n_CO2 = 1.0  # Moles of CO2 per cycle
n_H2 = 1.5  # Moles of H2 per cycle
x_CO2_0 = 0.998  # Initial CO2 mole fraction
x_H2O_0 = 0.005  # Initial H2O mole fraction
n_range = np.linspace(2, 40, 30)  # Range of oxide amounts

# === Run simulations
results = []
for material in materials:
    X_CO2_results = []
    X_CO2 = 0
    for n_oxide in n_range:
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
            print(f"{material}, n_oxide = {n_oxide:.2f}, cycles = {cycles}, X_CO2 = {X_CO2:.4f}")

        X_CO2_results.append(X_CO2)
    results.append(X_CO2_results)

# === Ideal Reactor Limits
membrane_max = [X_CO2_membrane(T, n_H2) for _ in n_range]
catalytic_max = [X_CO2_catalytic(T, n_H2) for _ in n_range]

# === Plot Results
fig, ax = plt.subplots(figsize=(4.1, 3.4), facecolor='white')
ax.set_xlabel('$n_{\mathrm{oxide}}$ [-]')
ax.set_ylabel('CO$_2$ Conversion $X_{\mathrm{CO_2}}$ [-]')

for i, material in enumerate(materials):
    ax.plot(n_range, results[i], ls='-', lw=1.2, label=labels[i])

#ax.plot(n_range, membrane_max, ls='-', lw=1.0, color='black', label='Ideal Counterflow')
ax.plot(n_range, catalytic_max, ls='--', lw=1.0, color='black', label='Cofeed')

#ax.set_title(f"$n_{{\mathrm{{H_2}}}}$ = {n_H2}, $T$ = {T} K")
ax.set_xlim(2,40)
ax.set_ylim(0.1,1.01)
ax.legend(bbox_to_anchor=(1.02, 1.0), loc='upper left')
ax.text(0.7, 0.1, '$T =$ '+ str(T) + ' $\, \mathrm{K}$' + '\n' + "$n_\mathrm{CO_2} = $" + str(n_H2) + '\n' + "$n_\mathrm{H_2} = $" + str(n_H2), color ="black",
        transform = ax.transAxes)

# Save and show
plt.savefig("materials_comparison_vs_n_oxide.png", dpi=300, bbox_inches='tight')
plt.show()
