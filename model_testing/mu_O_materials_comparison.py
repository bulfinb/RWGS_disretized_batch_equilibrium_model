import os
import sys
import matplotlib.pyplot as plt
from pathlib import Path

# Get the parent directory for importing model functions
parent_dir = Path(__file__).resolve().parent.parent

# Add parent directory to sys.path
sys.path.append(str(parent_dir))

from functions import build_arrays
from functions import build_mu_O_functions



# === SETUP
T = 1073
materials = ["CeO2", "LSF", "CeZr05", "CeZr15", "CeZr20"]

labels = [
    "CeO$_2$",
    "LSF",
    "CeZr05",
    "CeZr15",
    "CeZr20"
]
mu_O_delta_material = []
delta_ranges = []
for material in materials:
    delta_max, delta_range = build_arrays.set_delta_max_and_delta_range(material)
    delta_ranges.append(delta_range)
    # Generate mu_O(delta) function for the oxide and its inverse
    mu_O_delta_func, mu_O_delta_func_inv = build_mu_O_functions.set_mu_O_delta_functions(material, T, delta_range)
    mu_O_delta_material.append(mu_O_delta_func)
    # Generate X (gas composition) range
    X_range = build_arrays.X_array()

    # Generate mu_O(X) functions for gas streams
    mu_O_CO2_func = build_mu_O_functions.set_mu_O_CO2_function(T, X_range)
    mu_O_H2O_func = build_mu_O_functions.set_mu_O_H2O_function(T, X_range)


# === Plot Results
fig, ax = plt.subplots(figsize=(4.5, 4.0), facecolor='white')
ax.set_xlabel('$x_{\mathrm{CO_2}}$, $2\delta$   [-]')
ax.set_ylabel('Chemical potential $\mu_\\mathrm{O}$ [$\\mathrm{kJ \, mol^{-1} }$]')

for i, mu_O_delta_m in enumerate(mu_O_delta_material):
    if i == 1:
        ax.plot(2*(delta_ranges[i]-0.2), mu_O_delta_m(delta_ranges[i]), ls='--', lw=1.0, label=labels[i])
    else:
        ax.plot(2 * (delta_ranges[i]), mu_O_delta_m(delta_ranges[i]), ls='--', lw=1.0, label=labels[i])
# Plot the chemical potential profile in the CO2 flow
ax.plot(1-X_range, mu_O_CO2_func(X_range), ls='-', lw=1.2, color='black', label='CO$_2$')

# Plot the chemical potential profile in the H2 flow
ax.plot(1-X_range, mu_O_H2O_func(0.5*X_range), ls='-', lw=1.2, color='grey', label='H$_2$')


# Title and Legend
ax.legend( loc='upper right')

ax.set_xlim([0, 1])
ax.set_ylim([-380, -225])


os.makedirs("plots", exist_ok=True)
filename = Path("plots") / "mu_O_materials.png"
# Save and show plot
plt.savefig(filename, dpi=150, bbox_inches='tight')
plt.show()