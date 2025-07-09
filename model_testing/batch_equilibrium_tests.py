from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add project root to sys.path
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

from functions import build_arrays
from functions import build_mu_O_functions
from functions import batch_equilibrium_solver

material = "CeO2"
T= 1073
gas_mesh = 100
oxide_mesh = 100
n_oxide = 10
n_CO2 = 1.0
n_H2 = 1.0
x_H2O_0 = 0.002
x_CO2_0 = 0.998
delta_i_red = 0.001
delta_i_ox = 0.1

mbf_red = (n_oxide / n_H2) * (gas_mesh / oxide_mesh)
# Set delta_max and delta range
delta_max, delta_range = build_arrays.set_delta_max_and_delta_range(material)
mu_O_delta_func, mu_O_delta_func_inv = build_mu_O_functions.set_mu_O_delta_functions(material, T, delta_range)

# Generate mu_O(delta) function for the oxide and its inverse
mu_O_delta_func, mu_O_delta_func_inv = build_mu_O_functions.set_mu_O_delta_functions(material, T, delta_range)

# Generate X (gas composition) range
X_range = build_arrays.X_array()

# Generate mu_O(X) functions for gas streams
mu_O_CO2_func = build_mu_O_functions.set_mu_O_CO2_function(T, X_range)
mu_O_H2O_func = build_mu_O_functions.set_mu_O_H2O_function(T, X_range)


# Create differential change in oxygen mass balance arrays for batch equilibrium solver
d_delta, d_X = build_arrays.mass_balance_arrays(
    delta_max=delta_max, n_gas=n_H2,
    n_oxide=n_oxide, gas_mesh=gas_mesh, oxide_mesh=oxide_mesh
)

Delta_delta = batch_equilibrium_solver.reduction(mu_O_delta_func, mu_O_H2O_func, delta_i_red, x_H2O_0, d_delta, d_X)
Delta_delta_2 = batch_equilibrium_solver.oxidation(mu_O_delta_func, mu_O_CO2_func, delta_i_ox, x_CO2_0, d_delta, d_X)

# VISUALIZE REDUCTION EQUILIBRIUM - This plot is figure 3 in the paper
fig, ax = plt.subplots(figsize=(4.5, 3.5), facecolor='white')
ax.set_xlabel('$\Delta \delta$ [-]')
ax.set_ylabel('$\\mu_O$ [$\\mathrm{kJ \, mol^{-1}}$]')
ax.plot(d_delta, mu_O_H2O_func(x_H2O_0 + d_X), label = "$\\mu_O(x_\\mathrm{H2O})$ gas")
ax.plot(d_delta, mu_O_delta_func(delta_i_red + d_delta), label="$\\mu_O(\\delta)$ oxide")
ax.plot(Delta_delta[0], mu_O_delta_func(delta_i_red + Delta_delta[0]), lw= 0.0,
        color = 'black', ms = 9, marker ='x', label = '$\\Delta \\delta $ equilibrium')
ax.set_ylim((-350,-260))
ax.set_xlim(0, 0.105)
plt.legend()
filename =  Path("plots") / "mu_O_equilibrium_red.png"
plt.savefig(filename, dpi=300, bbox_inches='tight')
plt.show()



# VISUALIZE OXIDATION EQUILIBRIUM
fig, ax = plt.subplots(figsize=(4.5, 3.5), facecolor='white')
ax.set_xlabel('$\Delta \delta$ [-]')
ax.set_ylabel('$\\mu_O$ [$\\mathrm{kJ \, mol^{-1}}$]')
ax.plot(d_delta, mu_O_CO2_func(x_CO2_0 -  d_X), label = "$\\mu_O(x_\\mathrm{CO2})$ gas")
ax.plot(d_delta, mu_O_delta_func( delta_i_ox- d_delta), label="$\\mu_O(\\delta)$ oxide")
ax.plot(Delta_delta_2[0], mu_O_delta_func(delta_i_ox -  Delta_delta_2[0]), lw= 0.0,
        color = 'black', ms = 9, marker ='x', label = '$\\Delta \\delta $ equilibrium')
ax.set_ylim((-350,-260))
ax.set_xlim(0, 0.105)
plt.legend()
filename =  Path("plots") / "mu_O_equilibrium_ox.png"
plt.savefig(filename, dpi=300, bbox_inches='tight')
plt.show()




