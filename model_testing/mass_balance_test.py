import time
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import matplotlib as mpl
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color= plt.cm.viridis(np.linspace(0, 1, 21)))

from pathlib import Path

# Add project root to sys.path to import project modules
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

from functions.mass_balance import calculate_mass_balance
from functions.main import simulate_cycle

"""
This script runs a countercurrent chemical looping model to estimate thermodynamic limits.
It simulates multiple redox cycles and visualizes the evolution of non-stoichiometry (delta),
reactor conversions, and oxygen mass balances.
"""

# === Set material from options - "CeO2", "CeZr05", "CeZr15", "CeZr20", "LSF"
material = "CeZr15"

# === Define process conditions
T = 1023  # Temperature in Kelvin
n_CO2 = 1.0  # Moles of CO2 per oxidation cycle
n_H2 = 1.5  # Moles of H2 per reduction cycle
n_oxide = 20  # Moles of oxide material
cycles = 20  # Number of redox cycles
x_CO2_0 = 0.998  # Initial CO2 mole fraction
x_H2O_0 = 0.005  # Initial H2O mole fraction

# === Initialize storage for plotting variables
delta_red_end = []
delta_ox_end = []
X_CO2_values = []
X_H2_values = []
O_bal_gas_values = []
O_bal_oxide_values = []
O_bal_oxide_values_red = []

# === Initial state of the oxide bed (uniform delta profile)
delta_x_0 = np.zeros(100)
first_cycle = True

# === Run simulation for defined number of cycles
for cycle in range(cycles):
    tic = time.perf_counter()

    # Get previous oxidation profile or initial state
    if first_cycle:
        delta_ox_prev = delta_x_0
    else:
        delta_ox_prev = delta_t_x_ox[-1]

    # Run one redox cycle (reduction followed by oxidation)
    delta_t_x_red, x_H2O_t_x_red, delta_t_x_ox, x_CO2_t_x_ox = simulate_cycle(
        material=material,
        first_cycle=first_cycle,
        delta_x_0=delta_x_0,
        x_CO2_0=x_CO2_0,
        x_H2O_0=x_H2O_0,
        T=T,
        n_H2=n_H2,
        n_oxide=n_oxide
    )

    # Update flag and initial condition for next cycle
    first_cycle = False
    delta_x_0 = delta_t_x_ox[-1]

    # Compute conversions and oxygen balances
    X_CO2, X_H2, O_bal_gas, O_bal_oxide = calculate_mass_balance(
        delta_t_x_red, x_H2O_t_x_red,
        delta_t_x_ox, x_CO2_t_x_ox,
        x_CO2_0=x_CO2_0,
        x_H2O_0=x_H2O_0,
        n_H2=n_H2,
        n_oxide=n_oxide
    )

    # Change in delta during reduction step
    d_delta_CeO2 = delta_t_x_red[-1] - delta_ox_prev
    nO_CeO2_red = d_delta_CeO2.sum() * n_oxide / len(delta_t_x_red[-1])

    # H2O produced and compute oxygen balance during reduction
    n_H2O = (x_H2O_t_x_red.T[-1].mean() - x_H2O_0) * n_H2
    O_Bal_oxide_red = n_H2O / nO_CeO2_red

    # Store simulation results
    delta_red_end.append(delta_t_x_red[-1])
    delta_ox_end.append(delta_t_x_ox[-1])
    X_CO2_values.append(X_CO2)
    X_H2_values.append(X_H2)
    O_bal_gas_values.append(O_bal_gas)
    O_bal_oxide_values.append(O_bal_oxide)
    O_bal_oxide_values_red.append(O_Bal_oxide_red)

    toc = time.perf_counter()
    print(f"Cycle {cycle}: solve time {toc - tic:.2f} s, ", O_bal_oxide, O_bal_gas)

# the first reduction mass balance is not valid.
O_bal_oxide_values_red.pop(0)

print("Final CO2 conversion:", X_CO2_values[-1])


# === PLOTTING ===
os.makedirs("plots", exist_ok=True)

# === Plot 1: Delta profiles along the bed for all cycles
x_space = np.arange(len(delta_t_x_red[-1]))+1
fig, ax = plt.subplots(figsize=(4.0, 3.7))
ax.set_xlabel('Oxide element # [-]')
ax.set_ylabel('Non-stoichiometry $\\delta$ [-]')

ax.plot(x_space, delta_red_end[0], ls='-', lw=1.2,  label='End of reduction #1')


# Plot reduction delta profiles
for i in range(1, cycles-1):
    ax.plot(x_space, delta_red_end[i+1], ls='-', lw=1.0, label='__nolegend__')

# Reset color cycle for clarity
plt.gca().set_prop_cycle(None)

ax.plot(x_space, delta_ox_end[0], ls='--', lw=1.2,  label='End of oxidation #1')
# Plot oxidation delta profiles
for i in range(1, cycles-1):
    ax.plot(x_space, delta_ox_end[i+1], ls='--', lw=1.0, label='__nolegend__')

# Highlight final cycle
ax.plot(x_space, delta_red_end[-1], ls='-', lw=1.2,  label='End of reduction #20')
ax.plot(x_space, delta_ox_end[-1], ls='--', lw=1.2,  label='End of oxidation #20')
ax.set_xlim(0,100)
ax.set_ylim(0,0.2)
plt.legend()
filename = Path("plots") / "mass_balance_20_cycle_delta_profiles.png"
plt.savefig(filename, dpi=300, bbox_inches='tight')

plt.show()

cycles = np.linspace(1,20,20)
# === Plot 2: System mass balance
fig, ax = plt.subplots(figsize=(4.0, 3.7))
ax.set_xlabel('Cycle # [-]')
ax.set_ylabel('Mass Balance [-]')
ax.plot(np.linspace(2,20,19), O_bal_oxide_values_red, ls='--', lw=1.3, color='tab:orange', label='Reduction: H$_2$O/$\\Delta \delta$')
ax.plot(cycles, O_bal_oxide_values, ls=':', lw=2, color='tab:blue', label='Oxidation: CO/$\\Delta \delta$')
ax.plot(cycles, O_bal_gas_values, ls='--', lw=1.3, color = 'tab:green', label='Cycle: H$_2$O/CO')
ax.plot(cycles, X_H2_values, ls='-', color = 'grey', lw=1.3, label='$X_{\mathrm{H_2}}$')
ax.plot(cycles, X_CO2_values, ls='-', lw=1.3, color = 'black', label='$X_{\mathrm{CO_2}}$')
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_xlim(1,20)
ax.legend()
filename = Path("plots") / "mass_balance_20_cycles_O_balance.png"
plt.savefig(filename, dpi=300, bbox_inches='tight')
plt.show()
