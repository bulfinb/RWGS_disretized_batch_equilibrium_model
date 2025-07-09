from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sys
import time

# Add project root to sys.path
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

from functions.mass_balance import calculate_mass_balance
from functions.main import simulate_cycle


# === Set material from options - "CeO2" "CeZr05" "CeZr15" "CeZr20" "LSF"
material = "CeO2"

# === Define process conditions
T = 1063  # Temperature in Kelvin
n_CO2 = 0.00825  # Moles of CO2 per oxidation cycle
n_H2 = 0.0121  # Moles of H2 per reduction cycle
m_CeO2 = 83.2 # grams of oxide
n_oxide =  m_CeO2/172.11 # Moles of oxide material
cycles = 10  # Number of redox cycles
x_CO2_0 = 0.998  # Initial CO2 mole fraction
x_H2O_0 = 0.005  # Initial H2O mole fraction

# === Initialize storage for plotting variables
delta_red_end = []
delta_ox_end = []
X_CO2_values = []
X_H2_values = []
O_bal_gas_values = []
O_bal_oxide_values = []

# === Initial state of the oxide bed
delta_x_0 = np.zeros(100)
first_cycle = True

# === Run simulation for defined number of cycles
for cycle in range(cycles):
    tic = time.perf_counter()

    # Run redox simulation cycle
    delta_t_x_red, x_H2O_t_x_red, delta_t_x_ox, x_CO2_t_x_ox = simulate_cycle(
        material=material,
        first_cycle=first_cycle,
        delta_x_0=delta_x_0,
        x_CO2_0=x_CO2_0,
        x_H2O_0=x_H2O_0,
        T=T,
        n_CO2=n_CO2,
        n_H2=n_H2,
        n_oxide=n_oxide
    )

    # Update initial condition for next cycle
    first_cycle = False
    delta_x_0 = delta_t_x_ox[-1]

    # Compute mass balances and conversions
    X_CO2, X_H2, O_bal_gas, O_bal_oxide = calculate_mass_balance(
        delta_t_x_red, x_H2O_t_x_red,
        delta_t_x_ox, x_CO2_t_x_ox,
        x_CO2_0=x_CO2_0,
        x_H2O_0=x_H2O_0,
        n_H2=n_H2,
        n_CO2 = n_CO2,
        n_oxide=n_oxide
    )

    # Store results for plotting
    delta_red_end.append(delta_t_x_red[-1])
    delta_ox_end.append(delta_t_x_ox[-1])
    X_CO2_values.append(X_CO2)
    X_H2_values.append(X_H2)
    O_bal_gas_values.append(O_bal_gas)
    O_bal_oxide_values.append(O_bal_oxide)

    toc = time.perf_counter()
    print(f"Cycle {cycle}: solve time {toc - tic:.2f} s, ", O_bal_oxide, O_bal_gas)

print("Final CO2 conversion:", X_CO2)

# === Plot 1: Delta profiles over reactor bed
x_space = np.arange(len(delta_t_x_red[-1]))
fig, ax = plt.subplots(figsize=(4.5, 4))
ax.set_xlabel('Oxide element # [-]')
ax.set_ylabel('Non-stoichiometry $\\delta$ [-]')

# Plot each cycle
for i in range(cycles):
    ax.plot(x_space, delta_red_end[i], ls='-', lw=1.0, label=f'Reduction Cycle {i+1}')
plt.gca().set_prop_cycle(None)
for i in range(cycles):
    ax.plot(x_space, delta_ox_end[i], ls='--', lw=1.0, label=f'Oxidation Cycle {i+1}')

# Highlight final cycle in black
ax.plot(x_space, delta_red_end[-1], ls='-', lw=1.0, color='black', label='Final Reduction')
ax.plot(x_space, delta_ox_end[-1], ls='--', lw=1.0, color='black', label='Final Oxidation')

plt.tight_layout()
plt.show()

# === Plot 2: CO2 and H2 conversion vs cycle number
fig, ax = plt.subplots(figsize=(4.5, 4))
ax.set_xlabel('Cycle Number [-]')
ax.set_ylabel('Conversion $X_i$ [-]')
ax.plot(range(cycles), X_H2_values, ls='-', lw=1.0, label='$X_{\mathrm{H_2}}$')
ax.plot(range(cycles), X_CO2_values, ls='-', lw=1.0, label='$X_{\mathrm{CO_2}}$')
ax.legend()
plt.tight_layout()
plt.show()

# === Plot 3: Oxygen balance in gas and solid
fig, ax = plt.subplots(figsize=(4.5, 4))
ax.set_xlabel('Cycle Number [-]')
ax.set_ylabel('Oxygen Mass Balance [-]')
ax.plot(range(cycles), O_bal_gas_values, ls='-', lw=1.0, label='Gas Phase: H$_2$O/CO')
ax.plot(range(cycles), O_bal_oxide_values, ls='-', lw=1.0, label='Solid Phase: CO/$\\Delta \delta$')
ax.legend()
plt.tight_layout()
plt.show()
