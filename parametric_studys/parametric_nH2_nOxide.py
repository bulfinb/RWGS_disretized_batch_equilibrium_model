import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path

# Get the parent directory
parent_dir = Path(__file__).resolve().parent.parent

# Add parent directory to sys.path
sys.path.append(str(parent_dir))

from functions.mass_balance import calculate_mass_balance
from functions.main import cycle_untill_balanced

"""
This script runs a parametric sweep over H2-to-oxide ratios at constant temperature.
It generates and saves heatmaps for CO2 and H2 conversions and oxygen balance metrics.
"""

# === CONFIGURATION ===
material = "CeZr20"  # Choose from "CeO2", "CeZr05", "CeZr15", "LSF"
T = 973  # Temperature in Kelvin
x_H2O_0 = 0.005  # Initial mole fraction of H2O in H2 stream
x_CO2_0 = 0.998  # Initial mole fraction of CO2 in CO2 stream


# === PARAMETRIC SWEEP RANGES ===
n_oxide_r = np.arange(5, 50, 1)  # Moles of oxide
n_H2_r = np.arange(1.01, 3.01, 0.05)  # Moles of H2

# === INITIALIZE RESULT ARRAYS ===
Map_X_CO2 = np.zeros((len(n_H2_r), len(n_oxide_r)))
Map_X_H2 = np.zeros((len(n_H2_r), len(n_oxide_r)))
O_balance_gas = np.zeros((len(n_H2_r), len(n_oxide_r)))
O_balance_oxide = np.zeros((len(n_H2_r), len(n_oxide_r)))

# === SIMULATION LOOP ===
for i, n_H2 in enumerate(n_H2_r):
    for j, n_oxide in enumerate(n_oxide_r):
        tic = time.perf_counter()

        # Run redox cycles until gas-phase O balance converges
        delta_t_x_red, x_H2O_t_x_red, delta_t_x_ox, x_CO2_t_x_ox, cycles = cycle_untill_balanced(
            material=material,
            T=T,
            n_H2=n_H2,
            n_oxide=n_oxide,
            x_CO2_0=x_CO2_0,
            x_H2O_0=x_H2O_0
        )

        # Evaluate conversion and balance metrics
        X_CO2, X_H2, O_bal_gas, O_bal_oxide = calculate_mass_balance(
            delta_t_x_red, x_H2O_t_x_red,
            delta_t_x_ox, x_CO2_t_x_ox,
            n_H2=n_H2, n_oxide=n_oxide,
            x_CO2_0=x_CO2_0, x_H2O_0=x_H2O_0
        )

        # If oxide is in excess, fill forward and break out
        if cycles >= 20:
            Map_X_H2[i][j - 1:] = X_H2
            Map_X_CO2[i][j - 1:] = X_CO2
            O_balance_oxide[i][j - 1:] = O_bal_oxide
            O_balance_gas[i][j - 1:] = O_bal_gas
            print("Excess oxide detected, skipping remaining n_oxide for this n_H2")
            toc = time.perf_counter()
            print(f"n_ox: {n_oxide:.1f}, n_H2: {n_H2:.3f}, cycles: {cycles}, "
                  f"X_H2: {X_H2:.3f}, X_CO2: {X_CO2:.3f}, MB_gas: {O_bal_gas:.3f}, MB_oxide: {O_bal_oxide:.3f}, t: {toc - tic:.2f}s")
            break

        # Store results
        Map_X_H2[i][j] = X_H2
        Map_X_CO2[i][j] = X_CO2
        O_balance_oxide[i][j] = O_bal_oxide
        O_balance_gas[i][j] = O_bal_gas

        toc = time.perf_counter()
        print(f"n_ox: {n_oxide:.1f}, n_H2: {n_H2:.3f}, cycles: {cycles}, "
              f"X_H2: {X_H2:.3f}, X_CO2: {X_CO2:.3f}, MB_gas: {O_bal_gas:.3f}, MB_oxide: {O_bal_oxide:.3f}, t: {toc - tic:.2f}s")

# === SAVE OUTPUTS ===
os.makedirs("heat_map_data", exist_ok=True)
np.savetxt(Path("heat_map_data") / f"{material}_X_CO2_map_vs_nH2_nOxide_at_T_{T}.csv", Map_X_CO2, delimiter=",")
np.savetxt(Path("heat_map_data") / f"{material}_X_H2_map_vs_nH2_nOxide_at_T_{T}.csv", Map_X_H2, delimiter=",")
np.savetxt(Path("heat_map_data") / f"{material}_O_balance_gas_vs_nH2_nOxide_at_T_{T}.csv", O_balance_gas, delimiter=",")
np.savetxt(Path("heat_map_data") / f"{material}_O_balance_oxide_vs_nH2_nOxide_at_T_{T}.csv", O_balance_oxide, delimiter=",")

# === PLOTTING ===
n_oxide_r, n_H2_r = np.meshgrid(n_oxide_r, n_H2_r)


def plot_heatmap(data, label, filename, levels, fmt):
    fig = plt.figure(figsize=(3.8, 3.5), facecolor='white')
    ax = fig.add_subplot()
    ax.set_xlabel('$n_\\mathrm{oxide}$')
    ax.set_ylabel('$n_\\mathrm{H2}$')
    im = ax.pcolormesh(n_oxide_r, n_H2_r, data, rasterized=True)
    ct = ax.contour(n_oxide_r, n_H2_r, data, levels=levels, colors='black')
    ax.clabel(ct, fmt=fmt)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(label=label)
    filename = Path("heat_map_plots") / filename
    plt.savefig(filename , dpi=400, bbox_inches='tight')
    plt.show()

os.makedirs("heat_map_plots", exist_ok=True)
plot_heatmap(Map_X_CO2, '$X_\\mathrm{CO2}$', f"{material}_X_CO2_map_vs_nH2_nOxide_at_T_{T}.png",
             levels=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.98], fmt='%1.2f')
plot_heatmap(Map_X_H2, '$X_\\mathrm{H_2}$', f"{material}_X_H2_map_vs_nH2_nOxide_at_T_{T}.png",
             levels=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.98], fmt='%1.1f')
plot_heatmap(O_balance_gas, '$O_\\mathrm{balance}$', f"{material}_O_balance_gas_vs_nH2_nOxide_at_T_{T}.png",
             levels=[0.99, 0.999, 1.001, 1.01], fmt='%1.3f')
