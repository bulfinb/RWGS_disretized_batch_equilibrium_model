import sys
import os
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path

"""
This script plots 3 parametric sweeps used in the publication
"""



# === LOAD THE DATA ===
CeO2 = np.loadtxt(Path("heat_map_data") / f"CeO2_X_CO2_map_vs_nH2_nOxide_at_T_973.csv", delimiter=",")
LSF = np.loadtxt(Path("heat_map_data") / f"LSF_X_CO2_map_vs_nH2_nOxide_at_T_973.csv",  delimiter=",")
CeZr20 = np.loadtxt(Path("heat_map_data") / f"CeZr20_X_CO2_map_vs_nH2_nOxide_at_T_973.csv", delimiter=",")


# === PLOTTING ===
n_oxide_r = np.arange(5, 50, 1)  # Moles of oxide
n_H2_r = np.arange(1.01, 3.01, 0.05)  # Moles of H2
n_oxide_r, n_H2_r = np.meshgrid(n_oxide_r, n_H2_r)

levels=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.98]

# PLOT 1
fig = plt.figure(figsize=(3.4, 3.3), facecolor='white')
ax = fig.add_subplot()
ax.set_xlabel('$n_\\mathrm{oxide}$')
ax.set_ylabel('$n_\\mathrm{H2}$')
im = ax.pcolormesh(n_oxide_r, n_H2_r, CeO2, rasterized=True)
ct = ax.contour(n_oxide_r, n_H2_r, CeO2, levels=levels, colors='black')
ax.clabel(ct, fmt='%1.2f')
ax.text(0.68, 0.9, 'CeO$_2$',
        transform = ax.transAxes, fontsize=14)
filename = Path("publication_plots") / "parametric_CeO2.png"
plt.savefig(filename , dpi=400, bbox_inches='tight')
plt.show()

# PLOT 2
fig = plt.figure(figsize=(3.4, 3.3), facecolor='white')
ax = fig.add_subplot()
ax.set_xlabel('$n_\\mathrm{oxide}$')
ax.set_ylabel('$n_\\mathrm{H2}$')
im = ax.pcolormesh(n_oxide_r, n_H2_r, LSF, rasterized=True)
ct = ax.contour(n_oxide_r, n_H2_r, LSF, levels=levels, colors='black')
ax.clabel(ct, fmt='%1.2f')
ax.text(0.75, 0.9, 'LSF',
        transform = ax.transAxes, fontsize=14)
filename = Path("publication_plots") / "parametric_LSF.png"
plt.savefig(filename , dpi=400, bbox_inches='tight')
plt.show()

# PLOT 3 - with colorbar
fig = plt.figure(figsize=(4.1, 3.3), facecolor='white')
ax = fig.add_subplot()
ax.set_xlabel('$n_\\mathrm{oxide}$')
ax.set_ylabel('$n_\\mathrm{H2}$')
im = ax.pcolormesh(n_oxide_r, n_H2_r, CeZr20, rasterized=True)
ct = ax.contour(n_oxide_r, n_H2_r, CeZr20, levels=levels, colors='black')
ax.clabel(ct, fmt='%1.2f')
cb = fig.colorbar(im, ax=ax)
cb.set_label(label="$X_\\mathrm{CO_2}$")
ax.text(0.7, 0.9, 'CeZr20',
        transform = ax.transAxes, fontsize=14)
filename = Path("publication_plots") / "parametric_CeZr20.png"
plt.savefig(filename , dpi=400, bbox_inches='tight')
plt.show()
