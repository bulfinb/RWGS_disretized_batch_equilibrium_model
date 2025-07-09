from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os

# Add project root to sys.path
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

from functions.mass_balance import calculate_mass_balance
from functions.main import cycle_until_balanced

# === CONFIGURATION ===
material = "CeO2"
T = 1073                # Temperature [K]
p = 100000              # Pressure [Pa]
n_CO2 = 1.0             # Moles of CO2 per cycle
n_H2 = 1.5              # Moles of H2 per cycle
n_oxide = 40            # Moles of oxide in the bed
x_CO2_0 = 0.998         # Initial mole fraction CO2
x_H2O_0 = 0.002         # Initial mole fraction H2O
mesh_range = np.arange(2, 250, 2)

# === Initialize storage lists ===
solve_time = []
X_CO2_values = []
X_H2_values = []
O_bal_gas_values = []
O_bal_oxide_values = []

# === Mesh Sweep Simulation ===
for mesh in mesh_range:
    tic = time.perf_counter()
    delta_t_x_red, x_H2O_t_x_red, delta_t_x_ox, x_CO2_t_x_ox, cycles = cycle_until_balanced(
        material=material,
        T=T,
        n_H2=n_H2,
        n_oxide=n_oxide,
        x_CO2_0=x_CO2_0,
        x_H2O_0=x_H2O_0,
        oxide_mesh=mesh,
        gas_mesh=mesh
    )

    # Compute mass balance
    X_CO2, X_H2, O_bal_gas, O_bal_oxide = calculate_mass_balance(
        delta_t_x_red, x_H2O_t_x_red,
        delta_t_x_ox, x_CO2_t_x_ox,
        n_H2=n_H2, n_oxide=n_oxide,
        x_CO2_0=x_CO2_0, x_H2O_0=x_H2O_0
    )

    X_CO2_values.append(X_CO2)
    X_H2_values.append(X_H2)
    O_bal_gas_values.append(O_bal_gas)
    O_bal_oxide_values.append(O_bal_oxide)

    toc = time.perf_counter()
    solve_time.append(toc-tic)
    print(f"mesh: {mesh}, cycles: {cycles}, X_H2: {X_H2:.3f}, X_CO2: {X_CO2:.3f}, "
          f"MB_gas: {O_bal_gas:.3e}, MB_oxide: {O_bal_oxide:.3e}, t: {toc - tic:.2f}s")

# === Plot Results ===
os.makedirs("plots", exist_ok=True)

plt.figure(figsize=(3.4, 3.3), facecolor='white')
plt.plot(mesh_range, X_CO2_values, label="$X_{CO2}$", color = 'black')
plt.plot(mesh_range, X_H2_values, label="$X_{H2}$", color = 'grey')
plt.xlabel("$N_\\mathrm{oxide}$, $N_\\mathrm{gas}$  [-]")
plt.ylabel("Conversion extent $X_i$ [-]")
plt.axvline(x = 100, color = 'black', ls = '--', label = 'default value')
plt.legend()

plt.set_xlim(0, 250)
filename = Path("plots") / "mesh_test_conversion.png"
plt.savefig(filename, dpi=300, bbox_inches='tight')

plt.figure()
plt.plot(mesh_range, solve_time, label="$X_{CO2}$")
plt.xlabel("Number of mesh elements [-]")
plt.ylabel("Solve time [-]")
filename = Path("plots") / "mesh_test_solve_time.png"
plt.savefig(filename, dpi=300, bbox_inches='tight')

plt.figure()
plt.plot(mesh_range, O_bal_gas_values, label="O-balance CO/H$_2$O")
plt.plot(mesh_range, O_bal_oxide_values, label="O-balance oxide/CO$_2$")
plt.xlabel("Mesh Size")
plt.ylabel("Oxygen Balance Error []")
plt.legend()
filename = Path("plots") / "mesh_test_mass_balance.png"
plt.savefig(filename, dpi=300, bbox_inches='tight')
plt.show()

