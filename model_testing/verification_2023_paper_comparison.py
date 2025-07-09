from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import os

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color= plt.cm.viridis(np.linspace(0, 1, 11)))

# Add project root to sys.path
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

from functions.mass_balance import calculate_mass_balance
from functions.main import simulate_cycle

def calculate_moles(u, t, A, p, T):
    """Calculates moles of gases from flow rates and reduction/oxidation times"""
    R = 8.3145  # J/mol/K
    return u * t * A * p / (R * T)

def plot_profiles(delta_t_x, text_label, oxide_mesh, gas_mesh, filename):
    x_space = np.arange(oxide_mesh)
    fig, ax = plt.subplots(figsize=(4.5, 4))
    ax.set_xlabel('Oxide element # ($x_i$) [-]')
    ax.set_ylabel('Non-stoichiometry $\\delta$ [-]')
    ax.plot(x_space, delta_t_x[0], lw=1.0, label=f'Gas element # 1 ($t_1$)')
    for i in range(10, gas_mesh, max(round(gas_mesh / 10), 1)):
        ax.plot(x_space, delta_t_x[i], lw=1.0, label='__nolegend__')
    ax.plot(x_space, delta_t_x[-1], lw=1.0, label=f'Gas element # 100')
    #ax.legend()
    ax.set_xlim(0, oxide_mesh)
    ax.set_ylim(0, 0.08)
    ax.text(0.45, 0.5, text_label, transform=ax.transAxes)
    ax.text(0.45, 0.5, text_label, transform=ax.transAxes)
    plt.legend()
    filename = Path("plots") / filename
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

def plot_outlet_mole_fraction(data_t_x, label, gas_mesh, filename):
    t_space = np.arange(gas_mesh)
    fig, ax = plt.subplots(figsize=(4.5, 4))
    ax.set_xlabel('Time step [-]')
    ax.set_ylabel(label)
    ax.plot(t_space, data_t_x.T[-1], lw=1.0)
    filename = Path("plots") / filename
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

# === CONFIGURATION === From Bulfin et al. 2023 Intensification....
material = "CeO2_simple"
T = 1073  # Temperature in Kelvin
p = 1e5   # Pressure in Pascals (1 bar)
A_reactor = 1  # Reactor cross-sectional area in m^2
n_oxide = 10000  # Total moles of oxide in the reactor
u = 0.05  # Gas flow velocity in m/s

# Times from ESI
t_red1, t_ox1 = 521, 436.6 # s
t_red2, t_ox2 = 446.7, 404.5 # s

# Input molefractions
x_CO2_0, x_H2O_0 = 0.999, 0.001

# Mesh
gas_mesh, oxide_mesh = 100, 100

# === CYCLE 1 ===
n_H2_cycle_1 = calculate_moles(u, t_red1, A_reactor, p, T)
n_CO2_cycle_1 = calculate_moles(u, t_ox1, A_reactor, p, T)
print("n_CO2_1 = ", n_CO2_cycle_1, "n_H2_1 = ", n_H2_cycle_1)


delta_red_1, xH2O_red_1, delta_ox_1, xCO2_ox_1 = simulate_cycle(
    material=material,
    first_cycle=True,
    x_CO2_0=x_CO2_0,
    x_H2O_0=x_H2O_0,
    T=T,
    n_CO2=n_CO2_cycle_1,
    n_H2=n_H2_cycle_1,
    n_oxide=n_oxide,
    gas_mesh = gas_mesh,
    oxide_mesh = oxide_mesh
)

X_CO2, X_H2, O_bal_gas, O_bal_oxide = calculate_mass_balance(
    delta_red_1, xH2O_red_1, delta_ox_1, xCO2_ox_1,
    n_CO2=n_CO2_cycle_1, n_H2=n_H2_cycle_1, n_oxide=n_oxide,
    x_CO2_0=x_CO2_0, x_H2O_0=x_H2O_0
)
print("Cycle 1: X_CO2 =", X_CO2, ", X_H2 =", X_H2, ", O_bal_gas =", O_bal_gas, ", O_bal_oxide =", O_bal_oxide)

delta_x_0 = delta_ox_1[-1]

# === CYCLE # 2 ===
n_H2_cycle_2 = calculate_moles(u, t_red2, A_reactor, p, T)
n_CO2_cycle_2 = calculate_moles(u, t_ox2, A_reactor, p, T)
print("n_CO2_2 = ", n_CO2_cycle_2, "n_H2_2 = ", n_H2_cycle_2)

delta_red_2, xH2O_red_2, delta_ox_2, xCO2_ox_2 = simulate_cycle(
    material=material,
    first_cycle=False,
    delta_x_0=delta_x_0,
    x_CO2_0=x_CO2_0,
    x_H2O_0=x_H2O_0,
    T=T,
    n_CO2=n_CO2_cycle_2,
    n_H2=n_H2_cycle_2,
    n_oxide=n_oxide,
    gas_mesh=gas_mesh,
    oxide_mesh=oxide_mesh
)

X_CO2, X_H2, O_bal_gas, O_bal_oxide = calculate_mass_balance(
    delta_red_2, xH2O_red_2, delta_ox_2, xCO2_ox_2,
    n_CO2=n_CO2_cycle_2, n_H2=n_H2_cycle_2, n_oxide=n_oxide,
    x_CO2_0=x_CO2_0, x_H2O_0=x_H2O_0
)
print("Cycle 2: X_CO2 =", X_CO2, ", X_H2 =", X_H2, ", O_bal_gas =", O_bal_gas, ", O_bal_oxide =", O_bal_oxide)


# === PLOTTING ===
os.makedirs("plots", exist_ok=True)
plot_profiles(delta_red_1, 'H$_2$ flow $\\rightarrow$', oxide_mesh, gas_mesh, "2023_paper_model_reduction_delta_x_t_cycle1.png")
plot_outlet_mole_fraction(xH2O_red_1, '$X_{H2O}$ out', gas_mesh, "2023_paper_model_reduction_X_H2O_outflow_cycle1.png")
plot_profiles(delta_ox_1, '$\\leftarrow$ CO$_2$ flow', oxide_mesh, gas_mesh, "2023_paper_model_oxidation_delta_x_t_cycle1.png")
plot_profiles(delta_red_2, 'H$_2$ flow $\\rightarrow$', oxide_mesh, gas_mesh, "2023_paper_model_reduction_delta_x_t_cycle2.png")
plot_outlet_mole_fraction(xH2O_red_2, '$X_{H2O}$ out', gas_mesh,  "2023_paper_model_reduction_X_H2O_outflow_cycle2.png")
plot_profiles(delta_ox_2, '$\\leftarrow$ CO$_2$ flow', oxide_mesh, gas_mesh,  "2023_paper_model_oxidation_delta_x_t_cycle2.png")


