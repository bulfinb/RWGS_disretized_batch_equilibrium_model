import numpy as np
import matplotlib.pyplot as plt

from functions.main import simulate_cycle
from functions.mass_balance import calculate_mass_balance
from functions.energy_balance import calculate_energy_balance


import matplotlib as mpl
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color= plt.cm.viridis(np.linspace(0, 1, 11)))


# === Set material from options - "CeO2", "CeZr05", "CeZr15", "CeZr20", "LSF"
material = "CeZr15"

# === Define process conditions
T = 1023  # Temperature in Kelvin
n_CO2 = 1.0  # Moles of CO2 per oxidation cycle
n_H2 = 1.0  # Moles of H2 per reduction cycle
n_oxide = 10  # Moles of oxide material
x_CO2_0 = 0.998  # Initial CO2 mole fraction
x_H2O_0 = 0.005  # Initial H2O mole fraction
gas_mesh, oxide_mesh = 100, 100 # Discretization

# === Simulate 1 cycle
delta_red, xH2O_red, delta_ox, xCO2_ox = simulate_cycle(
    material=material,
    first_cycle=True,
    x_CO2_0=x_CO2_0,
    x_H2O_0=x_H2O_0,
    T=T,
    n_CO2=n_CO2,
    n_H2=n_H2,
    n_oxide=n_oxide,
    gas_mesh = gas_mesh,
    oxide_mesh = oxide_mesh
)

# === CALCULATE ENERGY AND MASS BALANCE
X_CO2, X_H2, O_bal_gas, O_bal_oxide = calculate_mass_balance(
    delta_red, xH2O_red, delta_ox, xCO2_ox,
    n_CO2=n_CO2, n_H2=n_H2, n_oxide=n_oxide,
    x_CO2_0=x_CO2_0, x_H2O_0=x_H2O_0
)
print("MASS BALANCE, X_CO2 =", X_CO2, ", X_H2 =", X_H2, ", O_bal_gas =", O_bal_gas, ", O_bal_oxide =", O_bal_oxide)

# calculate energy balance
Q_red_abs, Q_ox_rel = calculate_energy_balance(
        T, delta_red, material=material, X_CO2=X_CO2)
print("HEAT FLOWS (per mole CO2), Q_red =", Q_red_abs, " kJ/mol", ", Q_ox =", Q_ox_rel, " kJ/mol")


# === PLOTTING ===
def plot_delta_profiles(delta_t_x, text_label, oxide_mesh, gas_mesh):
    x_space = np.arange(oxide_mesh)
    fig, ax = plt.subplots(figsize=(4.0, 3.7))
    ax.set_xlabel('Oxide element # ($x_i$) [-]')
    ax.set_ylabel('Non-stoichiometry $\\delta$ [-]')
    ax.plot(x_space, delta_t_x[0], lw=1.0, label=f'Gas element # 1 ($t_1$)')
    for i in range(10, gas_mesh, max(round(gas_mesh / 10), 1)):
        ax.plot(x_space, delta_t_x[i], lw=1.0, label='__nolegend__')
    ax.plot(x_space, delta_t_x[-1], lw=1.0, label=f'Gas element # 100')
    #ax.legend()
    ax.set_xlim(0, oxide_mesh)
    ax.text(0.45, 0.5, text_label, transform=ax.transAxes)
    ax.text(0.45, 0.5, text_label, transform=ax.transAxes)
    plt.legend()
    plt.show()

def plot_outlet_mole_fraction(x_i_t_x, label, gas_mesh, reduction = True) :
    t_space = np.arange(gas_mesh)
    fig, ax = plt.subplots(figsize=(4.5, 4))
    ax.set_xlabel('Time step [-]')
    ax.set_ylabel(label)
    if reduction:
        ax.plot(t_space, x_i_t_x.T[-1], lw=1.0)
    else:
        ax.plot(t_space, x_i_t_x.T[0], lw=1.0)
    ax.set_xlim(0, gas_mesh)
    ax.set_ylim(0, 1.03)
    plt.legend()
    plt.show()



plot_delta_profiles(delta_red, 'H$_2$ flow $\\rightarrow$', oxide_mesh, gas_mesh)
plot_outlet_mole_fraction(xH2O_red, '$x_{H2O}$ out', gas_mesh)
plot_delta_profiles(delta_ox, '$\\leftarrow$ CO$_2$ flow', oxide_mesh, gas_mesh)
plot_outlet_mole_fraction(xCO2_ox, '$x_{CO_2}$ out', gas_mesh, reduction = False)
