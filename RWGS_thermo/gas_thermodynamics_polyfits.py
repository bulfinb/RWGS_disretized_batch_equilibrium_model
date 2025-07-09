import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

"""
Fits polynomial expressions for:
1. dG_ws: H2O -> H2 + 1/2 O2
2. dG_co: CO2 -> CO + 1/2 O2
3. mu_O2_std: Standard chemical potential of O2
4. dH_ws: Enthalpy change for H2O splitting
5. dH_co: Enthalpy change for CO2 splitting

The resulting polynomials are printed and optionally can be reused in simulation modules.
"""

Rgas = 8.3145  # J/mol·K
dG_ws = []
dG_co = []
dH_ws = []
dH_co = []
mu_O2 = []
Trange = np.arange(573, 1273, 10)
p = 100000  # Pa

# Define Cantera gas objects
O2 = ct.Solution('gri30.yaml')
H2 = ct.Solution('gri30.yaml')
H2O = ct.Solution('gri30.yaml')
CO2 = ct.Solution('gri30.yaml')
CO = ct.Solution('gri30.yaml')

O2.X = {'O2': 1}
H2.X = {'H2': 1}
H2O.X = {'H2O': 1}
CO2.X = {'CO2': 1}
CO.X = {'CO': 1}

for T in Trange:
    for gas in [O2, H2, H2O, CO2, CO]:
        gas.TP = T, p
        gas.equilibrate('TP')

    # Extract H and S [kJ/mol]
    h_O2, s_O2 = O2.enthalpy_mole/1000, O2.entropy_mole/1000
    h_H2, s_H2 = H2.enthalpy_mole/1000, H2.entropy_mole/1000
    h_H2O, s_H2O = H2O.enthalpy_mole/1000, H2O.entropy_mole/1000
    h_CO2, s_CO2 = CO2.enthalpy_mole/1000, CO2.entropy_mole/1000
    h_CO, s_CO = CO.enthalpy_mole/1000, CO.entropy_mole/1000

    # Reaction thermodynamics
    dH_rxn1 = -h_H2O + h_H2 + 0.5 * h_O2
    dS_rxn1 = -s_H2O + s_H2 + 0.5 * s_O2
    dG_rxn1 = dH_rxn1 - T * dS_rxn1

    dH_rxn2 = -h_CO2 + h_CO + 0.5 * h_O2
    dS_rxn2 = -s_CO2 + s_CO + 0.5 * s_O2
    dG_rxn2 = dH_rxn2 - T * dS_rxn2

    mu_O2_std = h_O2 - T * s_O2

    dH_ws.append(dH_rxn1)
    dH_co.append(dH_rxn2)
    dG_ws.append(dG_rxn1)
    dG_co.append(dG_rxn2)
    mu_O2.append(mu_O2_std)

# Fit polynomials (degree 3)
fit_dG_ws = np.polyfit(Trange, dG_ws, 3)
fit_dG_co = np.polyfit(Trange, dG_co, 3)
fit_dH_ws = np.polyfit(Trange, dH_ws, 3)
fit_dH_co = np.polyfit(Trange, dH_co, 3)
fit_mu_O2 = np.polyfit(Trange, mu_O2, 3)

# Print for direct use in simulation code
print("\ndG_WS polynomial coefficients:")
print(fit_dG_ws)
print("\ndG_CO polynomial coefficients:")
print(fit_dG_co)
print("\ndH_WS polynomial coefficients:")
print(fit_dH_ws)
print("\ndH_CO polynomial coefficients:")
print(fit_dH_co)
print("\nmu_O2_std polynomial coefficients:")
print(fit_mu_O2)

# Plot Gibbs free energy changes
plt.figure()
plt.plot(Trange, dG_ws, 'o', label='ΔG H2O splitting')
plt.plot(Trange, dG_co, 'o', label='ΔG CO2 splitting')
plt.plot(Trange, np.polyval(fit_dG_ws, Trange), '-', label='Fit ΔG H2O')
plt.plot(Trange, np.polyval(fit_dG_co, Trange), '-', label='Fit ΔG CO2')
plt.xlabel('Temperature [K]')
plt.ylabel('ΔG [kJ/mol]')
plt.legend()
plt.title('Gibbs Free Energy Change vs Temperature')
plt.grid(True)
plt.show()

# Plot enthalpy of reaction
plt.figure()
plt.plot(Trange, dH_ws, 'o', label='ΔH H2O splitting')
plt.plot(Trange, dH_co, 'o', label='ΔH CO2 splitting')
plt.plot(Trange, np.polyval(fit_dH_ws, Trange), '-', label='Fit ΔH H2O')
plt.plot(Trange, np.polyval(fit_dH_co, Trange), '-', label='Fit ΔH CO2')
plt.xlabel('Temperature [K]')
plt.ylabel('ΔH [kJ/mol]')
plt.legend()
plt.title('Enthalpy of Reaction vs Temperature')
plt.grid(True)
plt.show()

# Plot mu_O2_std
plt.figure()
plt.plot(Trange, mu_O2, 'o', label='μ_O2 std')
plt.plot(Trange, np.polyval(fit_mu_O2, Trange), '-', label='Fit μ_O2')
plt.xlabel('Temperature [K]')
plt.ylabel('μ_O2 [kJ/mol]')
plt.title('Standard Oxygen Chemical Potential')
plt.grid(True)
plt.legend()
plt.show()

# Optionally: define polynomials for import
mu_O2_std_poly = np.poly1d(fit_mu_O2)
dG_WS_poly = np.poly1d(fit_dG_ws)
dG_CO_poly = np.poly1d(fit_dG_co)
dH_WS_poly = np.poly1d(fit_dH_ws)
dH_CO_poly = np.poly1d(fit_dH_co)

