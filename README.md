# 🔥 Discretized Redox Reactor Simulation

This code simulates chemical-looping RWGS in a 1D plug-flow-like reactor with oxides including LSF and CeO2 based oxides.
It uses a discretizion of the oxide (space) and the gas phase (like time discritization) and then solves a sequential series 
of batch equilibriums over the discritised elements. It is useful for quickly simulating thermodynamic limits. 

---

## 📘 Project Summary

The simulation computes:
- **Reduction phase**: H₂ + oxide → H₂O + reduced oxide
- **Oxidation phase**: CO₂ + reduced oxide → CO + oxide with countercurrent

It tracks:
- Oxygen stoichiometry (`δ`)
- Gas compositions (`x_H2O`, `x_CO2`)
- Mass balances
- Oxygen chemical potentials in both solid and gas phases

---

## 🗂 File Overview

| File | Purpose                                                                              |
|------|--------------------------------------------------------------------------------------|
| `main.py` | Orchestrates one full chemical-looping cycle over the discretized oxide and gas grid |
| `discretized_reduction.py` | Models reduction phase across time and space                                         |
| `discretized_oxidation.py` | Models oxidation phase across time and space                                         |
| `batch_equilibrium_solver.py` | Solves for chemical equilibrium between solid and gas via μ_O matching               |
| `mass_balance.py` | Computes conversion and oxygen mass balance metrics                                  |
| `build_arrays.py` | Generates discretized δ and gas composition arrays                                   |
| `build_mu_O_functions.py` | Builds interpolation functions for μ_O (solid & gas)                                 |
| `gas_thermodynamics.py` | Provides μ_O and pO₂ equilibrium functions for H₂/H₂O and CO/CO₂                     |
| `oxide_thermodynamics.py` | μ_O equilibrium models for multiple oxide materials (CeO₂, CeZr, LSF, etc.)          |

---

## ⚙️ How It Works

1. **Set material and operating conditions** in `simulate_cycle()`
2. **Discretize** reactor space (oxide_mesh) and gas flow (gas_mesh)
3. **Run reduction** using `compute_reduction_step()`
4. **Feed reduction output** to `compute_oxidation_step()`
5. **Analyze results** using `calculate_mass_balance()` or visualization routines

---

## 🚀 Quick Example

```python
from main import simulate_cycle

delta_red, x_H2O_red, delta_ox, x_CO2_ox = simulate_cycle(
    material="CeO2",
    T=1073,
    x_H2O_0=0.005,
    x_CO2_0=0.995
)
```

---

## 📊 Output Format

- `delta_red`, `delta_ox`: shape = (gas_mesh, oxide_mesh)
- `x_H2O_red`, `x_CO2_ox`: gas mole fractions across reactor grid
- Use `.T[0]` and `.T[-1]` slices to inspect outlet behavior

---

## 🧪 Materials Supported

- `CeO2`, `CeO2_simple`
- `CeZr05`, `CeZr15`
- `LSF`
- `CeZr_lowH` (hypothetical low-enthalpy ceria based oxide)

---

## 📐 Thermodynamics

- All μ_O values returned in **kJ/mol**
- μ_O(gas) from polynomial fits of cantera data and equilibrium expressions
- μ_O(solid) from literature-derived enthalpy and entropy models

---

## 📄 License and Attribution

If using this framework in publications or reports, please cite the source model/data (e.g., Bulfin 2016, Kuhn 2011) and acknowledge the simulation code if applicable.
