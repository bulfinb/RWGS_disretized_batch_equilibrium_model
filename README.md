---

# 🔥 Countercurrent Chemical Looping RWGS - Discretized Batch Reactors Simulation

This code simulates chemical-looping RWGS in a 1D plug-flow-like reactor with oxides including LSF and CeO2 based oxides.
It uses a discretizion of the oxide (space) and the gas phase (like time discritization) and then solves a sequential series 
of batch equilibriums over the discritised elements. It is useful for quickly simulating thermodynamic limits. 


## 📄 Please Ackowledge

If you make use of this model please cite the published paper describing it [Add citation], and both the materials source data and model (e.g., Bulfin 2016 model, Kuhn 2011 data).

---

## 📘 Project Summary

The simulation computes:
- **Reduction phase**: H₂ + oxide → H₂O + reduced oxide
- **Oxidation phase**: CO₂ + reduced oxide → CO + oxide with countercurrent

It tracks:
- Oxygen non-stoichiometry (`δ`)
- Gas compositions (`x_H2O`, `x_CO2`)
- Mass balances over the whole cycle
- Energy balance over the whole cycle


---

## 🗂 File Overview

| File                            | Purpose                                                                 |
|---------------------------------|-------------------------------------------------------------------------|
| `main.py`                       | Orchestrates redox cycle using the `cycle_until_balanced` function     |
| `discretized_reduction.py`     | Models reduction phase across discretized gas and oxide elements       |
| `discretized_oxidation.py`     | Models oxidation phase across discretized gas and oxide elements       |
| `batch_equilibrium_solver.py`  | Solves gas–solid equilibrium by matching μ_O                           |
| `mass_balance.py`              | Computes mass balance metrics and oxygen transfer                      |
| `energy_balance.py`            | Estimates enthalpy and energy exchange across reactor phases           |
| `build_arrays.py`              | Generates discretized mesh grids for oxide and gas                     |
| `build_mu_O_functions.py`      | Builds interpolating functions for μ_O (solid and gas)                 |
| `gas_thermodynamics.py`        | Gas phase equilibrium properties (μ_O, pO₂) for H₂/H₂O and CO/CO₂      |
| `oxide_thermodynamics.py`      | μ_O equilibrium functions for solid oxides (CeO₂, CeZr, LSF, etc.)     |

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
from main import cycle_until_balanced

results = cycle_until_balanced(
    material="CeO2",
    T=1073,
    x_H2O_0=0.005,
    x_CO2_0=0.995,
    tol=1e-4
)
```

---

## 📊 Output Format

- `delta_red`, `delta_ox`: non-stoichiometry evolution during reduction and oxidation, arrray with shape = (gas_mesh, oxide_mesh)
- `x_H2O_red`, `x_CO2_ox`: gas mole fractions across reactor grid during reduction and oxidation,  arrray with shape = (gas_mesh, oxide_mesh)
- Use `.T[0]` and `.T[-1]` use slices of the mole fraction arrays to inspect outlet behavior for oxidation [0] and reduction [-1]

---

## 🧪 Materials Supported

- `CeO2`, `CeO2_simple`
- `CeZr05`, `CeZr15`, `CeZr20`
- `LSF`

---

## 📐 Thermodynamics

- All μ_O values returned in **kJ/mol**
- μ_O(gas) from polynomial fits of cantera data and equilibrium expressions
- μ_O(solid) from literature-derived enthalpy and entropy models
