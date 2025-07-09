import numpy as np
from functions.gas_thermodynamics import mu_O2_std

# === CeO2 ACCURATE MODEL ===

def Ce_Dh(d):
    """Partial molar enthalpy from Bulfin 2016 (PCCP)"""
    return 395000.0 - 31400 * np.log10(d)


def Ce_Ds(d):
    """Partial molar entropy from Bulfin 2016 (PCCP)"""
    return 160.5 + 2.94 * 8.3145 * (np.log(0.34 - d) - np.log(d))


def Ce_avgdH(d0, d1):
    """Average ΔH from d0 to d1 (for heat demand analysis)"""
    return (395000.0 * (d1 - d0) - 31400 * (d1 * (np.log10(d1) - 2.302) - d0 * (np.log10(d0) - 2.302))) / (d1 - d0)


def Ce_pO2(d, T):
    """Returns pO2 in bar, with 1 bar being the reference pressure"""
    return np.exp(-2 * Ce_Dh(d) / (8.3145 * T)) * np.exp(2 * Ce_Ds(d) / 8.3145)


def Ce_mu_O(d, T):
    """Returns μ_O in kJ/mol O"""
    return 0.5 * (mu_O2_std(T) + 8.3145 * T * np.log(Ce_pO2(d, T))) / 1000

def Cm_CeO2(T):
    """Molar specific heat capacity of CeO₂ [J/(mol·K)], T in K"""
    return 0.172 * (215.7733 + 0.658257 * T - 6.915805e-4 * T**2 + 3.30995e-7 * T**3 - 5.593514e-11 * T**4)


# === Ce₀.₉₅Zr₀.₀₅O₂ MODEL ===
# === Base data for all CeZr materials is from YANG et al. 2014 Ceria−Zirconia Solid Solutions dx.doi.org/10.1021/cm503131p

def CeZr05_Dh(d):
    """Partial molar enthalpy from Bulfin 2016 (PCCP).
    Independent of delta but we keep delta for consistance with other CeO2 materials"""
    return 395000.0


def CeZr05_Ds(d):
    """Partial molar entropy from Bulfin 2016 (PCCP)"""
    return 144.0 + 2.44 * 8.3145 * (np.log(0.41 - d) - np.log(d))


def CeZr05_pO2(d, T):
    """Returns pO2 in bar, with 1 bar being the reference pressure"""
    return np.exp(-2 * CeZr05_Dh(d) / (8.3145 * T)) * np.exp(2 * CeZr05_Ds(d) / 8.3145)


def CeZr05_mu_O(d, T):
    """Returns μ_O in kJ/mol O"""
    return 0.5 * (mu_O2_std(T) + 8.3145 * T * np.log(CeZr05_pO2(d, T))) / 1000

# === Ce₀.₈₅Zr₀.₁₅O₂ MODEL ===

def CeZr15_Dh(d):
    """Partial molar enthalpy from Bulfin 2016 (PCCP)"""
    return 392000.0 + 18220 * np.log10(d)


def CeZr15_Ds(d):
    """Partial molar entropy from Bulfin 2016 (PCCP)"""
    return 137.5 + 2.353 * 8.3145 * (np.log(0.425 - d) - np.log(d))


def CeZr15_pO2(d, T):
    """Returns pO2 in bar, with 1 bar being the reference pressure"""
    return np.exp(-2 * CeZr15_Dh(d) / (8.3145 * T)) * np.exp(2 * CeZr15_Ds(d) / 8.3145)


def CeZr15_mu_O(d, T):
    """Returns μ_O in kJ/mol O"""
    return 0.5 * (mu_O2_std(T) + 8.3145 * T * np.log(CeZr15_pO2(d, T))) / 1000



# === Ce0.8Zr0.2O₂ MODEL ===
def CeZr20_Dh(d):
    """Partial molar enthalpy from Yang et al. 2014 Cer-Zirconia"""
    return 392000.0 + 27220 * np.log10(d)


def CeZr20_Ds(d):
    """Partial molar entropy from Bulfin 2016 (PCCP) but d_max =0.4 and n=1/0.4 modified to fit Zr20 data"""
    return 135.5 + 2.5 * 8.3145 * (np.log(0.4 - d) - np.log(d))


def CeZr20_pO2(d, T):
    """Returns pO2 in bar, with 1 bar being the reference pressure"""
    return np.exp(-2 * CeZr20_Dh(d) / (8.3145 * T)) * np.exp(2 * CeZr20_Ds(d) / 8.3145)

def CeZr20_mu_O(d, T):
    """Returns μ_O in kJ/mol O"""
    return 0.5 * (mu_O2_std(T) + 8.3145 * T * np.log(CeZr20_pO2(d, T))) / 1000


# === CeZr LOW-H MODEL ===

def CeZr_lowH_Dh(d):
    """Made up material with a lower than CeZr materials normally have"""
    return 300000.0


def CeZr_lowH_pO2(d, T):
    """Returns pO2 in bar, with 1 bar being the reference pressure"""
    return np.exp(-2 * CeZr_lowH_Dh(d) / (8.3145 * T)) * np.exp(2 * CeZr15_Ds(d) / 8.3145)


def CeZr_lowH_mu_O(d, T):
    """Returns μ_O in kJ/mol O"""
    return 0.5 * (mu_O2_std(T) + 8.3145 * T * np.log(CeZr_lowH_pO2(d, T))) / 1000


# === CeO₂ SIMPLE MODEL ===

def Ce_simple_Dh(d):
    """Partial molar enthalpy from Bulfin 2016 (PCCP) - simple model"""
    return 430000.0


def Ce_simple_Ds(d):
    """Partial molar entropy from Bulfin 2016 (PCCP) - simple model"""
    return 165 + 2.31 * 8.3145 * (np.log(0.35 - d) - np.log(d))


def Ce_simple_pO2(d, T):
    """Returns pO2 in bar, with 1 bar being the reference pressure"""
    return np.exp(-2 * Ce_simple_Dh(d) / (8.3145 * T)) * np.exp(2 * Ce_simple_Ds(d) / 8.3145)


def Ce_simple_mu_O(d, T):
    """Returns μ_O in kJ/mol O"""
    return 0.5 * (mu_O2_std(T) + 8.3145 * T * np.log(Ce_simple_pO2(d, T))) / 1000


# === La0.6Sr0.4FeO3 − δ (LSF) MODEL (M. Kuhn et al. Solid State Ionics 195 (2011) 7–15) ===

def LSF_K_i(T):
    """K_i from Kuhn et al. 2011 (eq. 15, table 1)"""
    return np.exp(-95750 / (8.3145 * T) - 21.63 / 8.3145)


def LSF_K_ox(T):
    """K_ox from Kuhn et al. 2011 (eq. 15, table 1)"""
    return np.exp(95620 / (8.3145 * T) - 54.27 / 8.3145)


def LSF_pO2(delta, T):
    """pO₂ from Kuhn et al. 2011 (eq. 9)"""
    # Term 1
    K_ox = LSF_K_ox(T)
    K_i = LSF_K_i(T)
    term1 = -(1 / np.sqrt(K_ox)) * (2 * delta - 0.4) * np.sqrt(3 - delta) / ((2 * delta + 0.6) * np.sqrt(delta))

    # Term 2
    term2a = (1 / K_ox) * ((3 - delta) * (2 * delta - 0.4)**2 / (delta * (2 * delta + 0.6)**2))
    term2b = (K_i  / K_ox) * (4 * (3 - delta) * (1.4 - 2 * delta) / (delta * (2 * delta + 0.6)))
    inner = term1 + np.sqrt(term2a + term2b)

    return (1 / 16) * inner**4

def LSF_Dh(delta, T):
    """0.5 R ln(pO2) = \Delta H(\delta, T)/ T + \Delta S(\delta, T)
    Use slope of this equation with dT = 1 to get dH"""
    dy = -0.5 * 8.3145 * (np.log(LSF_pO2(delta,T+1))-np.log(LSF_pO2(delta,T)))
    dx = 1/(T+1) - 1/T
    dH = dy / dx
    return dH


def LSF_mu_O(d, T):
    """Returns μ_O in kJ/mol O"""
    return 0.5 * (mu_O2_std(T) + 8.3145 * T * np.log(LSF_pO2(d, T))) / 1000

def Debye_heat_capacity(T_debye, T):
        """
        Calculates the heat capacity at a given temperature from the Debye Temperature
        :param T_debye:          debye temperatute
        :param T:        temperature
        :return:            heat capacity per atom
        """
        t_ratio = T / T_debye

        def integrand(x):
            return (x ** 4 * np.exp(x)) / (np.exp(x) - 1) ** 2

        if isinstance(t_ratio, int) or isinstance(t_ratio, float):
            cv = 9 * 8.314 * (t_ratio ** 3) * quad(integrand, 0, t_ratio ** -1)[0]

        else:
            cv = []
            for i in range(len(t_ratio)):
                cv_i = 9 * 8.314 * (t_ratio[i] ** 3) * quad(integrand, 0, t_ratio[i] ** -1)[
                    0]  # quad is numerical integrate
                cv = np.append(cv, cv_i)
        return cv


def Cm_LSF(T, delta):
    """
    Calculates the heat capacity of SrFeO_3-delta at a given temperature. Uses the Debye model with the temperature calculated from materials project SrFeO3
    To include the delta dependence we use Dulong petit law, i.e. Cm \propto (5-delta) which is the number of atoms per molecule.
    :param T:          Temperature [K]
    :param delta:        Non stoicheometry delta
    :return:            heat capacity per mole of SrFeO3-delta
    """
    T_debye_LSF = 572.22
    Cm_result = (5 - delta) * Debye_heat_capacity(T_debye_SrFeO3, T)
    return Cm_result