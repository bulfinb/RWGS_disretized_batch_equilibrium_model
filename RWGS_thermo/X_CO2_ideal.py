import numpy as np


"""In this file we calculate the thermodynamic limit on CO2 conversion extent for the reverse 
water-gas shift process performed in a countercurrent oxygen permeable membrane reactor. 
The calculation follows the methodology developed in previous work DOI:10.1039/C8CP07077F"""

# From gas thermodynamics polyfits
dG_WS = np.poly1d([2.70345921e-06, -1.13933543e-02, -4.06918554e+01,  2.41896010e+05])  # H2 + 0.5O2 = H2O
dG_CS = np.poly1d([2.32796741e-07,  1.06605461e-03, -8.99506394e+01,  2.84170742e+05])  # CO + 0.5O2 = CO2

def X_CO2_membrane(T, nH2):
    """Determine the maximum conversion of CO2 for countercurrent flows of CO2 and H2 with
     flow rates of F_CO2 = 1, F_H2 = nH2, at temperature T [K], and presure p [Pa].
     kappa is the exchange coefficient which corresponds to the amount of exchanged O2.
     kappa = 0 means no reaction, kappa = 0.5 complete reaction.
     Paramaters
     -T   (temperature [K])
     -p   (pressure [Pa])
     -nH2 (hydrogen excess n_H2/n_CO2 [-])
     returns
     -X_CO2 (conversion extent of CO2 [-])"""
    # Make arrays of the pO2(kappa) in each flow (note \mu_O2(p_O2))
    pO2_flow2 = []
    pO2_flow1 = []
    # complete conversion of CO2 to CO gives kappa = 0.5, kappa in range 0-0.5
    kappa_range = np.arange(0, 0.502, 0.001)
    for kappa in kappa_range:
        X_CO2 = 0.9999 - 2 * kappa
        X_H2O = 2 * kappa / nH2
        pO2_f1 = 1 * np.exp(-2 * dG_CS(T) / (8.3145 * T)) * (X_CO2 / (1 - X_CO2)) ** 2
        pO2_f2 = 1 * np.exp(-2 * dG_WS(T) / (8.3145 * T)) * ( X_H2O / (1 - X_H2O)) ** 2
        # add the new pO2 values to the arrays
        pO2_flow1.append(pO2_f1)
        # for the second flow add pO2 value to the start of array to reverse kappa for countercurrent
        pO2_flow2.insert(0, pO2_f2)
        # if the arrays meet at any point, then the chemical potentials of oxygen are equal at that point
        # and we have reached the maximum exchange extent kappa_max
        if np.any(np.asarray(pO2_flow1) - np.asarray(pO2_flow2) <= 0) or kappa > 0.499:
            X_CO2 = 2 * kappa
            break
    return X_CO2

def X_CO2_catalytic(T, nH2):
    """Determine the maximum conversion of CO2 for countercurrent flows of CO2 and H2 with
     flow rates of F_CO2 = 1, F_H2 = nH2, at temperature T [K], and presure p [Pa].
     kappa is the exchange coefficient which corresponds to the amount of exchanged O2.
     kappa = 0 means no reaction, kappa = 0.5 complete reaction.
     Paramaters
     -T   (temperature [K])
     -p   (pressure [Pa])
     -nH2 (hydrogen excess n_H2/n_CO2 [-])
     returns
     -X_CO2 (conversion extent of CO2 [-])"""
    # Make arrays of the pO2(kappa) in each flow (note \mu_O2(p_O2))
    pO2_flow2 = []
    pO2_flow1 = []
    # complete conversion of CO2 to CO gives kappa = 0.5, kappa in range 0-0.5
    kappa_range = np.arange(0, 0.502, 0.001)
    for kappa in kappa_range:
        X_CO2 = 0.9999 - 2 * kappa
        X_H2O = 2 * kappa / nH2
        pO2_f1 = 1 * np.exp(-2 * dG_CS(T) / (8.3145 * T)) * (X_CO2 / (1 - X_CO2)) ** 2
        pO2_f2 = 1 * np.exp(-2 * dG_WS(T) / (8.3145 * T)) * ( X_H2O / (1 - X_H2O)) ** 2
        # add the new pO2 values to the arrays
        pO2_flow1.append(pO2_f1)
        # for the second flow add pO2 value to the start of array to reverse kappa for countercurrent
        pO2_flow2.append(pO2_f2)
        # if the arrays meet at any point, then the chemical potentials of oxygen are equal at that point
        # and we have reached the maximum exchange extent kappa_max
        if np.any(np.asarray(pO2_flow1) - np.asarray(pO2_flow2) <= 0) or kappa > 0.499:
            X_CO2 = 2 * kappa
            break
    return X_CO2