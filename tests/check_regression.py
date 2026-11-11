"""Regression checks for the SSME model (run with `make test`).

Reference values are the model's own output after the thrust-equation fix (bug 1),
the mixture molar-mass fix (bug 2), the energy-balance chamber temperature (bug 3),
the Gibbs-energy sign fix (bug 4), the full atom balance in the chamber
equilibrium (bug 5) and the completed H2/F2 model (bug 6).
When a later fix changes the physics on purpose, update REFERENCE here in the same commit.
"""
import contextlib
import io
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from exchanger import (FEE, OCC, OxNozzlePerf, ChamberTemperature, MOLAR_MASS,  # noqa: E402
                       FEE_F, FCC, FluorineChamberTemperature)
from Iterater import (h_table_g, oh_table_g, h2o_table_g, h2_table_g,  # noqa: E402
                      sensible_enthalpy, total_enthalpy, table_interp)

REL_TOL = 1e-6
REFERENCE = {
    'T_FEE': 4617.467833179157,
    'T_chamber': 3948.63358830588,
    'gamma': 1.1994127657131972,
    'M_e': 4.702669411808684,
    'm_dot': 457.9079251531902,
    'F_vac': 2150900.144665112,
    'Isp_vac': 478.820880943043,
    'F_sl': 1728532.8924128348,
    'F2_T_FEE': 7106.693298278559,
    'F2_T_chamber': 5189.731769416725,
    'F2_gamma': 1.3185388840633805,
    'F2_m_dot': 445.6611685885348,
    'F2_Isp_vac': 459.45088712620844,
    'F2_Isp_sl': 362.84214306328806,
    'Isp_sl': 384.7959396613474,
}

p_c = 204 * 101325
phi = (0.166 * MOLAR_MASS['O2'] / MOLAR_MASS['H2']) / 2
R_t, R_e = 5.15, 45.35

with contextlib.redirect_stdout(io.StringIO()):
    T_fee, dh_fee, _, _, h_react = FEE(850, p_c, 530, phi)
    T_c = ChamberTemperature(T_fee, p_c, 0.166, h_react)
    gamma, R, _, _, lma = OCC(T_c, p_c, 0.166)
    m_dot, v_e, p_e, F_vac, c_T, Isp_vac, T_e, M_e = OxNozzlePerf(R_t, R_e, gamma, T_c, R, p_c)
    _, _, _, F_sl, _, Isp_sl, _, _ = OxNozzlePerf(R_t, R_e, gamma, T_c, R, p_c, 101325)
    # H2/F2 at the same equivalence ratio and injection temperatures
    T_fee_F, dh_fee_F, h_react_F = FEE_F(850, p_c, 530, phi)
    T_c_F = FluorineChamberTemperature(T_fee_F, p_c, phi, h_react_F)
    gamma_F, R_F, _, _, lma_F = FCC(T_c_F, p_c, phi)
    m_dot_F, _, _, _, _, Isp_vac_F, _, _ = OxNozzlePerf(R_t, R_e, gamma_F, T_c_F, R_F, p_c)
    _, _, _, _, _, Isp_sl_F, _, _ = OxNozzlePerf(R_t, R_e, gamma_F, T_c_F, R_F, p_c, 101325)

got = dict(T_FEE=T_fee, T_chamber=T_c, gamma=float(gamma), M_e=M_e, m_dot=m_dot,
           F_vac=F_vac, Isp_vac=Isp_vac, F_sl=F_sl, Isp_sl=Isp_sl,
           F2_T_FEE=T_fee_F, F2_T_chamber=T_c_F, F2_gamma=float(gamma_F), F2_m_dot=m_dot_F,
           F2_Isp_vac=Isp_vac_F, F2_Isp_sl=Isp_sl_F)

failures = []
for key, ref in REFERENCE.items():
    ok = math.isclose(got[key], ref, rel_tol=REL_TOL)
    print(f"{'ok  ' if ok else 'FAIL'} {key:10s} {got[key]:.6g} (ref {ref:.6g})")
    if not ok:
        failures.append(key)

# Physics identities that must hold whatever the reference values are
A_t = math.pi * (R_t / 39.37) ** 2
A_e = math.pi * (R_e / 39.37) ** 2
# Equilibrium constants rebuilt from the Gibbs tables (products minus reactants, sigma = +1).
n_h2o, n_h2, n_h, n_oh = (float(x) for x in lma['mols'][0:4])
n_tot = n_h2o + n_h2 + n_h + n_oh
l = MOLAR_MASS['H2'] / (0.166 * MOLAR_MASS['O2'])  # mol O2 per mol H2, as in OCC
RT = 8.314 * T_c
K_h2o = math.exp(-(h_table_g(T_c) + oh_table_g(T_c) - h2o_table_g(T_c)) / RT) / (p_c / 101325)
K_h2 = math.exp(-(2 * h_table_g(T_c) - h2_table_g(T_c)) / RT) / (p_c / 101325)

# Energy balances (per mol O2). Complete combustion: 2 H2O + (2phi - 2) H2 at T_FEE.
h_fee_products = 2 * sensible_enthalpy('H2O', T_fee) + (2 * phi - 2) * sensible_enthalpy('H2', T_fee)
# Equilibrium: the dissociated products at T_c carry exactly the reactant enthalpy.
h_eq_products = sum(n / l * total_enthalpy(sp, T_c)
                    for sp, n in zip(['H2O', 'H2', 'H', 'OH'], [n_h2o, n_h2, n_h, n_oh]))

# H2/F2 chamber: atoms, equilibrium constants and energy balance (per mol F2)
f_hf, f_h2, f_h, f_f = (float(x) for x in lma_F['mols'][0:4])
f_tot = f_hf + f_h2 + f_h + f_f
RT_F = 8.314 * T_c_F
K_h2_F = math.exp(-(2 * h_table_g(T_c_F) - h2_table_g(T_c_F)) / RT_F) / (p_c / 101325)
K_hf_F = math.exp(-(h_table_g(T_c_F) + table_interp('F_table.csv', 'g_hat', T_c_F)
                    - table_interp('HydrogenFluoride_HF_table.csv', 'g_hat', T_c_F)) / RT_F) / (p_c / 101325)
h_fee_F = 2 * sensible_enthalpy('HF', T_fee_F) + (phi - 1) * sensible_enthalpy('H2', T_fee_F)
h_eq_F = sum(n * total_enthalpy(sp, T_c_F) for sp, n in zip(['HF', 'H2', 'H', 'F'], [f_hf, f_h2, f_h, f_f]))

identities = {
    'F2: F atoms conserved (2 in)': math.isclose(f_hf + f_f, 2, rel_tol=1e-9),
    'F2: H atoms conserved (2 phi in)': math.isclose(f_hf + 2*f_h2 + f_h, 2*phi, rel_tol=1e-9),
    'F2: K(H2 -> 2H) matches dG table': math.isclose(f_h**2 / (f_h2 * f_tot), K_h2_F, rel_tol=1e-9),
    'F2: K(HF -> H + F) matches dG table': math.isclose(f_h * f_f / (f_hf * f_tot), K_hf_F, rel_tol=1e-9),
    'F2: FEE energy balance (complete combustion)': abs(h_fee_F - dh_fee_F) < 1.0,
    'F2: chamber energy balance (equilibrium)': abs(h_eq_F - h_react_F) < 1.0,
    'F2: dissociation lowers T (T_c < T_FEE)': T_c_F < T_fee_F,
    'FEE energy balance (complete combustion)': abs(h_fee_products - dh_fee) < 1.0,
    'chamber energy balance (equilibrium)': abs(h_eq_products - h_react) < 1.0,
    'dissociation lowers T (T_c < T_FEE)': T_c < T_fee,
    'enthalpies are zero at 298 K': all(abs(sensible_enthalpy(sp, 298)) < 1e-9
                                        for sp in ['H2', 'O2', 'H2O', 'H', 'OH', 'HF', 'F2', 'F']),
    'K(H2O -> H + OH) matches dG table': math.isclose(n_h * n_oh / (n_tot * n_h2o), K_h2o, rel_tol=1e-9),
    'H atoms conserved (2 in)': math.isclose(2*n_h2o + 2*n_h2 + n_h + n_oh, 2, rel_tol=1e-9),
    'O atoms conserved (2l in)': math.isclose(n_h2o + n_oh, 2*l, rel_tol=1e-9),
    'K(H2 -> 2H) matches dG table': math.isclose(n_h**2 / (n_tot * n_h2), K_h2, rel_tol=1e-9),
    'F_vac = m_dot*v_e + p_e*A_e': math.isclose(F_vac, m_dot * v_e + p_e * A_e, rel_tol=1e-12),
    'F_vac - F_sl = p_atm*A_e': math.isclose(F_vac - F_sl, 101325 * A_e, rel_tol=1e-12),
    'c_T = F/(p_c*A_t)': math.isclose(c_T, F_vac / (p_c * A_t), rel_tol=1e-12),
}
for name, ok in identities.items():
    print(f"{'ok  ' if ok else 'FAIL'} {name}")
    if not ok:
        failures.append(name)

if failures:
    sys.exit(f'{len(failures)} check(s) failed: {failures}')
print('All checks passed.')
