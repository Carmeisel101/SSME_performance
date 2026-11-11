"""Regression checks for the SSME model (run with `make test`).

Reference values are the model's own output after the thrust-equation fix (bug 1),
the mixture molar-mass fix (bug 2), the Gibbs-energy sign fix (bug 4) and the
full atom balance in the chamber equilibrium (bug 5).
When a later fix changes the physics on purpose, update REFERENCE here in the same commit.
"""
import contextlib
import io
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from exchanger import FEE, OCC, OxNozzlePerf  # noqa: E402
from Iterater import h_table_g, oh_table_g, h2o_table_g, h2_table_g  # noqa: E402

REL_TOL = 1e-6
REFERENCE = {
    'T_chamber': 3907.4016344537067,
    'gamma': 1.1984035316417352,
    'M_e': 4.695198983375653,
    'm_dot': 461.81693456133354,
    'F_vac': 2152338.8160133893,
    'Isp_vac': 475.08550086708595,
    'F_sl': 1729971.5637611125,
    'Isp_sl': 381.85642554994064,
}

p_c = 204 * 101325
phi = (0.166 * 31.998 / 2.016) / 2
R_t, R_e = 5.15, 45.35

with contextlib.redirect_stdout(io.StringIO()):
    T_c, *_ = FEE(850, p_c, 530, phi)
    gamma, R, _, _, lma = OCC(T_c, p_c, 0.166)
    m_dot, v_e, p_e, F_vac, c_T, Isp_vac, T_e, M_e = OxNozzlePerf(R_t, R_e, gamma, T_c, R, p_c)
    _, _, _, F_sl, _, Isp_sl, _, _ = OxNozzlePerf(R_t, R_e, gamma, T_c, R, p_c, 101325)

got = dict(T_chamber=T_c, gamma=float(gamma), M_e=M_e, m_dot=m_dot,
           F_vac=F_vac, Isp_vac=Isp_vac, F_sl=F_sl, Isp_sl=Isp_sl)

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
l = (2.02 / 0.166) / (2 * 15.99)  # mol O2 per mol H2, as in OCC
RT = 8.314 * T_c
K_h2o = math.exp(-(h_table_g(T_c) + oh_table_g(T_c) - h2o_table_g(T_c)) / RT) / (p_c / 101325)
K_h2 = math.exp(-(2 * h_table_g(T_c) - h2_table_g(T_c)) / RT) / (p_c / 101325)

identities = {
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
