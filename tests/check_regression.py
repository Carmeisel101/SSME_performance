"""Regression checks for the SSME model (run with `make test`).

Reference values are the model's own output after the thrust-equation fix (bug 1),
the mixture molar-mass fix (bug 2) and the Gibbs-energy sign fix (bug 4).
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
    'gamma': 1.1954325736053695,
    'M_e': 4.673298555622012,
    'm_dot': 470.1795573822487,
    'F_vac': 2156600.0834952,
    'Isp_vac': 467.55948058757104,
    'F_sl': 1734232.8312429232,
    'Isp_sl': 375.98857943086904,
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
# The current species model normalizes by N_H2O + N_H2 (see bug 5), so the check uses that too.
n_h2o, n_h2, n_h, n_oh = (float(x) for x in lma['mols'][0:4])
RT = 8.314 * T_c
K_h2o = math.exp(-(h_table_g(T_c) + oh_table_g(T_c) - h2o_table_g(T_c)) / RT) / (p_c / 101325)
K_h2 = math.exp(-(2 * h_table_g(T_c) - h2_table_g(T_c)) / RT) / (p_c / 101325)

identities = {
    'K(H2O -> H + OH) matches dG table': math.isclose(n_h * n_oh / ((n_h2o + n_h2) * n_h2o), K_h2o, rel_tol=1e-9),
    'K(H2 -> 2H) matches dG table': math.isclose(n_h**2 / ((n_h2o + n_h2) * n_h2), K_h2, rel_tol=1e-9),
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
