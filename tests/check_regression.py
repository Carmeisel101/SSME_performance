"""Regression checks for the SSME model (run with `make test`).

Reference values are the model's own output after the thrust-equation fix.
When a later fix changes the physics on purpose, update REFERENCE here in the same commit.
"""
import contextlib
import io
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from exchanger import FEE, OCC, OxNozzlePerf  # noqa: E402

REL_TOL = 1e-6
REFERENCE = {
    'T_chamber': 3907.4016344537067,
    'gamma': 1.16215449646827,
    'M_e': 4.437020690819335,
    'm_dot': 509.436188232153,
    'F_vac': 2207117.6879859576,
    'Isp_vac': 441.63826154075105,
    'F_sl': 1784750.4357336806,
    'Isp_sl': 357.1237202311504,
}

p_c = 204 * 101325
phi = (0.166 * 31.998 / 2.016) / 2
R_t, R_e = 5.15, 45.35

with contextlib.redirect_stdout(io.StringIO()):
    T_c, *_ = FEE(850, p_c, 530, phi)
    gamma, R, *_ = OCC(T_c, p_c, 0.166)
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
identities = {
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
