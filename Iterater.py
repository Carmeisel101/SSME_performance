from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import fsolve

# Thermo tables live next to this file, so the code runs from any working directory
TABLES = Path(__file__).resolve().parent / 'tables'


def h2_tableH(Temp):
    '''
    :param Temp: temperature in K
    :return h_hat: enthalpy in J/mol also known as heat of formation
    '''
    h2_table = pd.read_csv(TABLES / 'Hydrogen_H2_table.csv')

    idx = (h2_table['T']-Temp).abs().idxmin()
    T_cloest = h2_table['T'][idx]
    h_hat_close = h2_table['H*-H*_298'][idx]
    if Temp > T_cloest:
        T_upper = h2_table['T'][idx+1]
        h_hat_upper = h2_table['H*-H*_298'][idx+1]
        h_hat = h_hat_close + (h_hat_upper-h_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = h2_table['T'][idx-1]
        h_hat_lower = h2_table['H*-H*_298'][idx-1]
        h_hat = h_hat_close + (h_hat_lower-h_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)
    h_hat = h_hat*4184

    return h_hat

def f2_tableH(Temp):
    '''
    :param Temp: temperature in K
    :return h_hat: enthalpy in J/mol also known as heat of formation
    '''
    f2_table = pd.read_csv(TABLES / 'Fluorine_F2_table.csv')

    idx = (f2_table['T']-Temp).abs().idxmin()
    T_cloest = f2_table['T'][idx]
    h_hat_close = f2_table['H*-H*_298'][idx]
    if Temp > T_cloest:
        T_upper = f2_table['T'][idx+1]
        h_hat_upper = f2_table['H*-H*_298'][idx+1]
        h_hat = h_hat_close + (h_hat_upper-h_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = f2_table['T'][idx-1]
        h_hat_lower = f2_table['H*-H*_298'][idx-1]
        h_hat = h_hat_close + (h_hat_lower-h_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)
    h_hat = h_hat*4184

    return h_hat

def o2_tableH(Temp):
    '''
    :param Temp: temperature in K
    :return h_hat: enthalpy in J/mol also known as heat of formation
    '''
    o2_table = pd.read_csv(TABLES / 'o2_table.csv')

    idx = (o2_table['T']-Temp).abs().idxmin()
    T_cloest = o2_table['T'][idx]
    h_hat_close = o2_table['H*-H*_298'][idx]
    if Temp > T_cloest:
        T_upper = o2_table['T'][idx+1]
        h_hat_upper = o2_table['H*-H*_298'][idx+1]
        h_hat = h_hat_close + (h_hat_upper-h_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = o2_table['T'][idx-1]
        h_hat_lower = o2_table['H*-H*_298'][idx-1]
        h_hat = h_hat_close + (h_hat_lower-h_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)
    h_hat = h_hat

    return h_hat

def hf_tableH(Temp):
    '''
    :param Temp: temperature in K
    :return h_hat: enthalpy in J/mol also known as heat of formation
    '''
    hf_table = pd.read_csv(TABLES / 'HydrogenFluoride_HF_table.csv')

    idx = (hf_table['T']-Temp).abs().idxmin()
    T_cloest = hf_table['T'][idx]
    h_hat_close = hf_table['h_hat'][idx]
    if Temp > T_cloest:
        T_upper = hf_table['T'][idx+1]
        h_hat_upper = hf_table['h_hat'][idx+1]
        h_hat = h_hat_close + (h_hat_upper-h_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = hf_table['T'][idx-1]
        h_hat_lower = hf_table['h_hat'][idx-1]
        h_hat = h_hat_close + (h_hat_lower-h_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)
    h_hat = h_hat*4184

    return h_hat

def h2o_table_f(Temp):
    '''
    :param Temp: temperature in K
    :return h_hat: enthalpy in J/mol also known as heat of formation
    '''
    h2o_table = pd.read_csv(TABLES / 'H2O_table.csv')

    idx = (h2o_table['T']-Temp).abs().idxmin()
    T_cloest = h2o_table['T'][idx]
    h_hat_close = h2o_table['h_hat'][idx]
    if Temp > T_cloest:
        T_upper = h2o_table['T'][idx+1]
        h_hat_upper = h2o_table['h_hat'][idx+1]
        h_hat = h_hat_close + (h_hat_upper-h_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = h2o_table['T'][idx-1]
        h_hat_lower = h2o_table['h_hat'][idx-1]
        h_hat = h_hat_close + (h_hat_lower-h_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return h_hat

def oh_table_f(Temp):
    '''
    :param Temp: temperature in K
    :return h_hat: enthalpy in J/mol also known as heat of formation
    '''
    oh_table = pd.read_csv(TABLES / 'OH_table.csv')

    idx = (oh_table['T']-Temp).abs().idxmin()
    T_cloest = oh_table['T'][idx]
    h_hat_close = oh_table['h_hat'][idx]
    if Temp > T_cloest:
        T_upper = oh_table['T'][idx+1]
        h_hat_upper = oh_table['h_hat'][idx+1]
        h_hat = h_hat_close + (h_hat_upper-h_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = oh_table['T'][idx-1]
        h_hat_lower = oh_table['h_hat'][idx-1]
        h_hat = h_hat_close + (h_hat_lower-h_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return h_hat

def h_table_f(Temp):
    '''
    :param Temp: temperature in K
    :return h_hat: enthalpy in J/mol also known as heat of formation
    '''
    h_table = pd.read_csv(TABLES / 'H_table.csv')

    idx = (h_table['T']-Temp).abs().idxmin()
    T_cloest = h_table['T'][idx]
    h_hat_close = h_table['h_hat'][idx]
    if Temp > T_cloest:
        T_upper = h_table['T'][idx+1]
        h_hat_upper = h_table['h_hat'][idx+1]
        h_hat = h_hat_close + (h_hat_upper-h_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = h_table['T'][idx-1]
        h_hat_lower = h_table['h_hat'][idx-1]
        h_hat = h_hat_close + (h_hat_lower-h_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return h_hat

def h2o_tableT(h_hat):
    '''
    :param h_hat: enthalpy in J/mol also known as heat of formation
    :return Temp: temperature in K
    '''
    h2o_table = pd.read_csv(TABLES / 'H2O_table.csv')

    idx = (h2o_table['h_f']-h_hat).abs().idxmin()
    h_hat_cloest = h2o_table['h_f'][idx]
    T_close = h2o_table['T'][idx]
    if h_hat > h_hat_cloest:
        h_hat_upper = h2o_table['h_f'][idx+1]
        T_upper = h2o_table['T'][idx+1]
        T = T_close + (T_upper-T_close)/(h_hat_upper-h_hat_cloest)*(h_hat-h_hat_cloest)
    else:
        h_hat_lower = h2o_table['h_f'][idx-1]
        T_lower = h2o_table['T'][idx-1]
        T = T_close + (T_lower-T_close)/(h_hat_lower-h_hat_cloest)*(h_hat-h_hat_cloest)

    return T

def h2_tableT(h_hat):
    '''
    :param h_hat: enthalpy in J/mol also known as heat of formation
    :return Temp: temperature in K
    '''
    h2_table = pd.read_csv(TABLES / 'Hydrogen_H2_table.csv')
    h_hat = h_hat/4184
    idx = (h2_table['H*-H*_298']-h_hat).abs().idxmin()
    h_hat_cloest = h2_table['H*-H*_298'][idx]
    T_close = h2_table['T'][idx]
    if h_hat > h_hat_cloest:
        h_hat_upper = h2_table['H*-H*_298'][idx+1]
        T_upper = h2_table['T'][idx+1]
        T = T_close + (T_upper-T_close)/(h_hat_upper-h_hat_cloest)*(h_hat-h_hat_cloest)
    else:
        h_hat_lower = h2_table['H*-H*_298'][idx-1]
        T_lower = h2_table['T'][idx-1]
        T = T_close + (T_lower-T_close)/(h_hat_lower-h_hat_cloest)*(h_hat-h_hat_cloest)

    return T

def hf_tableT(h_hat):
    '''
    :param h_hat: enthalpy in J/mol also known as heat of formation
    :return Temp: temperature in K
    '''
    f2_table = pd.read_csv(TABLES / 'HydrogenFluoride_HF_table.csv')
    h_hat = h_hat/4184
    idx = (f2_table['H*-H*_298']-h_hat).abs().idxmin()
    h_hat_cloest = f2_table['H*-H*_298'][idx]
    T_close = f2_table['T'][idx]
    if h_hat > h_hat_cloest:
        h_hat_upper = f2_table['H*-H*_298'][idx+1]
        T_upper = f2_table['T'][idx+1]
        T = T_close + (T_upper-T_close)/(h_hat_upper-h_hat_cloest)*(h_hat-h_hat_cloest)
    else:
        h_hat_lower = f2_table['H*-H*_298'][idx-1]
        T_lower = f2_table['T'][idx-1]
        T = T_close + (T_lower-T_close)/(h_hat_lower-h_hat_cloest)*(h_hat-h_hat_cloest)

    return T


def oh_table_g(Temp):
    '''
    :param Temp: temperature in K
    :return g_hat: Gibbs free energy in J/mol
    '''
    oh_table = pd.read_csv(TABLES / 'OH_table.csv')

    idx = (oh_table['T']-Temp).abs().idxmin()
    T_cloest = oh_table['T'][idx]
    g_hat_close = oh_table['g_hat'][idx]
    if Temp > T_cloest:
        T_upper = oh_table['T'][idx+1]
        g_hat_upper = oh_table['g_hat'][idx+1]
        g_hat = g_hat_close + (g_hat_upper-g_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = oh_table['T'][idx-1]
        g_hat_lower = oh_table['g_hat'][idx-1]
        g_hat = g_hat_close + (g_hat_lower-g_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return g_hat

def h_table_g(Temp):
    '''
    :param Temp: temperature in K
    :return g_hat: Gibbs free energy in J/mol
    '''
    h_table = pd.read_csv(TABLES / 'H_table.csv')

    idx = (h_table['T']-Temp).abs().idxmin()
    T_cloest = h_table['T'][idx]
    g_hat_close = h_table['g_hat'][idx]
    if Temp > T_cloest:
        T_upper = h_table['T'][idx+1]
        g_hat_upper = h_table['g_hat'][idx+1]
        g_hat = g_hat_close + (g_hat_upper-g_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = h_table['T'][idx-1]
        g_hat_lower = h_table['g_hat'][idx-1]
        g_hat = g_hat_close + (g_hat_lower-g_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return g_hat

def h2o_table_g(Temp):
    '''
    :param Temp: temperature in K
    :return g_hat: Gibbs free energy in J/mol
    '''
    h2o_table = pd.read_csv(TABLES / 'H2O_table.csv')

    idx = (h2o_table['T']-Temp).abs().idxmin()
    T_cloest = h2o_table['T'][idx]
    g_hat_close = h2o_table['g_hat'][idx]
    if Temp > T_cloest:
        T_upper = h2o_table['T'][idx+1]
        g_hat_upper = h2o_table['g_hat'][idx+1]
        g_hat = g_hat_close + (g_hat_upper-g_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = h2o_table['T'][idx-1]
        g_hat_lower = h2o_table['g_hat'][idx-1]
        g_hat = g_hat_close + (g_hat_lower-g_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return g_hat

def h2_table_g(Temp):
    '''
    :param Temp: temperature in K
    :return g_hat: Gibbs free energy in J/mol
    '''
    h2_table = pd.read_csv(TABLES / 'Hydrogen_H2_table.csv')

    idx = (h2_table['T']-Temp).abs().idxmin()
    T_cloest = h2_table['T'][idx]
    g_hat_close = h2_table['g_hat'][idx]
    if Temp > T_cloest:
        T_upper = h2_table['T'][idx+1]
        g_hat_upper = h2_table['g_hat'][idx+1]
        g_hat = g_hat_close + (g_hat_upper-g_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = h2_table['T'][idx-1]
        g_hat_lower = h2_table['g_hat'][idx-1]
        g_hat = g_hat_close + (g_hat_lower-g_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return g_hat

def h2o_table_cp(Temp):
    '''
    :param Temp: temperature in K
    :return cp_hat: heat capacity in J/mol-K
    '''
    h2o_table = pd.read_csv(TABLES / 'H2O_table.csv')

    idx = (h2o_table['T']-Temp).abs().idxmin()
    T_cloest = h2o_table['T'][idx]
    cp_hat_close = h2o_table['cp_hat'][idx]
    if Temp > T_cloest:
        T_upper = h2o_table['T'][idx+1]
        cp_hat_upper = h2o_table['cp_hat'][idx+1]
        cp_hat = cp_hat_close + (cp_hat_upper-cp_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = h2o_table['T'][idx-1]
        cp_hat_lower = h2o_table['cp_hat'][idx-1]
        cp_hat = cp_hat_close + (cp_hat_lower-cp_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return cp_hat

def h2_table_cp(Temp):
    '''
    :param Temp: temperature in K
    :return cp_hat: heat capacity in J/mol-K
    '''
    h2_table = pd.read_csv(TABLES / 'Hydrogen_H2_table.csv')

    idx = (h2_table['T']-Temp).abs().idxmin()
    T_cloest = h2_table['T'][idx]
    cp_hat_close = h2_table['cp_hat'][idx]
    if Temp > T_cloest:
        T_upper = h2_table['T'][idx+1]
        cp_hat_upper = h2_table['cp_hat'][idx+1]
        cp_hat = cp_hat_close + (cp_hat_upper-cp_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = h2_table['T'][idx-1]
        cp_hat_lower = h2_table['cp_hat'][idx-1]
        cp_hat = cp_hat_close + (cp_hat_lower-cp_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return cp_hat

def oh_table_cp(Temp):
    '''
    :param Temp: temperature in K
    :return cp_hat: heat capacity in J/mol-K
    '''
    oh_table = pd.read_csv(TABLES / 'OH_table.csv')

    idx = (oh_table['T']-Temp).abs().idxmin()
    T_cloest = oh_table['T'][idx]
    cp_hat_close = oh_table['cp_hat'][idx]
    if Temp > T_cloest:
        T_upper = oh_table['T'][idx+1]
        cp_hat_upper = oh_table['cp_hat'][idx+1]
        cp_hat = cp_hat_close + (cp_hat_upper-cp_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = oh_table['T'][idx-1]
        cp_hat_lower = oh_table['cp_hat'][idx-1]
        cp_hat = cp_hat_close + (cp_hat_lower-cp_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return cp_hat

def h_table_cp(Temp):
    '''
    :param Temp: temperature in K
    :return cp_hat: heat capacity in J/mol-K
    '''
    h_table = pd.read_csv(TABLES / 'H_table.csv')

    idx = (h_table['T']-Temp).abs().idxmin()
    T_cloest = h_table['T'][idx]
    cp_hat_close = h_table['cp_hat'][idx]
    if Temp > T_cloest:
        T_upper = h_table['T'][idx+1]
        cp_hat_upper = h_table['cp_hat'][idx+1]
        cp_hat = cp_hat_close + (cp_hat_upper-cp_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = h_table['T'][idx-1]
        cp_hat_lower = h_table['cp_hat'][idx-1]
        cp_hat = cp_hat_close + (cp_hat_lower-cp_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return cp_hat

def MachSolve(gamma, r_A, guess):
    '''
    Solve for Mach number given gamma, r_A, and guess
    :param gamma: the gamma of the gas
    :param r_A: the ratio of the areas
    :param guess: the initial guess of the Mach number
    :return: Mach number
    '''

    # M = var('M')
    front_cont = ((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))
    sec_cont = (gamma-1)/2
    func = lambda M: (((front_cont)*((1+sec_cont*M**2)**((gamma+1)/(2*(gamma-1)))))/M) - r_A
    sol = fsolve(func, guess)

    return sol

def h2_tableS(Temp):
    '''
    :param Temp: temperature in K
    :return S_hat: entropy in J/mol-K
    '''
    h2_table = pd.read_csv(TABLES / 'Hydrogen_H2_table.csv')

    idx = (h2_table['T']-Temp).abs().idxmin()
    T_cloest = h2_table['T'][idx]
    S_hat_close = h2_table['S'][idx]
    if Temp > T_cloest:
        T_upper = h2_table['T'][idx+1]
        S_hat_upper = h2_table['S'][idx+1]
        S_hat = S_hat_close + (S_hat_upper-S_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = h2_table['T'][idx-1]
        S_hat_lower = h2_table['S'][idx-1]
        S_hat = S_hat_close + (S_hat_lower-S_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return S_hat

def hf_tableS(Temp):
    '''
    :param Temp: temperature in K
    :return S_hat: entropy in J/mol-K
    '''
    hf_table = pd.read_csv(TABLES / 'HydrogenFluoride_HF_table.csv')

    idx = (hf_table['T']-Temp).abs().idxmin()
    T_cloest = hf_table['T'][idx]
    S_hat_close = hf_table['S'][idx]
    if Temp > T_cloest:
        T_upper = hf_table['T'][idx+1]
        S_hat_upper = hf_table['S'][idx+1]
        S_hat = S_hat_close + (S_hat_upper-S_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = hf_table['T'][idx-1]
        S_hat_lower = hf_table['S'][idx-1]
        S_hat = S_hat_close + (S_hat_lower-S_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return S_hat

def f2_tableS(Temp):
    '''
    :param Temp: temperature in K
    :return S_hat: entropy in J/mol-K
    '''
    f2_table = pd.read_csv(TABLES / 'Fluorine_F2_table.csv')

    idx = (f2_table['T']-Temp).abs().idxmin()
    T_cloest = f2_table['T'][idx]
    S_hat_close = f2_table['S'][idx]
    if Temp > T_cloest:
        T_upper = f2_table['T'][idx+1]
        S_hat_upper = f2_table['S'][idx+1]
        S_hat = S_hat_close + (S_hat_upper-S_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = f2_table['T'][idx-1]
        S_hat_lower = f2_table['S'][idx-1]
        S_hat = S_hat_close + (S_hat_lower-S_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)

    return S_hat
# Me = MachSolve(1.237115, 77.5426, 5)
# print('Me = ', Me)

def table_interp(filename, column, Temp):
    '''
    Linear interpolation of one column of a thermo table in temperature
    :param filename: table file name inside ./tables
    :param column: column to interpolate
    :param Temp: temperature in K (must lie inside the table)
    :return: interpolated value
    '''
    table = pd.read_csv(TABLES / filename)
    T_min, T_max = table['T'].min(), table['T'].max()
    if not T_min <= Temp <= T_max:
        raise ValueError(f'T = {Temp} K is outside {filename} ({T_min}-{T_max} K)')
    return float(np.interp(Temp, table['T'], table[column]))


# Enthalpy of formation at 298 K, J/mol, from the tables (HF: -64.8 kcal/mol; F: NIST-JANAF)
H_F_298 = {'H2O': -241845.0, 'H2': 0.0, 'H': 217977.0, 'OH': 38985.0, 'O2': 0.0,
           'HF': -64.8*4184, 'F': 79390.0, 'F2': 0.0}


def sensible_enthalpy(species, Temp):
    '''
    Sensible enthalpy h(T) - h(298 K) in J/mol, built from the tables:
      H2, O2 : tabulated H - H_298 (H2 in kcal/mol)
      H2O    : tabulated h_f column
      H      : constant cp (monatomic ideal gas)
      OH     : from its formation enthalpy, since dh_f(OH) = h_OH - (1/2)h_H2 - (1/2)h_O2
      HF, F2 : tabulated H - H_298 (kcal/mol)
      F      : tabulated H - H298 (J/mol, NIST-JANAF)
    '''
    if species == 'H2':
        return table_interp('Hydrogen_H2_table.csv', 'H*-H*_298', Temp)*4184
    if species == 'O2':
        return table_interp('o2_table.csv', 'H*-H*_298', Temp)
    if species == 'H2O':
        return table_interp('H2O_table.csv', 'h_f', Temp)
    if species == 'H':
        return table_interp('H_table.csv', 'cp_hat', Temp)*(Temp - 298)
    if species == 'OH':
        return (table_interp('OH_table.csv', 'h_hat', Temp) - H_F_298['OH']
                + 0.5*sensible_enthalpy('H2', Temp) + 0.5*sensible_enthalpy('O2', Temp))
    if species == 'HF':
        return table_interp('HydrogenFluoride_HF_table.csv', 'H*-H*_298', Temp)*4184
    if species == 'F2':
        return table_interp('Fluorine_F2_table.csv', 'H*-H*_298', Temp)*4184
    if species == 'F':
        return table_interp('F_table.csv', 'H-H298', Temp)
    raise ValueError(f'no enthalpy data for {species}')


def total_enthalpy(species, Temp):
    '''
    Formation + sensible enthalpy in J/mol (elements H2, O2 are zero at 298 K, the tables' reference)
    '''
    return H_F_298[species] + sensible_enthalpy(species, Temp)
