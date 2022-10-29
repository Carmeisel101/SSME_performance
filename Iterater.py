import numpy as np
import pandas as pd
from sympy import Eq, var, solve, diff
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


def h2_tableH(Temp):
    '''
    :param Temp: temperature in K
    :return h_hat: enthalpy in J/mol also known as heat of formation
    '''
    h2_table = pd.read_csv('./tables/Hydrogen_H2_table.csv')

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
    f2_table = pd.read_csv('./tables/Fluorine_F2_table.csv')

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
    o2_table = pd.read_csv('./tables/o2_table.csv')

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
    hf_table = pd.read_csv('./tables/HydrogenFluoride_HF_table.csv')

    idx = (hf_table['T']-Temp).abs().idxmin()
    T_cloest = hf_table['T'][idx]
    h_hat_close = hf_table['H*-H*_298'][idx]
    if Temp > T_cloest:
        T_upper = hf_table['T'][idx+1]
        h_hat_upper = hf_table['H*-H*_298'][idx+1]
        h_hat = h_hat_close + (h_hat_upper-h_hat_close)/(T_upper-T_cloest)*(Temp-T_cloest)
    else:
        T_lower = hf_table['T'][idx-1]
        h_hat_lower = hf_table['H*-H*_298'][idx-1]
        h_hat = h_hat_close + (h_hat_lower-h_hat_close)/(T_lower-T_cloest)*(Temp-T_cloest)
    h_hat = h_hat*4184

    return h_hat

def h2o_table_f(Temp):
    '''
    :param Temp: temperature in K
    :return h_hat: enthalpy in J/mol also known as heat of formation
    '''
    h2o_table = pd.read_csv('./tables/H2O_table.csv')

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
    oh_table = pd.read_csv('./tables/OH_table.csv')

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
    h_table = pd.read_csv('./tables/H_table.csv')

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
    h2o_table = pd.read_csv('./tables/H2O_table.csv')

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
    h2_table = pd.read_csv('./tables/Hydrogen_H2_table.csv')
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

def oh_table_g(Temp):
    '''
    :param Temp: temperature in K
    :return g_hat: Gibbs free energy in J/mol
    '''
    oh_table = pd.read_csv('./tables/OH_table.csv')

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
    h_table = pd.read_csv('./tables/H_table.csv')

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
    h2o_table = pd.read_csv('./tables/H2O_table.csv')

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
    h2_table = pd.read_csv('./tables/Hydrogen_H2_table.csv')

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
    print('front_cont = ', front_cont)
    sec_cont = (gamma-1)/2
    print('sec_cont = ', sec_cont)
    func = lambda M: (((front_cont)*((1+sec_cont*M**2)**((gamma+1)/(2*(gamma-1)))))/M) - r_A
    sol = fsolve(func, guess)

    return sol


