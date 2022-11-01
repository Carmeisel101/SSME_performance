from Iterater import *
import pandas as pd
import numpy as np
from sympy import Eq, solve, var


def FEE(Hydr_inject_temp, p_starting, Ox_inject_temp, phi):
    '''
    This is a function that returns the calculations of the first enthalpy exchanger
    :param Hydr_inject_temp: temperature of hydrogen injection in K
    :param p_starting: starting pressure in Pa
    :param Ox_inject_temp: temperature of oxygen injection in K
    :param phi: equivalence ratio
    :return T the temperature of the first enthalpy exchanger
    :return del_h_total the total enthalpy change in the first enthalpy exchanger
    :return b_hat_table the table of b_hat values
    :return del_h2_hat the enthalpy change of H2
    '''


    h_f_h2o = abs(h2o_table_f(298)) #J/mol
    h_f_h = h_table_f(298)/1000
    h_oh_f = oh_table_f(298)/1000  # KJ/mol

    b_hat_table = pd.DataFrame(columns=['substance', 'b_hat'])
    b_hat_table = b_hat_table.append({'substance': 'H2', 'b_hat': 0.0266}, ignore_index=True)
    b_hat_table = b_hat_table.append({'substance': 'H2O', 'b_hat': 0.03049}, ignore_index=True)
    b_hat_table = b_hat_table.append({'substance': 'O2', 'b_hat': 0.0318}, ignore_index=True)
    b_hat_table = b_hat_table.append({'substance': 'N2', 'b_hat': 0.0391}, ignore_index=True)

    b_hat_h2 = b_hat_table['b_hat'][0]
    b_hat_o2 = b_hat_table['b_hat'][2]

    del_h2_hat = h2_tableH(Hydr_inject_temp) + b_hat_h2*(p_starting - 101325)/1000
    del_o2_hat = o2_tableH(Ox_inject_temp) + b_hat_o2*(p_starting - 101325)/1000
    mols_h2 = phi*2
    mols_o2 = 1

    del_h_total = del_h2_hat*mols_h2 + del_o2_hat*mols_o2
    # At this point the chemical equation is O2 + mols*H2 -> (phi-2)*H2 +2H2O
    del_h_total = del_h_total + 2*h_f_h2o

    # weighting function
    r = 0.17947269916082
    q = (1- r)/2

    T_h2o = h2o_tableT(q*del_h_total)
    T_h2 = h2_tableT(r*del_h_total)

    T = (T_h2o + T_h2)/2

    return T, del_h_total, del_h2_hat, b_hat_table

def OCC(FEE_T, p_starting, FO_ratio):
    '''
    This is a function that returns the calculations of the oxygen combustion chamber
    :param FEE_T: temperature of the first enthalpy exchanger in K
    :param p_starting: starting pressure in Pa
    :param FO_ratio: fuel to oxidizer ratio
    :return gamma: the ratio of specific heats
    :return R: the gas constant
    :return cp_total: the total specific heat capacity
    :return cv : the specific heat capacity at constant volume
    :return LawOfMassAction_df: the dataframe of the law of mass action
    '''

    # The chemical equation is H2O -> H + OH
    molfrac_h2o_react = 1
    molfrac_h_react = 0
    molfrac_oh_react = 0
    molfrac_h2o_prod = 0
    molfrac_h_prod = 1
    molfrac_oh_prod = 1

    # The chemical equation is H2 -> 2H
    molfrac_h2_react = 1
    molfrac_2h_react = 0
    molfrac_2h_prod = 2
    molfrac_h2_prod = 0

    molfrac_h2o_total = molfrac_h2o_prod - molfrac_h2o_react
    molfrac_h_total = molfrac_h_prod - molfrac_h_react
    molfrac_oh_total = molfrac_oh_prod - molfrac_oh_react
    sigma_nu = molfrac_h2o_total + molfrac_h_total + molfrac_oh_total
    # print('sigma_nu = ', sigma_nu)

    molfrac_h2_total = molfrac_h2_prod - molfrac_h2_react
    molfrac_2h_total = molfrac_2h_prod - molfrac_2h_react
    sigma_nu = molfrac_h2_total + molfrac_2h_total
    # print('sigma_nu = ', sigma_nu)

    # Gibbs free energy of formation
    g_f_h2o = h2o_table_g(FEE_T)
    g_f_h = h_table_g(FEE_T)
    g_f_oh = oh_table_g(FEE_T)
    g_f_h2 = h2_table_g(FEE_T)

    # Gibbs free energy of reaction
    g_r_h2 = 2*(g_f_h - g_f_h2)
    g_r_h2o = g_f_h + g_f_h2o - g_f_oh

    k_h2 = ((p_starting/101325)**-1)*np.exp(-g_r_h2/(8.314*FEE_T))
    k_h2o = ((p_starting/101325)**-1)*np.exp(-g_r_h2o/(8.314*FEE_T))


    # k_h2o = 0.007149
    # k_h2 = 0.011874839

    mass_fuel = 2.02 #H2
    mass_oxider = mass_fuel/FO_ratio
    l = mass_oxider/ (2*15.99)

    N_h = var('N_h')
    N_oh = var('N_oh')
    N_h2 = var('N_h2')
    N_h2o = var('N_h2o')


    # solve the system of equations - Law of mass action
    sol = solve([Eq(N_h2o, 2*l), Eq(N_h2o+N_h2, 1), Eq((N_h*N_oh)/((N_h2+N_h2o)*N_h2o), k_h2o),
                 Eq((N_h**2)/((N_h2+N_h2o)*N_h2), k_h2)], [N_h, N_oh, N_h2, N_h2o])
    # split sol into two separate lists
    sol1 = sol[0]
    sol2 = sol[1]

    # if all the solutions are real and positive then assign mols_sol to the solution
    if all(i > 0 for i in sol1) and all(i > 0 for i in sol2):
        mols_sol = sol1
    else:
        mols_sol = sol2

    # print('mols_sol = ', mols_sol)
    N_h = mols_sol[0]
    N_oh = mols_sol[1]
    N_h2 = mols_sol[2]
    N_h2o = mols_sol[3]
    N_total = N_h + N_oh + N_h2 + N_h2o

    molar_mass_h2o = 18.01528
    molar_mass_h = 1.00794
    molar_mass_oh = 17.00734
    molar_mass_h2 = 2.01588
    molar_mass_total = molar_mass_h2o + molar_mass_h + molar_mass_oh + molar_mass_h2

    LawOfMassAction_df = pd.DataFrame({'substance': ['N_H2O', 'N_H2', 'N_H', 'N_OH', 'Total'],
                                       'mols': [N_h2o, N_h2, N_h, N_oh, N_total]})
    LawOfMassAction_df['X_(i)'] = LawOfMassAction_df['mols']
    LawOfMassAction_df['m_hat_(i)'] = pd.Series([molar_mass_h2o, molar_mass_h2, molar_mass_h, molar_mass_oh, molar_mass_total])\
                                      * LawOfMassAction_df['X_(i)']

    # sum m_hat_(i)[i=0 to i=3]
    m_hat_total = sum(LawOfMassAction_df['m_hat_(i)'])-LawOfMassAction_df['m_hat_(i)'][4]
    LawOfMassAction_df['m_hat_(i)'][4] = m_hat_total
    LawOfMassAction_df['Y_(i)'] = LawOfMassAction_df['m_hat_(i)'] / m_hat_total

    h2_cp = h2_table_cp(FEE_T)
    h2o_cp = h2o_table_cp(FEE_T)
    h_cp = h_table_cp(FEE_T)
    oh_cp = oh_table_cp(FEE_T)

    cp_hats = pd.Series([h2o_cp, h2_cp, h_cp, oh_cp])
    cps = cp_hats / pd.Series([molar_mass_h2o, molar_mass_h2, molar_mass_h, molar_mass_oh])
    cp_df = pd.DataFrame({'cp_hat': cp_hats, 'cp': cps})
    cp_df['Y_(i)*cp_(i)'] = cp_df['cp'] * LawOfMassAction_df['Y_(i)'][0:4]

    cp_total = sum(cp_df['Y_(i)*cp_(i)'])

    LawOfMassAction_df = pd.concat([LawOfMassAction_df, cp_df], axis=1)

    R = 8.314/m_hat_total
    cv = cp_total - R

    gamma = cp_total/cv

    return gamma, R, cp_total, cv, LawOfMassAction_df

def FEE_F(phi, p_starting, del_h2_hat, b_hat_table, Ox_inject_temp):
    '''
    This is a function that returns the calculations of the fuel enthalpy exchanger - flourine
    :param phi the equivalence ratio
    :param p_starting combustion pressure in Pa
    :param del_h2_hat heat of formation of H2 from FEE - Oxygen
    :param b_hat_table Van der Waals table of constants
    :param Ox_inject_temp: Oxidizer injection temperature in K
    :return T the temperature of the FEE in K
    :return b_hat_table the Van der Waals table of constants - edited
    '''

    b_hat_table.append({'substance': 'F2', 'b_hat': 0.02896}, ignore_index = True)
    b_hat_table.append({'substance': 'HF', 'b_hat': 0.0739}, ignore_index = True)

    b_hat_f2 = b_hat_table['b_hat'][3]
    F2_inject_temp = Ox_inject_temp

    del_f2_hat = f2_tableH(F2_inject_temp) + b_hat_f2*(p_starting - 101325)/1000
    del_hf_hat = abs(hf_tableH(298))
    del_h_total = phi*del_h2_hat + del_f2_hat
    del_h_total = del_h_total + del_hf_hat


    r = phi - 1
    q = 1

    r_perc = r/(r+q)
    q_perc = 1-r_perc

    T_h2 = h2_tableT(del_h_total)
    T_hf = hf_tableT(q_perc*del_h_total)

    T = (T_hf + T_h2) / 2


    return b_hat_table, T


