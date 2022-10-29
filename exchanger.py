from Iterater import *
import pandas as pd


def FEE(Hydr_inject_temp, p_starting, Ox_inject_temp, phi):
    '''
    This is a function that returns the calculations of the first enthalpy exchanger
    :param Hydr_inject_temp: temperature of hydrogen injection in K
    :param p_starting: starting pressure in Pa
    :param Ox_inject_temp: temperature of oxygen injection in K
    :param phi: equivalence ratio
    :return T the temperature of the first enthalpy exchanger
    :return del_h_total the total enthalpy change in the first enthalpy exchanger
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

    return T, del_h_total

def OCC(FEE_T):
    '''
    This is a function that returns the calculations of the oxygen combustion chamber
    :return T the temperature of the oxygen cooling chamber
    :return del_h_total the total enthalpy change in the oxygen cooling chamber
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
    print('sigma_nu = ', sigma_nu)

    molfrac_h2_total = molfrac_h2_prod - molfrac_h2_react
    molfrac_2h_total = molfrac_2h_prod - molfrac_2h_react
    sigma_nu = molfrac_h2_total + molfrac_2h_total
    print('sigma_nu = ', sigma_nu)

    # Gibbs free energy of formation
    g_f_h2o = h2o_table_g(FEE_T)
    g_f_h = h_table_g(FEE_T)
    g_f_oh = oh_table_g(FEE_T)
    g_f_h2 = h2_table_g(FEE_T)

    print('g_f_h2o = ', g_f_h2o)
    print('g_f_h = ', g_f_h)
    print('g_f_oh = ', g_f_oh)
    print('g_f_h2 = ', g_f_h2)





