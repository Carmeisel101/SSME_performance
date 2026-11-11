from Iterater import *
from scipy.optimize import fsolve
import pandas as pd
import numpy as np
from scipy.optimize import brentq

# Molar masses in g/mol, used everywhere (SSMEMain, FEE, OCC)
MOLAR_MASS = {'H2': 2.01588, 'O2': 31.9988, 'H2O': 18.01528, 'OH': 17.00734, 'H': 1.00794}
R_HAT = 8.314  # J/mol-K


def FEE(Hydr_inject_temp, p_starting, Ox_inject_temp, phi):
    '''
    This is a function that returns the calculations of the first enthalpy exchanger
    :param Hydr_inject_temp: temperature of hydrogen injection in K
    :param p_starting: starting pressure in Pa
    :param Ox_inject_temp: temperature of oxygen injection in K
    :param phi: equivalence ratio
    :return T the complete-combustion (no dissociation) temperature in K: the first estimate
    :return del_h_total enthalpy available to heat the products, J per mol O2
    :return del_h2_hat the enthalpy of the injected H2, J/mol
    :return b_hat_table the table of b_hat values
    :return h_reactants enthalpy of the injected propellants, J per mol O2
    '''


    # Van der Waals co-volume b_hat (L/mol)
    b_hat_table = pd.DataFrame({'substance': ['H2', 'H2O', 'O2', 'N2'],
                                'b_hat': [0.0266, 0.03049, 0.0318, 0.0391]})

    b_hat_h2 = b_hat_table['b_hat'][0]
    b_hat_o2 = b_hat_table['b_hat'][2]

    # injected propellant enthalpy (J/mol): sensible + real-gas pressure correction b*(p - p_ref)
    del_h2_hat = h2_tableH(Hydr_inject_temp) + b_hat_h2*(p_starting - 101325)/1000
    del_o2_hat = o2_tableH(Ox_inject_temp) + b_hat_o2*(p_starting - 101325)/1000
    mols_h2 = phi*2
    mols_o2 = 1

    # Complete combustion, per mole of O2: O2 + 2phi H2 -> 2 H2O + (2phi - 2) H2
    # Energy available to heat the products = reactant enthalpy + heat released forming 2 H2O
    del_h_total = del_h2_hat*mols_h2 + del_o2_hat*mols_o2 - 2*H_F_298['H2O']

    # First estimate of the chamber temperature: the single T at which the products'
    # sensible enthalpy equals that energy (no dissociation yet, so this is an upper bound)
    def energy_residual(T):
        return 2*sensible_enthalpy('H2O', T) + (mols_h2 - 2)*sensible_enthalpy('H2', T) - del_h_total

    T = brentq(energy_residual, 300, 4900)
    h_reactants = del_h2_hat*mols_h2 + del_o2_hat*mols_o2

    return T, del_h_total, del_h2_hat, b_hat_table, h_reactants


def ChamberTemperature(FEE_T, p_starting, FO_ratio, h_reactants):
    '''
    Adiabatic chamber temperature with dissociation (the second enthalpy exchanger).
    Finds the T at which the equilibrium products from OCC carry the same enthalpy
    as the injected propellants. Dissociation absorbs energy, so the answer lies
    below the complete-combustion estimate FEE_T.
    :param FEE_T: complete-combustion temperature from FEE in K (upper bracket)
    :param p_starting: chamber pressure in Pa
    :param FO_ratio: fuel to oxidizer mass ratio
    :param h_reactants: enthalpy of the injected propellants, J per mol O2
    :return T_c: chamber temperature in K
    '''
    l = MOLAR_MASS['H2']/(FO_ratio*MOLAR_MASS['O2'])  # mol O2 per mol H2, as in OCC
    species = ['H2O', 'H2', 'H', 'OH']                 # order of OCC's mole table

    def energy_residual(T):
        mols = OCC(T, p_starting, FO_ratio)[4]['mols'][0:4]
        # OCC works per mol H2; dividing by l puts the products per mol O2
        h_products = sum(float(n)/l*total_enthalpy(sp, T) for sp, n in zip(species, mols))
        return h_products - h_reactants

    return brentq(energy_residual, 2000, FEE_T, xtol=1e-6)

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

    # Gibbs free energy of reaction: dG = sum_i (nu''_(i) - nu'_(i)) g_f(i), i.e. products minus reactants
    g_r_h2 = 2*g_f_h - g_f_h2             # H2  -> 2H
    g_r_h2o = g_f_h + g_f_oh - g_f_h2o    # H2O -> H + OH

    k_h2 = ((p_starting/101325)**-1)*np.exp(-g_r_h2/(R_HAT*FEE_T))
    k_h2o = ((p_starting/101325)**-1)*np.exp(-g_r_h2o/(R_HAT*FEE_T))


    # k_h2o = 0.007149
    # k_h2 = 0.011874839

    # mol O2 per mol H2 fed: (1/FO_ratio) kg O2 per kg H2, converted to moles
    l = MOLAR_MASS['H2']/(FO_ratio*MOLAR_MASS['O2'])

    # Products per mole of H2 fed: H2 + l O2 -> N_H2O H2O + N_H2 H2 + N_H H + N_OH OH
    #   H atoms:  2N_H2O + 2N_H2 + N_H + N_OH = 2
    #   O atoms:   N_H2O + N_OH               = 2l
    #   H2O -> H + OH:  X_H X_OH / X_H2O = N_H N_OH / (N_H2O N_total) = k_h2o
    #   H2  -> 2H:      X_H^2 / X_H2     = N_H^2   / (N_H2 N_total)   = k_h2
    # Solved in log(N) so every mole number stays positive (no root picking needed).
    def residuals(log_n):
        n_h2o, n_h2, n_h, n_oh = np.exp(log_n)
        n_tot = n_h2o + n_h2 + n_h + n_oh
        return [2*n_h2o + 2*n_h2 + n_h + n_oh - 2,
                n_h2o + n_oh - 2*l,
                np.log(n_h*n_oh/(n_h2o*n_tot)) - np.log(k_h2o),
                np.log(n_h**2/(n_h2*n_tot)) - np.log(k_h2)]

    # start from complete combustion with a little dissociation
    guess = np.log([2*l, 1 - 2*l, 0.01, 0.01])
    log_n, info, ier, msg = fsolve(residuals, guess, full_output=True, xtol=1e-12)
    if ier != 1 or max(abs(r) for r in residuals(log_n)) > 1e-9:
        raise RuntimeError(f'OCC equilibrium solve did not converge: {msg}')

    N_h2o, N_h2, N_h, N_oh = (float(n) for n in np.exp(log_n))
    N_total = N_h + N_oh + N_h2 + N_h2o

    molar_mass_h2o = MOLAR_MASS['H2O']
    molar_mass_h = MOLAR_MASS['H']
    molar_mass_oh = MOLAR_MASS['OH']
    molar_mass_h2 = MOLAR_MASS['H2']

    LawOfMassAction_df = pd.DataFrame({'substance': ['N_H2O', 'N_H2', 'N_H', 'N_OH', 'Total'],
                                       'mols': [N_h2o, N_h2, N_h, N_oh, N_total]})
    # Mole fractions must be normalized by the total moles: once H and OH form,
    # N_total = N_H2O + N_H2 + N_H + N_OH > 1, so X_(i) != N_(i).
    LawOfMassAction_df['X_(i)'] = LawOfMassAction_df['mols'] / N_total
    # X_(i)*m_hat_(i): the species' contribution to the mixture molar mass (g/mol)
    LawOfMassAction_df['m_hat_(i)'] = pd.Series([molar_mass_h2o, molar_mass_h2, molar_mass_h, molar_mass_oh, 0.0])\
                                      * LawOfMassAction_df['X_(i)']

    # mixture molar mass m_hat = sum_i X_(i) m_hat_(i), over the four species
    m_hat_total = sum(LawOfMassAction_df['m_hat_(i)'][0:4])
    LawOfMassAction_df.loc[4, 'm_hat_(i)'] = m_hat_total
    # mass fractions Y_(i) = X_(i) m_hat_(i) / m_hat
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

    R = R_HAT/m_hat_total
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

    b_hat_table = pd.concat([b_hat_table,
                             pd.DataFrame({'substance': ['F2', 'HF'], 'b_hat': [0.02896, 0.0739]})],
                            ignore_index=True)

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

def OxNozzlePerf(ThroatRadius, ExitRadius, gamma, FEE_T, R, p_starting, p_ambient=0.0):
    '''
    This is a function that returns the performance of the oxidizer nozzle
    :param ThroatRadius: throat radius in inches
    :param ExitRadius: exit radius in inches
    :param gamma: ratio of specific heats (frozen, from OCC)
    :param FEE_T: chamber temperature in K
    :param R: gas constant in kJ/kg-K
    :param p_starting: chamber pressure in Pa
    :param p_ambient: ambient pressure in Pa (0 = vacuum, 101325 = sea level)
    :return: m_dot, v_e, p_e, Thrust, c_T, I_sp, T_e, M_e
    '''

    Area_ratio = (ExitRadius/ThroatRadius)**2
    # change np.float64 values to a float
    gamma = float(gamma)
    FEE_T = float(FEE_T)
    R = float(R)
    R = R*1000 # convert to J/kg-K

    M_e = float(MachSolve(gamma, Area_ratio, 5.0)[0])
    T_e = FEE_T*((1+((gamma-1)/2)*M_e**2)**(-1))
    v_e = M_e * np.sqrt(gamma * R * T_e)
    p_e = p_starting * (1 + ((gamma - 1) / 2) * M_e ** 2) ** (-gamma / (gamma - 1))
    ExitRadius = ExitRadius/39.37 # convert to meters
    ThroatRadius = ThroatRadius/39.37 # convert to


    m_dot = (p_starting * np.pi * (ThroatRadius)**2) / (np.sqrt(FEE_T)) * \
        np.sqrt(gamma/R)*((2/(gamma+1))**((gamma+1)/(2*(gamma-1))))

    A_e = np.pi * ExitRadius**2
    A_t = np.pi * ThroatRadius**2

    # Control-volume momentum balance: F = m_dot*v_e + (p_e - p_a)*A_e.
    # Chamber pressure only enters through c_T's normalization (p_c*A_t), not the thrust itself.
    Thrust = m_dot * v_e + (p_e - p_ambient) * A_e
    c_T = Thrust / (p_starting * A_t)
    I_sp = Thrust / (m_dot * 9.81)

    return m_dot, v_e, p_e, Thrust, c_T, I_sp, T_e, M_e

def FCC(FEE_F):
    '''
    This is a function that returns the performance of the flourine combustion chamber
    :return:
    '''


    # The chemical equation is H2 +F2-> 2HF
    molfrac_h2_react = 1
    molfrac_f2_react = 1
    molfrac_hf_react = 0
    molfrac_h2_prod = 0
    molfrac_f2_prod = 0
    molfrac_hf_prod = 2


    molfrac_h2_total = molfrac_h2_react - molfrac_h2_prod
    molfrac_f2_total = molfrac_f2_react - molfrac_f2_prod
    molfrac_hf_total = molfrac_hf_react - molfrac_hf_prod
    sigma_nu = molfrac_hf_total + molfrac_h2_total + molfrac_f2_total
    # print('sigma_nu = ', sigma_nu)

    # Gibbs free energy of formation of H2, F2, and HF
    del_h_h2 = h2_tableH(FEE_F)
    del_h_f2 = f2_tableH(FEE_F)
    del_h_hf = hf_tableH(FEE_F)

    s_h2 = h2_tableS(FEE_F)
    s_f2 = f2_tableS(FEE_F)
    s_hf = hf_tableS(FEE_F)

    # Gibbs
    del_g_h2 = del_h_h2 - FEE_F*s_h2
    del_g_f2 = del_h_f2 - FEE_F*s_f2
    del_g_hf = del_h_hf - FEE_F*s_hf

    del_g_total = -del_g_f2 - del_g_h2 + del_g_hf
    print('del_g_total = ', del_g_total)