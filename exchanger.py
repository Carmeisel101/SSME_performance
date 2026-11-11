from Iterater import *
from scipy.optimize import fsolve
import pandas as pd
import numpy as np
from scipy.optimize import brentq

# Molar masses in g/mol, used everywhere (SSMEMain, FEE, OCC)
MOLAR_MASS = {'H2': 2.01588, 'O2': 31.9988, 'H2O': 18.01528, 'OH': 17.00734, 'H': 1.00794,
              'F2': 37.99681, 'HF': 20.00634, 'F': 18.99840}
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


def EquilibriumTemperature(T_upper, h_reactants, products):
    '''
    Second enthalpy exchanger: the temperature at which the equilibrium products carry
    the same enthalpy as the injected propellants. Dissociation absorbs energy, so the
    answer lies below the complete-combustion estimate T_upper.
    :param T_upper: complete-combustion temperature in K (upper bracket)
    :param h_reactants: enthalpy of the injected propellants, J per mol oxidizer
    :param products: function T -> {species: mol per mol oxidizer} at equilibrium
    :return: temperature in K
    '''
    def energy_residual(T):
        return sum(n*total_enthalpy(sp, T) for sp, n in products(T).items()) - h_reactants

    return brentq(energy_residual, 2000, T_upper, xtol=1e-6)


def ChamberTemperature(FEE_T, p_starting, FO_ratio, h_reactants):
    '''
    Adiabatic chamber temperature of the H2/O2 chamber, with dissociation
    :param FEE_T: complete-combustion temperature from FEE in K (upper bracket)
    :param p_starting: chamber pressure in Pa
    :param FO_ratio: fuel to oxidizer mass ratio
    :param h_reactants: enthalpy of the injected propellants, J per mol O2
    :return T_c: chamber temperature in K
    '''
    l = MOLAR_MASS['H2']/(FO_ratio*MOLAR_MASS['O2'])  # mol O2 per mol H2, as in OCC
    species = ['H2O', 'H2', 'H', 'OH']                 # order of OCC's mole table

    def products(T):
        mols = OCC(T, p_starting, FO_ratio)[4]['mols'][0:4]
        # OCC works per mol H2; dividing by l puts the products per mol O2
        return {sp: float(n)/l for sp, n in zip(species, mols)}

    return EquilibriumTemperature(FEE_T, h_reactants, products)

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

def FEE_F(Hydr_inject_temp, p_starting, F2_inject_temp, phi):
    '''
    First enthalpy exchanger for the H2/F2 chamber (same equivalence ratio as the oxygen case)
    :param Hydr_inject_temp: temperature of hydrogen injection in K
    :param p_starting: chamber pressure in Pa
    :param F2_inject_temp: temperature of fluorine injection in K
    :param phi: equivalence ratio; stoichiometric is H2 + F2 -> 2HF, so phi = mol H2 per mol F2
    :return T: complete-combustion (no dissociation) temperature in K
    :return del_h_total: enthalpy available to heat the products, J per mol F2
    :return h_reactants: enthalpy of the injected propellants, J per mol F2
    '''
    b_hat_h2 = 0.0266   # Van der Waals co-volume, L/mol
    b_hat_f2 = 0.02896

    del_h2_hat = h2_tableH(Hydr_inject_temp) + b_hat_h2*(p_starting - 101325)/1000
    del_f2_hat = sensible_enthalpy('F2', F2_inject_temp) + b_hat_f2*(p_starting - 101325)/1000
    h_reactants = phi*del_h2_hat + del_f2_hat

    # Complete combustion, per mole of F2: F2 + phi H2 -> 2 HF + (phi - 1) H2
    del_h_total = h_reactants - 2*H_F_298['HF']

    def energy_residual(T):
        return 2*sensible_enthalpy('HF', T) + (phi - 1)*sensible_enthalpy('H2', T) - del_h_total

    T = brentq(energy_residual, 300, 9900)

    return T, del_h_total, h_reactants


def FluorineChamberTemperature(FEE_T, p_starting, phi, h_reactants):
    '''
    Adiabatic chamber temperature of the H2/F2 chamber, with dissociation
    :param FEE_T: complete-combustion temperature from FEE_F in K
    :param p_starting: chamber pressure in Pa
    :param phi: equivalence ratio (mol H2 per mol F2)
    :param h_reactants: enthalpy of the injected propellants, J per mol F2
    :return T_c: chamber temperature in K
    '''
    species = ['HF', 'H2', 'H', 'F']   # order of FCC's mole table

    def products(T):
        mols = FCC(T, p_starting, phi)[4]['mols'][0:4]
        return {sp: float(n) for sp, n in zip(species, mols)}

    # the H and F tables stop at 6000 K
    return EquilibriumTemperature(min(FEE_T, 6000), h_reactants, products)

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

def FCC(T, p_starting, phi):
    '''
    Equilibrium composition and gas properties of the H2/F2 chamber
    Products per mole of F2 fed: F2 + phi H2 -> N_HF HF + N_H2 H2 + N_H H + N_F F
      F atoms:  N_HF + N_F              = 2
      H atoms:  N_HF + 2N_H2 + N_H      = 2 phi
      H2 -> 2H:     N_H^2     / (N_H2 N_total) = K_n
      HF -> H + F:  N_H N_F   / (N_HF N_total) = K_n
    F2 is left out of the products: it is negligible in a fuel-rich chamber at these temperatures.
    :param T: temperature in K
    :param p_starting: chamber pressure in Pa
    :param phi: equivalence ratio (mol H2 per mol F2)
    :return gamma, R (kJ/kg-K), cp_total (kJ/kg-K), cv, LawOfMassAction_df
    '''
    p_ratio = p_starting/101325
    g_f_h = h_table_g(T)
    g_f_h2 = h2_table_g(T)
    g_f_f = table_interp('F_table.csv', 'g_hat', T)
    g_f_hf = table_interp('HydrogenFluoride_HF_table.csv', 'g_hat', T)

    # Gibbs free energy of reaction, products minus reactants; both reactions have sigma = +1
    k_h2 = p_ratio**-1*np.exp(-(2*g_f_h - g_f_h2)/(R_HAT*T))       # H2 -> 2H
    k_hf = p_ratio**-1*np.exp(-(g_f_h + g_f_f - g_f_hf)/(R_HAT*T))  # HF -> H + F

    def residuals(log_n):
        n_hf, n_h2, n_h, n_f = np.exp(log_n)
        n_tot = n_hf + n_h2 + n_h + n_f
        return [n_hf + n_f - 2,
                n_hf + 2*n_h2 + n_h - 2*phi,
                np.log(n_h**2/(n_h2*n_tot)) - np.log(k_h2),
                np.log(n_h*n_f/(n_hf*n_tot)) - np.log(k_hf)]

    # start from complete combustion with a little dissociation
    guess = np.log([2, phi - 1, 0.01, 0.01])
    log_n, info, ier, msg = fsolve(residuals, guess, full_output=True, xtol=1e-12)
    if ier != 1 or max(abs(r) for r in residuals(log_n)) > 1e-9:
        raise RuntimeError(f'FCC equilibrium solve did not converge: {msg}')

    N_hf, N_h2, N_h, N_f = (float(n) for n in np.exp(log_n))
    N_total = N_hf + N_h2 + N_h + N_f

    molar_masses = np.array([MOLAR_MASS['HF'], MOLAR_MASS['H2'], MOLAR_MASS['H'], MOLAR_MASS['F']])
    mols = np.array([N_hf, N_h2, N_h, N_f])
    X = mols/N_total
    m_hat_total = float(np.sum(X*molar_masses))   # mixture molar mass, g/mol
    Y = X*molar_masses/m_hat_total                # mass fractions

    cp_hats = np.array([table_interp('HydrogenFluoride_HF_table.csv', 'c_p', T)*4.184,   # cal -> J
                        h2_table_cp(T),
                        h_table_cp(T),
                        table_interp('F_table.csv', 'cp_hat', T)])
    cps = cp_hats/molar_masses                    # kJ/kg-K
    cp_total = float(np.sum(Y*cps))

    LawOfMassAction_df = pd.DataFrame({'substance': ['N_HF', 'N_H2', 'N_H', 'N_F', 'Total'],
                                       'mols': list(mols) + [N_total],
                                       'X_(i)': list(X) + [1.0],
                                       'Y_(i)': list(Y) + [1.0],
                                       'cp': list(cps) + [cp_total]})

    R = R_HAT/m_hat_total
    cv = cp_total - R
    gamma = cp_total/cv

    return gamma, R, cp_total, cv, LawOfMassAction_df
