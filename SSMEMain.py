import pandas as pd
from exchanger import *


if __name__ == '__main__':

    p_starting = 204 #starting pressure in atm
    p_starting = p_starting * 101325 #convert to pascals
    FuelOxidizer_ratio = 0.166
    # equivalence ratio: (mol H2 per mol O2) / 2, stoichiometric being 2 H2 + O2
    phi = (FuelOxidizer_ratio * MOLAR_MASS['O2'] / MOLAR_MASS['H2'])/2
    Hydr_inject_temp = 850 #K
    Ox_inject_temp = 530   #K
    ThroatRadius = 5.15 #inches
    ExitRadius = 45.35 #inches
    FEE_T, FEE_del_h_total, del_h2_hat, b_hat_table, h_reactants = FEE(Hydr_inject_temp, p_starting, Ox_inject_temp, phi)
    print('phi = ', phi)

    print('Temperature of the first enthalpy exchanger (complete combustion) = ', FEE_T, 'K')
    print('Total enthalpy change in the first enthalpy exchanger = ', FEE_del_h_total, 'J')

    # second enthalpy exchanger: equilibrium (dissociated) products carry the reactant enthalpy
    T_c = ChamberTemperature(FEE_T, p_starting, FuelOxidizer_ratio, h_reactants)
    print('Chamber temperature (equilibrium, with dissociation) = ', T_c, 'K')

    gamma, R, cp_total, cv, LawOfMassAction_df = OCC(T_c, p_starting, FuelOxidizer_ratio)
    print('gamma = ', gamma)

    p_sea_level = 101325 # Pa
    m_dotOx, v_eOx, p_eOx, ThrustOx, c_TOx, I_spOx, T_e_Ox, M_eOx = OxNozzlePerf(ThroatRadius, ExitRadius, gamma, T_c, R, p_starting)
    _, _, _, ThrustOx_SL, c_TOx_SL, I_spOx_SL, _, _ = OxNozzlePerf(ThroatRadius, ExitRadius, gamma, T_c, R, p_starting, p_sea_level)

    print('Ox exit Mach number = ', M_eOx)
    print('Ox mass flow rate = ', m_dotOx, 'kg/s')
    print('Ox exit velocity = ', v_eOx, 'm/s')
    print('Ox exit pressure = ', p_eOx, 'Pa')
    print('Ox exit temperature = ', T_e_Ox, 'K')
    print('Ox thrust (vacuum) = ', ThrustOx, 'N')
    print('Ox thrust coefficient (vacuum) = ', c_TOx)
    print('Ox specific impulse (vacuum) = ', I_spOx, 's')
    print('Ox thrust (sea level) = ', ThrustOx_SL, 'N')
    print('Ox thrust coefficient (sea level) = ', c_TOx_SL)
    print('Ox specific impulse (sea level) = ', I_spOx_SL, 's')


    print('Fluorine Problem')
    # Same injection temperatures and equivalence ratio, with F2 in place of O2:
    # stoichiometric H2 + F2 -> 2HF, so phi is also the mol H2 per mol F2
    F2_inject_temp = Ox_inject_temp
    FEE_T_F, del_h_total_F, h_reactants_F = FEE_F(Hydr_inject_temp, p_starting, F2_inject_temp, phi)
    print('Oxidizer/fuel mass ratio (Fluorine) = ', MOLAR_MASS['F2']/(phi*MOLAR_MASS['H2']))
    print('Temperature of the first enthalpy exchanger (Fluorine, complete combustion) = ', FEE_T_F, 'K')

    T_c_F = FluorineChamberTemperature(FEE_T_F, p_starting, phi, h_reactants_F)
    print('Chamber temperature (Fluorine, equilibrium, with dissociation) = ', T_c_F, 'K')

    gamma_F, R_F, cp_total_F, cv_F, LawOfMassAction_df_F = FCC(T_c_F, p_starting, phi)
    print('gamma (Fluorine) = ', gamma_F)

    m_dotF, v_eF, p_eF, ThrustF, c_TF, I_spF, T_e_F, M_eF = OxNozzlePerf(ThroatRadius, ExitRadius, gamma_F, T_c_F, R_F, p_starting)
    _, _, _, ThrustF_SL, c_TF_SL, I_spF_SL, _, _ = OxNozzlePerf(ThroatRadius, ExitRadius, gamma_F, T_c_F, R_F, p_starting, p_sea_level)

    print('F exit Mach number = ', M_eF)
    print('F mass flow rate = ', m_dotF, 'kg/s')
    print('F exit velocity = ', v_eF, 'm/s')
    print('F exit pressure = ', p_eF, 'Pa')
    print('F exit temperature = ', T_e_F, 'K')
    print('F thrust (vacuum) = ', ThrustF, 'N')
    print('F thrust coefficient (vacuum) = ', c_TF)
    print('F specific impulse (vacuum) = ', I_spF, 's')
    print('F thrust (sea level) = ', ThrustF_SL, 'N')
    print('F thrust coefficient (sea level) = ', c_TF_SL)
    print('F specific impulse (sea level) = ', I_spF_SL, 's')
