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


    print('Flourine Problem')
    b_hat_table, FEE_T_flourine = FEE_F(phi, p_starting, del_h2_hat, b_hat_table, Ox_inject_temp)

    print('Temperature of the first enthalpy exchanger (Flourine) = ', FEE_T_flourine, 'K')

    FCC(FEE_T_flourine)