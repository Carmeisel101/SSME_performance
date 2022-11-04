import pandas as pd
from exchanger import *


if __name__ == '__main__':

    p_starting = 204 #starting pressure in atm
    p_starting = p_starting * 101325 #convert to pascals
    FuelOxidizer_ratio = 0.166
    mol_O2 = 31.998 # g/mol #single molecule mass
    mol_H2 = 2.016 # g/mol #single molecule mass
    phi = (FuelOxidizer_ratio* mol_O2 / mol_H2)/2
    Hydr_inject_temp = 850 #K
    Ox_inject_temp = 530   #K
    ThroatRadius = 5.15 #inches
    ExitRadius = 45.35 #inches
    FEE_T, FEE_del_h_total, del_h2_hat, b_hat_table = FEE(Hydr_inject_temp, p_starting, Ox_inject_temp, phi)
    print('phi = ', phi)

    print('Temperature of the first enthalpy exchanger = ', FEE_T, 'K')
    print('Total enthalpy change in the first enthalpy exchanger = ', FEE_del_h_total, 'J')


    gamma, R, cp_total, cv, LawOfMassAction_df = OCC(FEE_T, p_starting, FuelOxidizer_ratio)
    print('gamma = ', gamma)

    m_dotOx, v_eOx, p_eOx, ThrustOx, c_TOx, I_spOx, T_e_Ox = OxNozzlePerf(ThroatRadius, ExitRadius, gamma, FEE_T, R, p_starting)

    print('Ox mass flow rate = ', m_dotOx, 'kg/s')
    print('Ox exit velocity = ', v_eOx, 'm/s')
    print('Ox exit pressure = ', p_eOx, 'Pa')
    print('Ox thrust = ', ThrustOx, 'N')
    print('Ox thrust coefficient = ', c_TOx)
    print('Ox specific impulse = ', I_spOx, 's')
    print('Ox exit temperature = ', T_e_Ox, 'K')


    print('Flourine Problem')
    b_hat_table, FEE_T_flourine = FEE_F(phi, p_starting, del_h2_hat, b_hat_table, Ox_inject_temp)

    print('Temperature of the first enthalpy exchanger (Flourine) = ', FEE_T_flourine, 'K')

    FCC(FEE_T_flourine)