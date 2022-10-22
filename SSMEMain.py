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
    FEE = FEE(Hydr_inject_temp, p_starting, Ox_inject_temp, phi)

