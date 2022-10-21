import pandas as pd


if __name__ == '__main__':


    FuelOxidizer_ratio = 0.166
    mol_O2 = 31.998 # g/mol
    mol_H2 = 2.016 # g/mol
    phi = (FuelOxidizer_ratio* mol_O2 / mol_H2)/2
    print('phi = ', phi)

