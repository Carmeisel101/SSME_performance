import pandas as pd


def h2_tableH(Temp):
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
    print('h_hat = ', h_hat)

    return h_hat

def f2_tableH(Temp):
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

def hf_tableH(Temp):
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

# T = h2_tableH(850)