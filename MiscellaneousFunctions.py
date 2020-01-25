# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 18:54:55 2020

@author: parha
"""
#%% Import
import Parameters as prms
import numpy      as np
from Converters import mol_m3_to_mol_lit


#%% Calcium Intracellular Density Calculating Function

a_d         = prms.A/(prms.F* prms.z_Ca* prms.Vol)

def Ca_Dens_Calc(D_Ca_i_prev= 0, I_Ca_curr= 0):
    """
        Calculate the current intracellular Density of Ca2+ by the relation below:
            [Ca2+_i][n] = [Ca2+_i][n-1] + a_d* Ts* I_Ca[n]
            
        - D_Ca_i_prev:
            [Ca2+_i][n-1]
            in 'mol/lit'
        - I_Ca_curr:
            I_Ca[n]
            in 'A/m^2'
        - returns:
            [Ca2+_i][n]
            in 'mol/m3'
    """
    
    return D_Ca_i_prev + (a_d* prms.Ts* I_Ca_curr)

#%% Transmembrane Ca2+ Current Calculating Function

ZF_RT   = 2* 37.45318352                        # (Zp*F)/(R*T)
a_i     = (4* 0.0777* 92.16* 1e8) / (2563.2) 

def Trans_Mem_Current_Calc(V_m= 0, m= 0, D_Ca_i= 0):
    """
        Calculates the Transmembrane Ca2+ Current due to GHK Model
        
        - V_m:
            Vm[n]
            in 'V'
        - m:
            m[n]
        - D_Ca_i:
            [Ca2+_i][n]
            in 'mol/lit'
        - returns:
            Im[n]
            in 'pA/cm^2'
    """
    if V_m == 0:
        return 0
    else:
        tmp_var = np.exp(ZF_RT* V_m)
        return (a_i* V_m* m) * ( (prms.D_Ca_o - (D_Ca_i* tmp_var) ) / (1 - tmp_var) )
    

#%% Test

if __name__=='__main__':
    # --> Ca_Dens_Calc()
    print( Ca_Dens_Calc(prms.D_Ca_i_0, 10e-2) )
    
    # --> Trans_Mem_Current_Calc()
    import matplotlib.pyplot as plt
    Im = []
    m = np.linspace(0, 1, 100)
    Vm = np.linspace(-60e-3, 250e-3, 6)
    
    for j in Vm:
        for i in m:
            Im.append( Trans_Mem_Current_Calc(j, i, prms.D_Ca_i_0*10) )
            
        plt.plot(m,Im,'.-')
        Im = []
    
    plt.legend(Vm, loc=[0,0])
    plt.hlines(0,-0.1,1.1, colors='fuchsia')
    plt.show()
    