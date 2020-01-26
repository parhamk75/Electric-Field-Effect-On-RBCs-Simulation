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

def Ca_Dens_Calc_R(D_Ca_i_prev= 0, I_Ca_curr= 0):
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

def Trans_Mem_Current_Calc_R(V_m= 0, m= 0, D_Ca_i= 0):
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
    
#%% m Calculating Function and its essential functions

def a_m_Calc_R(V_m= 0):
    """
        Calculates the alpha
        
        - V_m:
            Vm[n]
            in 'V'
        - returns:
            alpha[n]
            in ????
            #TODO: What is the unit of alpha?
    """
    return 8.5 / (1 + np.exp( ((1e3* V_m)-8) / -12.5 ) ) 
       
def b_m_Calc_R(V_m= 0):
    """
        Calculates the beta
        
        - V_m:
            Vm[n]
            in 'V'
        - returns:
            beta[n]
            in ????
            #TODO: What is the unit of beta?
    """
    return 35 / (1 + np.exp( ((1e3* V_m)+74) / 14.5 ) )
    
def tau_Calc_R(a_m= 0, b_m= 0):
    """
        Calculates the beta
        
        - a_m:
            alpha[n]
            in ????
        - b_m:
            beta[n]
            in ????
        - returns:
            tau[n]
            in ????
            #TODO: What is the unit of tau?
    """
    return 1 / (a_m + b_m)

c_m = 1e3 / ( (4* 0.0777* 92.16* 1e8) / (2563.2) )

def m_Inf_Calc_R(V_m= 0, D_Ca_i= 0, a_m= 0, b_m= 0):
    """
    """
    tmp_var_1 = np.exp(ZF_RT* V_m)
    tmp_var_2 = (1 - tmp_var_1) / (prms.D_Ca_o - (D_Ca_i* tmp_var_1) )
    tmp_var_3 = ( V_m - 0.135 ) / V_m
    tmp_var_4 = a_m / ( a_m + b_m )
    return c_m* tmp_var_1* tmp_var_2* tmp_var_3* tmp_var_4

def m_Calc_R(tau_curr= 0, m_prev= 0, m_inf_curr= 0):
    """
    """
    
    return ((1 - (prms.Ts/tau_curr))* m_prev) + ((prms.Ts/tau_curr)* m_inf_curr)
    
#%% Test

if __name__=='__main__':
    # --> Ca_Dens_Calc_R()
    print( Ca_Dens_Calc_R(prms.D_Ca_i_0, 10e-2) )
    
    # --> Trans_Mem_Current_Calc_R()
    import matplotlib.pyplot as plt
    Im = []
    m = np.linspace(0, 1, 100)
    Vm = np.linspace(-60e-3, 250e-3, 6)
    
    for j in Vm:
        for i in m:
            Im.append( Trans_Mem_Current_Calc_R(j, i, prms.D_Ca_i_0*10) )
            
        plt.plot(m,Im,'.-')
        Im = []
    
    plt.legend(Vm, loc=[0,0])
    plt.hlines(0,-0.1,1.1, colors='fuchsia')
    plt.show()
    