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
    a_d     = prms.A/(prms.F* prms.z_Ca* prms.Vol)
    
    return D_Ca_i_prev + (a_d* prms.Ts* I_Ca_curr)

#%% Transmembrane Ca2+ Current through Ca2+ Channel Calculating Function

def Ca_Channel_Current_Calc(V_m= 0, m= 0, D_Ca_i= 0):
    """
        Calculates the Transmembrane Ca2+ Current through the Ca2+ Channels due to GHK Model
        
        - V_m:
            Vm[n]
            in 'V'
        - m:
            m[n]
        - D_Ca_i:
            [Ca2+_i][n]
            in 'mol/lit'
        - returns:
            Ich_Ca[n]
            in 'pA/cm^2'
    """
    
    ZF_RT   = 2* 37.45318352                        # (Zp*F)/(R*T)
    a_i     = (4* 0.0777* 92.16* 1e8) / (2563.2) 

    if V_m == 0:
        return 0
    else:
        tmp_var = np.exp(ZF_RT* V_m)
        return (a_i* V_m* m) * ( (prms.D_Ca_o - (D_Ca_i* tmp_var) ) / (1 - tmp_var) )
    
#%% m Calculating Function and its essential functions

def a_m_Calc(V_m= 0):
    """
        Calculates the alpha
        
        - V_m:
            Vm[n]
            in 'V'
        - returns:
            alpha[n]
            in '1/ms'
    """
    return 8.5 / (1 + np.exp( ((1e3* V_m)-8) / -12.5 ) ) 
       
def b_m_Calc(V_m= 0):
    """
        Calculates the beta
        
        - V_m:
            Vm[n]
            in 'V'
        - returns:
            beta[n]
            in '1/ms'
    """
    return 35 / (1 + np.exp( ((1e3* V_m)+74) / 14.5 ) )
    
def tau_Calc(a_m= 0, b_m= 0):
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
            in 'ms'
    """
    return 1 / (a_m + b_m)

def m_Inf_Calc(V_m= 0, D_Ca_i= 0, a_m= 0, b_m= 0):
    """
    """
    c_m = 1e3 / ( (4* 0.0777* 92.16* 1e8) / (2563.2) )
    ZF_RT   = 2* 37.45318352                        # (Zp*F)/(R*T)
    
    tmp_var_1 = np.exp(ZF_RT* V_m)
    tmp_var_2 = (1 - tmp_var_1) / (prms.D_Ca_o - (D_Ca_i* tmp_var_1) )
    tmp_var_3 = ( V_m - 0.135 ) / V_m
    tmp_var_4 = a_m / ( a_m + b_m )
    
    return c_m* tmp_var_1* tmp_var_2* tmp_var_3* tmp_var_4

def m_Calc(tau_curr= 0, m_prev= 0, m_inf_curr= 0):
    """
    """
    
    return ((1 - (prms.Ts/tau_curr))* m_prev) + ((prms.Ts/tau_curr)* m_inf_curr)
    
#%% Calcium Pump Current Calculating Function

def Ca_Pump_Current_Calc(D_Ca_i= 0):
    """
        Calculates the Transmembrane Ca2+ Current through the Ca2+ Pumps
        - D_Ca_i:
            [Ca2+_i][n]:
            in ????
            #TODO: What is the unit of Density in this relation
        - returns:
            Ip_Ca[n]
            in 'μmol/(s.cm2)'
    """
    Fmax    = 5e-6          # in 'μmol/(s.cm2)'
    H       = 0.75
    Km      = 260
    alpha   = Km / D_Ca_i
    
    return Fmax / ( 1 + (alpha**H) )

#%% Test

if __name__=='__main__':
    # --> Ca_Dens_Calc()
    print( Ca_Dens_Calc(prms.D_Ca_i_0, 10e-2) )
    
    # --> Ca_Channel_Current_Calc()
    import matplotlib.pyplot as plt
    Ich_Ca = []
    m = np.linspace(0, 1, 100)
    Vm = np.linspace(-60e-3, 250e-3, 6)
    
    for j in Vm:
        for i in m:
            Ich_Ca.append( Ca_Channel_Current_Calc(j, i, prms.D_Ca_i_0*10) )
            
        plt.plot(m,Ich_Ca,'.-')
        Ich_Ca = []
    
    plt.legend(Vm, loc=[0,0])
    plt.hlines(0,-0.1,1.1, colors='fuchsia')
    plt.show()
    