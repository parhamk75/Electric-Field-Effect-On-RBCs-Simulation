# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 18:54:55 2020

@author: parha
"""
#%% Import
import Parameters as prms
import numpy      as np

from Converters import mol_m3_to_mol_lit
from Converters import pA_cm2_to_A_m2
from Converters import umol_s_to_pA
from Converters import uV_to_V

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
            in 'pA/Cm^2'
        - returns:
            [Ca2+_i][n]
            in 'mol/m3'
    """
    I_Ca    = pA_cm2_to_A_m2( I_Ca_curr )
    a_d     = prms.A/(prms.F* prms.z_Ca* prms.Vol)
    
    return D_Ca_i_prev + (a_d* prms.Ts* I_Ca)

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
        Calculates the tau
        
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

def m_Inf_Calc(V_m= 0, a_m= 0, b_m= 0):
    """
    """
    c_m     = 1e3 / ( (4* 0.0777* 92.16* 1e8) / (2563.2) )
    ZF_RT   = 2* 37.45318352                        # (Zp*F)/(R*T)
    
    Semi_D_Ca_i_0   = prms.D_Ca_i_0
    Semi_D_Ca_o     = prms.D_Ca_o
    
    tmp_var_1 = np.exp(ZF_RT* V_m)
    tmp_var_2 = (1 - tmp_var_1) / (Semi_D_Ca_o - (Semi_D_Ca_i_0* tmp_var_1) )
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
            #FIXME: What is the unit of Density in this relation
        - returns:
            Ip_Ca[n]
            in 'μmol/(s.cm2)'
    """
    Fmax    = 5e-6          # in 'μmol/(s.cm2)'
    H       = 0.75
    Km      = 260
    alpha   = Km / D_Ca_i
    
    return Fmax / ( 1 + (alpha**H) )

#%% Ca2+ Total Transmembrane Current Calculating Function
    
def Ca_Total_Current_Calc(Ip_Ca= 0, Ich_Ca= 0):
    """
        Calculates the Total Ca2+ Transmembrane Current Passing through the Ca2+ Channel and Pump
        
        - Ip_Ca:
            Ip_Ca[n]\n
            in 'μmol/(s.cm2)'
        - Ich_Ca:
            Ich_Ca[n]\n
            in 'pA/cm^2'
        - returns:
            I_Ca[n]\n
            in 'pA/Cm^2'
    """
    
    Ip_Ca_converted = umol_s_to_pA(Ip_Ca)
    
    return Ich_Ca + Ip_Ca_converted

#%% Total Transmembrane Current Calculating Function

def Total_Transmem_Current_Calc(*I):
    """
    """
    return np.sum(I)

#%% Stimulation Transmembrane Voltage Calculating Function
    
def Stim_Transmem_Voltage_Calc(n= 0, A_S= 1e-3, f_S= 50, T_ES= 1e-3, T_BS= 0):
    """
        Calculates the n'th Value of Stimulation
        
        - n:
            n
        - A_S:
            A_S
            in 'V'
        - f_S:
            f_S
            in 'Hz'
        - T_ES:
            T_ES
            in 's'
        - T_BS:
            T_BS
            in 's'
        - returns:
            V_S[n]
            in 'V'
    """
    
    if n* prms.Ts < T_BS:
        return 0
    if n* prms.Ts > T_ES:
        return 0
    
    return A_S* np.sin(2* np.pi* f_S* n* prms.Ts)

#%% Transmembrane Voltage Caused by Accumulation/Exodus of Ca2+ ions (Capacitive Voltage) Calculating Function

def Cap_Transmem_Voltage_Calc(V_C_Prev= 0, I_m_Curr= 0):
    """
        Calculates the Transmembrane Voltage Caused by Accumulation/Exodus of Ca2+ ions (Capacitive Voltage)
        
        - V_C_Prev:
            V_C[n-1]
            in 'V'
        - I_m_Curr:
            Im[n]
            in 'pA/Cm^2'
        - returns:
            V_C[n]
            in 'V'
    """
    V_C_Curr_uV = ( -1* (prms.Ts/prms.Cap)* I_m_Curr) + V_C_Prev
    
    return uV_to_V(V_C_Curr_uV)

#%% Transmembrane Voltage Calculating Function

def Transmem_Voltage_Calc(Vs= 0, Vc= 0, Vr= prms.V_r):
    """
    """
    
    return Vc + Vs + Vr

#%% Test

if __name__=='__main__':
    # --> Ca_Dens_Calc()
    print( Ca_Dens_Calc(prms.D_Ca_i_0, 10e-2) )
    
    # --> Ca_Channel_Current_Calc()
    import matplotlib.pyplot as plt
    Ich_Ca = []
    m = np.linspace(0, 1, 100)
    Vm = np.linspace(105e-3, 106e-3, 6)
    
    for j in Vm:
        for i in m:
            Ich_Ca.append( Ca_Channel_Current_Calc(j, i, prms.D_Ca_i_0*10) )
            
        plt.plot(m,Ich_Ca,'.-')
        Ich_Ca = []
    
    plt.legend(Vm, loc=[0,0])
    plt.hlines(0,-0.1,1.1, colors='fuchsia')
    plt.show()
    
    # --> Total_Transmem_Current_Calc()
    print( Total_Transmem_Current_Calc(3,4,5,6,7) )