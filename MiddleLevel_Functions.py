# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 09:18:51 2020

@author: parha
"""

#%% Import

from Essential_Functions import Ca_Pump_Current_Calc    as _ICaPmp
from Essential_Functions import Ca_Channel_Current_Calc as _ICaChan
from Essential_Functions import m_Inf_Calc              as _m_inf
from Essential_Functions import a_m_Calc                as _alpha
from Essential_Functions import b_m_Calc                as _beta
from Essential_Functions import tau_Calc                as _tau
from Essential_Functions import m_Calc                  as _m


import Parameters as prms
import numpy      as np

from Converters import umol_s_to_pA

#%% Rest Voltage and m_inf Finder by Considering just Ca2+

def Rest_V_and_m_inf_Finder_Ca(D_Ca_Rest= 0, error= 0.05, V_r_range=(-200e-3, 200e-3), m_inf_range= (0, 1)):
    """
        Finds Rest Voltage and m_inf for a Cell by Considering just Ca2+ (Pump and Channel Currents)
        
        - D_Ca_Rest:
            [Ca2+_i] which rest voltage and m_inf should hold it stable\n
            in 'mol/lit'
        - error:
            error of the answer.\n
            It's considered  as the percent of error for the pump current calculated
            in '%'
        - m_inf_range:
            range of m_inf which the search should include\n
            in the form:
                (smaller, larger)
        - V_r_range:
            range of Rest Voltage which the search should include\n
            in the form:
                (smaller, larger)
        - returns:
            (v_r, m_inf_r)\n
            in ('V', 'mol/lit')
    """
    Ip_Ca   = umol_s_to_pA( _ICaPmp(D_Ca_Rest) )
    
    v_r, m_inf_r, err = (0, 0, Ip_Ca* 1)
    """
    #TODO: make the resolution interactive with the error
    
    req_err = error* Ip_Ca          # Required Error
    while(err > req_err):
    """
    Resolution = 1000
    v_rs    = np.linspace(V_r_range[0], V_r_range[1], Resolution)
    m_infs  = np.linspace(m_inf_range[0], m_inf_range[1], Resolution)
    
    prgrs_cnt = 0
    for m in m_infs:
        prgrs_cnt += 1
        prgrs = 100* (prgrs_cnt / len(m_infs))
        print('{:.2f}%'.format(prgrs))
        for v in v_rs:
            tmp_var_1 = _ICaChan(v, m, D_Ca_Rest)            
            cur_err = np.abs( Ip_Ca + tmp_var_1 )
            if cur_err < err:
                v_r     = v
                m_inf_r = m
                err     = cur_err
            elif cur_err == err:
                """
                    #TODO: make a list of all possible answers
                """
                pass
    
    print('\n Ip_Ca = {} pA/Cm2\n V_r = {} mV\n m_inf_r = {}\n err = {}\n err = {} pA/Cm2'.format(Ip_Ca, v_r* 1e3, m_inf_r, err/Ip_Ca, err))
    return v_r, m_inf_r

#%% Calculating the m
    
def ML_m(Vm= 0, m_prev= 0):
    """
    """
    
    alpha   = _alpha(Vm)
    beta    = _beta(Vm)
    tau     = _tau(alpha, beta)
    m_inf   = _m_inf(Vm, alpha, beta)
    

#%% Test

if __name__ == '__main__':
    
    
    # Finding V_rest
    import numpy             as np
    import matplotlib.pyplot as plt
    
    D_Ca    = np.linspace(0.001, 0.5, 50)
    Ip_Ca   = []
    
    for D in D_Ca:
        Ip_Ca.append( _ICaPmp(D) )
        
    plt.plot(D_Ca, Ip_Ca, 'r')
    plt.show()
    
    Rest_V_and_m_inf_Finder_Ca(D_Ca_Rest= prms.D_Ca_i_0,
                               error= None,
                               V_r_range= (105e-3, 106e-3),
                               m_inf_range= (0.4, 0.5)
                               )
        
    
