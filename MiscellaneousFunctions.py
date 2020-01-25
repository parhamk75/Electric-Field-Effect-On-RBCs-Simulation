# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 18:54:55 2020

@author: parha
"""
#%% Import
import Parameters as prms
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
            in 'mol/lit'
    """
    
    return mol_m3_to_mol_lit( D_Ca_i_prev + (prms.a_d* prms.Ts* I_Ca_curr) )


#%% Test

if __name__=='__main__':
    print( Ca_Dens_Calc(prms.D_Ca_i_0, 10e-2) )