# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 18:10:16 2020

@author: parha
"""
#%% Import

import numpy as np

#%% Global Computation Parameter

Ts          = 1e-3                      # Infinitesimal time step duration in 's' (Sampling Time)
D_Ca_o      = 0.0024                    # Extracellular Ca2+ density

#%% Dimention

d           = 6 * 1e-6                  # Cell diametere in 'm'
A           = d ** 2                    # One face area in 'm^2'
Vol         = d ** 3                    # Cell volume in 'm^3'

#%% Initial Value

D_Ca_i_0    = 0.0024/np.exp(270/26.7)   # Initial Ca2+ intracellular density in 'mol/lit'


#%% Constant

z_Ca        = 2                         # Ca2+ electric charge
F           = 96487                     # Faraday Constant in 'C/mol'
R           = 8.314                     # Gas Constant in 'J/K' @27°C
T           = 300                       # Temperature in 'K' @27°C

#%% Derived Coefficient

RT_F        = (R*T)/F                   # in 'V' @27°C


#%% Test

if __name__ == '__main__':
    print( R* T )