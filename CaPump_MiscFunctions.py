# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 00:55:50 2020

@author: parha
"""
#%% Import

import Parameters as prms
import numpy      as np
import matplotlib.pyplot as plt

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
