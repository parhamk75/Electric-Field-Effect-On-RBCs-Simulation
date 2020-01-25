# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 20:01:37 2020

@author: parha
"""

def m3_to_lit(m3= 0):
    """
        Converts 'm^3' to 'lit'
        -m3
    """
    return m3* 1e3

def mol_m3_to_mol_lit(mol_m3= 0):
    """
        Converts 'mol/m^3' to 'mol/lit'
        - mol/m^3
    """
    return mol_m3* 1e-3

def pA_cm2_to_A_m2(pA_cm2= 0):
    """
        Converts 'pA/cm^2' to 'A/m^2'
        - pA/cm^2
    """
    return pA_cm2* 1e-8