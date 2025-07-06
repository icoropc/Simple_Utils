#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 17:26:47 2020

@author: icoropc
"""

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from math import sqrt

bulk_dens={}
bulk_dens['Au']=19
bulk_dens['PbS']=7.6
bulk_dens['FeO']=5.24
bulk_dens['ZnSe']=5.27


#Calculating the ligand length for two cases: fully saturated alipathic chains
# and aliphatic chains a chains with one cis double bond 
def calc_ligLength(nC, unsat=0):
    if nC==0:
        return 0
    elif unsat==0:
         ligLength = 0.12*(nC+1)
    elif unsat==1:
            ligLength = sqrt(3)/2*0.12*(nC/2+1)*2
    return ligLength

  

#Calculating the effective hard-sphere particle size given: R - the size of
#the inorganic coore, nC, the number of carbon atoms and the degree of 
#unsaturation (the default is zero, i.e fully saturated chains)
def calc_OPM(R, nC, unsat=0):
    ligLength=calc_ligLength(nC, unsat)
    eff_size=R*(1+3*ligLength/R)**(1/3)
    return eff_size

#Calculating all ratios in terms of B/A particle (A being larger)


def P_Mass(typeP, R,nC, unsat):
    graft_Dens=4.4
    A0=0.2
    dens_Ligand=0.9
    vol_Core=4/3*pi*R**3
    # print(vol_Core)
    mass_Core=vol_Core*bulk_dens[typeP]
    # print(mass_Core)
    ligand_Length=calc_ligLength(nC)
    vol_Ligand = 4*pi*R**2*graft_Dens*A0*ligand_Length
    mass_Ligand=vol_Ligand*dens_Ligand
    # print(mass_Core, mass_Ligand)
    total_Mass=mass_Ligand+mass_Core
    print(mass_Ligand/total_Mass)
    return total_Mass


def rel_Mass(ID1, R1, nC1, dsat1, ID2, R2, nC2, dsat2):
    mass_Ratio_part=P_Mass(ID2, R2, nC2, dsat2)/P_Mass(ID1, R1, nC1, dsat1)
    return mass_Ratio_part
    
conc_A=20
conc_B=20
typeA='FeO'
typeB='Au'
stoic=1

def rel_Mass_total(ID1, R1, nC1, dsat1, R2, nC2, dsat2, ID2, stoic):
    mass_Ratio_total=rel_Mass(ID1, R1, nC1, dsat1, R2, nC2, dsat2, ID2)*stoic
    return mass_Ratio_total

#Final calc of amounts:

def volumes_AB(ID1, R1, nC1, dsat1, R2, nC2, dsat2, ID2, stoic, conc_A, conc_B, mass_A):
        mass_Ratio_total=rel_Mass_total(ID1, R1, nC1, dsat1, R2, nC2, dsat2, ID2, stoic)
        vol_A = mass_A/conc_A
        vol_B = vol_A*mass_Ratio_total*conc_A/conc_B
        return [vol_A, vol_B]


# print(rel_Mass("FeO", 9.5/2, 18, 1, "Au", 4.4/2, 18, 0))
print(rel_Mass("FeO", 19./2., 0, 0, "Au", 14./2., 0, 0))



    