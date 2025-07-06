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
bulk_dens['ZnS']=4.09
bulk_dens['ZnSeZnS']=(4.09+5.27)/2

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


def P_Mass_Sphere(typeP, R,nC, unsat):
    graft_Dens=2.7
    # graft_Dens=2.9
    A0=0.22
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


def P_Mass_Cube(typeP, s,nC, unsat):
    graft_Dens=2.5
    # graft_Dens=2.9
    A0=0.22
    dens_Ligand=0.9
    vol_Core=s**3
    # print(vol_Core)
    mass_Core=vol_Core*bulk_dens[typeP]
    # print(mass_Core)
    ligand_Length=calc_ligLength(nC)
    vol_Ligand = 6*s**2*graft_Dens*A0*ligand_Length
    mass_Ligand=vol_Ligand*dens_Ligand
    # print(mass_Core, mass_Ligand)
    total_Mass=mass_Ligand+mass_Core
    print(mass_Ligand/total_Mass)
    return total_Mass

    
P_Mass_Sphere('ZnSeZnS', 3, 18, 1)
P_Mass_Cube('ZnSeZnS', 6.0, 18, 1)
# P_Mass('Au', 2.05, 18, 1)


    