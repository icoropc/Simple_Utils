#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 25 22:36:34 2019

@author: icoropc
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from numpy import pi
from numpy import log
from numpy import tanh
from numpy import exp 


#sns.set()
sns.set_style("white")
sns.set_style("ticks")
sns.set_palette("rocket",4)

A_Au=1.6
e=1.602*10**-19
eps_0=8.854*10**-12
kb=1.381*10**-23
J_to_eV=6.242*10**18
phi_test=0.0256

A_eff=A_Au
Av=6.0221*10**23
nm=10**-9


def NP_Vol(r):
    return 4/3*pi*r**3

def U_VDW_gen(A,R1,R2,r):
    U=-A/3*(R1*R2/(r**2-(R1+R2)**2)+R1*R2/(r**2-(R1-R2)**2)+1/2*log((r**2-(R1+R2)**2)/(r**2-(R1-R2)**2)))
    return U


def U_VDW_genmod(A,R1,R2,r):
    r=r+R1+R2
    A/=J_to_eV
    U=-A/3*(R1*R2/(r**2-(R1+R2)**2)+R1*R2/(r**2-(R1-R2)**2)+1/2*log((r**2-(R1+R2)**2)/(r**2-(R1-R2)**2)))
    print(U)
    return U

def U_elect_as(R,z1,z2,eps,nb,phi_s,T,r):
#    z=1/2*(z1+z2)
    z=1
    r=r+2*R
    lb=e**2/(4*pi*eps*eps_0*kb*T)
    nb*=Av*1000
    l_deb=(eps_0*eps*kb*T/(e**2*(z1**2+z2**2)*nb))**0.5
#    U=8*R*nm/(z**2*lb)*(tanh(e*z*phi_s/(4*kb*T)))**2*exp(-(r-2*R)*nm/l_deb)
    U = 64 * pi * R * nm * nb * (tanh(e*z*phi_s/(4*kb*T)))**2 * exp(-(r - 2 * R) * nm / l_deb) * l_deb ** 2
    print(l_deb)
    return U

def U_correction1(amp, gamma_a, x0, x):
    U_corr_1 = amp / (pi * gamma_a * (1 + ((x - x0) / gamma_a) ** 2))
    return U_corr_1

def U_correction2(amp, x0, x):
    U_corr = amp*(x - x0)**-8
    return U_corr

def U_total_sr(R,z1,z2,eps,nb,phi_s,T,A,r):
    Uel=U_elect_as(R,z1,z2,eps,nb,phi_s,T,r)
    Udisp=U_VDW_genmod(A,R,R,r)
    U_sr = U_correction1(370, 0.01, 0, r)
    U=Uel *1 + Udisp/(kb*T) *1 + U_sr *(R/2.5)
    return U

#Plotting    
xvals=np.linspace(0.02,10,10000)

plt.figure(figsize=[5,5])

legend=[]

# conc_list=[5]
size_list=[2,2.5,3, 3.5]
conc_list=np.asarray([40,75,100, 150])/1000

cond_list=[[50,172],[100,172],[150,172],[200,172]]
cond_list=[[50,172],[50,160],[50,150],[50,140]]

#for i in cond_list:
#    legend.append("[1:3 salt]: " + str(i[0])+" mM, diel: " + str(i[1]))

for i in size_list:
#    legend.append("[1:3 salt]: " + str(50)+" mM, diel: " + str(172))
    legend.append("NP radius: " + str(i)+" nm")

for i in size_list:
    size=i
    yvals=U_total_sr(i,1,3,172,50/1000,.3,298,A_Au,xvals)
    plt.plot(xvals,yvals)
    plt.legend(legend, fontsize=12, loc="lower right")


plt.axhline(-4, linestyle="--",linewidth=1, alpha=0.5)
plt.axvline(0.3, linestyle="--",linewidth=1, alpha=0.5)


plt.xlim(0,5)
plt.ylim(-10,3)




plt.xlabel('Interparticle separation (nm)',fontsize=12)
plt.ylabel('U/kT',fontsize=12)
#
plt.tight_layout()
#plt.savefig("Kubo.png",dpi=300)