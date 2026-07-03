#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 17:01:41 2025

@author: Igor Coropceanu
"""

import numpy as np
import estimate_luminous_efficiency_from_PWL as elc

h = 6.62607015 * 10**-34
c = 299792458 
e = 1.602176634 * 10**-19
K = 683.002
wl = 550 * 10**-9

EQE=0.4

def current_efficiency_from_EQE(EQE, PWL = 550):
    lcor = elc.calculate_luminous_efficiency_PWL(PWL)
    current_efficiency = EQE/100 * lcor* (K*h*c)/(np.pi*wl*e) 
    return current_efficiency

cur_eff=current_efficiency_from_EQE(20, 635)
print(cur_eff)