#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 17:01:41 2025

@author: Igor Coropceanu
"""

import numpy as np

h = 6.62607015 * 10**-34
c = 299792458 
e = 1.602176634 * 10**-19
K = 683.002
wl = 550 * 10**-9

EQE=0.4

current_efficiency = EQE * (K*h*c)/(np.pi*wl*e) 
