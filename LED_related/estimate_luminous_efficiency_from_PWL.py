#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 18:36:49 2025

@author: Igor Coropceanu

Requires Numpy 2.0+

"""



import numpy as np
import colour as clr
import matplotlib.pyplot as plt
from colour.plotting import *

def calculate_luminous_efficiency_PWL(PWL):
    sd = clr.sd_zeros()
    sd[PWL]=1
    color_efficiency=clr.luminous_efficiency(sd)
    return color_efficiency



# zzz=calculate_luminous_efficiency_PWL(460)

