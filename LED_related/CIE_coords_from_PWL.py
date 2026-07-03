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

sd = clr.sd_zeros()
sd[470]=1


# sd2 = clr.SpectralDistribution(nd_array[:,1], nd_array[:,0])


color_eff=clr.luminous_efficiency(sd)

coords=clr.sd_to_XYZ(sd)
xy = clr.XYZ_to_xy(coords)

