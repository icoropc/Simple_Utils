# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 16:44:59 2021

@author: icoropceanu
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def aging_scale(old_lum, new_lum, old_life, slope):
    newlife=old_life*(old_lum/new_lum)**slope
    return newlife

newlife=aging_scale(625, 1000, 22, 1.8)
print(newlife)

# newlife=aging_scale(545, 1000, 0.61, 1.6)
# print(newlife)