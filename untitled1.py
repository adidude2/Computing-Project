# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 18:25:52 2026

@author: adidu
"""

import numpy as np

a = np.array([[1,0,0],[0,1,0]])
b = np.array([[1,1,0],[0,1,1]])
c = a * b
print(c)
print(np.sum(a * b, axis = 1))
    