# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 18:25:52 2026

@author: adidu
"""

import numpy as np

data = np.array([10, 20, 30, 40, 50])

for t in range(0, len(data), 2):
    value = data[t]
    print(value)
    