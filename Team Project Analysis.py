# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 16:03:49 2025

@author: adidu
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
import pandas as pd

df = pd.read_excel(r"C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Stuff.xlsx")
print(df.head())
array_all = df.to_numpy()
arrayvy = df['VelocityY'].to_numpy()
print(arrayvy)