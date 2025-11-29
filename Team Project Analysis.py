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
import os

#df = pd.read_excel(r"C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Stuff.xlsx")
#print(df.head())
#array_all = df.to_numpy()
#arrayvy = df['VelocityY'].to_numpy()
#print(arrayvy)

os.chdir("C:/Users/adidu/Documents/waveforms/Logging/Wind Speed Logs")
directory=os.listdir()
print(directory)

da = pd.read_csv('15 2 1 138 0.csv')
db = pd.read_csv('18 2 1 138 0.csv')
dc = pd.read_csv('20 2 1 138 0.csv')
dda = pd.read_csv('18 2 1 119 0.csv')
ddb = pd.read_csv('20 2 1 119 0.csv')
ddc = pd.read_csv('21 2 1 119 0.csv')
ddda = pd.read_csv('20 2 1 96 0.csv')
dddb = pd.read_csv('24 2 1 96 0.csv')
dddc = pd.read_csv('27 2 1 96 0.csv')
dddda = pd.read_csv('22 2 1 72 0.csv')
ddddb = pd.read_csv('25 2 1 72 0.csv')
ddddc = pd.read_csv('31 2 1 72 0.csv')
ddddda = pd.read_csv('23 2 1 54 0.csv')
dddddb = pd.read_csv('20 2 1 54 0.csv')
dddddc = pd.read_csv('28 2 1 54 0.csv')

dar = da.to_numpy()
dbr = db.to_numpy()
dcr = dc.to_numpy()
ddar = dda.to_numpy()
ddbr = ddb.to_numpy()
ddcr = ddc.to_numpy()
dddar = ddda.to_numpy()
dddbr = dddb.to_numpy()
dddcr = dddc.to_numpy()
ddddar = dddda.to_numpy()
ddddbr = ddddb.to_numpy()
ddddcr = ddddc.to_numpy()
dddddar = ddddda.to_numpy()
dddddbr = dddddb.to_numpy()
dddddcr = dddddc.to_numpy()











#print(dda.head())