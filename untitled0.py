# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 16:17:37 2025

@author: adidu
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
import pandas as pd
import os


#os.chdir("C:/Users/adidu/Documents/waveforms/Logging/Wind Speed Logs")
#directory=os.listdir()
#print(directory)


df  = pd.read_csv(r"C:\Users\adidu\Documents\Code Stuff\Waveforms\Logging\Latest Fan Setting data\Rotating Fan\645 3 27 Spin.csv")
#df = pd.read_csv(r"C:\Users\adidu\Documents\Waveforms\Logging\Latest Fan Setting data\Rotating Fan\645 3 27 Spin.csv")

#df = pd.read_csv(r"C:\Users\adidu\Documents\Waveforms\Logging\Latest Fan Setting data\Rotating Fan\Trial Spin.csv")
print(df.head())

dt = 1.0 / 37   # in seconds
t = np.arange(len(df)) * dt
print(len(t))
#print(len(array_all[:,1]))
array_all = df.to_numpy()
plt.figure(figsize=(30, 5), dpi = 100)
#plt.spines['top'].set_visible(False)
#plt.spines['right'].set_visible(False)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude (V)')
#np.arange(len(array_all[:,1]))
plt.plot(t, array_all[:,1])
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#print(np.arange(len(array_all[:,1])))
plt.savefig(r"C:\Users\adidu\Documents\Code Stuff\Waveforms\Logging\Latest Fan Setting data\Rotating Fan\Fan Oscillation graph Filled in.png", transparent=False)


