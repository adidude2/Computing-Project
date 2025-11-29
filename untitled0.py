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

df = pd.read_csv(r"C:\Users\adidu\Documents\Waveforms\Logging\Latest Fan Setting data\Rotating Fan\645 3 27 Spin.csv")
print(df.head())
array_all = df.to_numpy()
plt.figure(figsize=(20, 5), dpi = 150)
plt.xlabel('Time (ms)')
plt.ylabel('Amplitude (V)')
plt.plot(np.arange(len(array_all[:,1])), array_all[:,1])
#print(np.arange(len(array_all[:,1])))
plt.savefig(r'C:\Users\adidu\Documents\Waveforms\Logging\Latest Fan Setting data\Rotating Fan\Fan Oscillation graph', transparent=False)


