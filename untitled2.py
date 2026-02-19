# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 21:11:38 2026

@author: adidu
"""
import numpy as np
#import cupy as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap
import scipy.integrate as integrate
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import pandas as pd
import datetime
import matplotlib.gridspec as gridspec


Msol = 1.989e+30
Kpc = 3.086e+19
Mhalo = 0.97e12 * Msol
Rvir  = 200 * Kpc
c     = 9.4

Rs = Rvir / c

r_min = 0.0001 * Kpc
r_max = 200 * Rs
Npts  = 1000

r = np.logspace(np.log10(r_min),
                np.log10(r_max),
                Npts)

Rs = Rvir / c

V = 4 * np.pi * Rs**3 * (
    
    np.log(1+c) - c/(1+c)
    
    )

rho0 = Mhalo / V

rhor = rho0/( 
    
    (r/Rs)*(1+(r/Rs))**2
    
    )

rho_ast = rhor * (Kpc**3 / Msol)
r_kpc   = r / Kpc

x = r / Rs

plt.figure(figsize=(6,5))
#plt.loglog(r_kpc, rho_ast, color = 'black')
plt.loglog(x, rho_ast, color='black')

#plt.xlabel("r [kpc]")
plt.xlabel(r"$r/R_s$")
plt.ylabel(r"$\rho(r)$ [$M_\odot$ kpc$^{-3}$]")
plt.title("NFW Density Profile")

plt.grid(True, which="both", ls="--", alpha=0.4)
plt.tight_layout()
plt.show()