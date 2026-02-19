# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 00:09:04 2026

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

# def Halo():
    
    
    
    
# def Posallocate(N, R, M):
#     u = np.random.rand(N)
#     radius = R * np.sqrt(u)
#     #Random Direction for vectors normalized
#     phi = np.random.rand(N) * 2 * np.pi
#     x = radius * np.cos(phi)
#     y = radius * np.sin(phi)
#     #z = np.random.normal(scale=0.01, size=N)
#     z= np.zeros(N)
    
#     Points = np.column_stack((x,y,z))
    
#     masses =  np.full(N, M/N)# 1e25/N)
#     return Points, masses

def Posallocate(N, Rd, Rmax, M, z0):
    radius = []

    while len(radius) < N:
        u1 = np.random.rand(N)
        u2 = np.random.rand(N)
        candidate = -Rd * np.log(u1 * u2)
        valid = candidate[candidate < Rmax]
        radius.extend(valid)

    radius = np.array(radius[:N])
    #Random Direction for vectors normalized
    phi = np.random.rand(N) * 2 * np.pi
    x = radius * np.cos(phi)
    y = radius * np.sin(phi)
    #z= np.zeros(N)
    z = np.random.normal(0, z0, N)
    
    
    Points = np.column_stack((x,y,z))
    
    masses =  np.full(N, M/N)# 1e25/N)
    return Points, masses

Kpc = 3.086e+19
Msol = 1.989e+30


Mstel = 5.04e10 * Msol
Rstel = 26.8 * Kpc


Rd    = 3 * 2.15 * Kpc
Rmax  = Rstel #5 * Rd
z0    = 0.3 * Kpc

z,mass = Posallocate(1000, Rd, Rmax, Mstel, z0)
# z,mass = Posallocate(1000, Rstel , Mstel)
print(z)#,len(mass))
print(len(z[:,0]),len(z[:,1]),len(z[:,2]))
a = Xmatrix = z[:,0]
b = Ymatrix = z[:,1]
c = Zmatrix = z[:,2]
print(type(len(Xmatrix)))
print(len(Ymatrix))
print(len(Zmatrix))

# Xmatrixn = Xmatrix - 3000
# Xmatrix = np.concatenate((Xmatrixn, Xmatrix))
# Ymatrix = np.concatenate((Ymatrix, Ymatrix))
# Zmatrix = np.concatenate((Zmatrix, Zmatrix))




fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')
ax.set_box_aspect([1, 1, 1])
L = 1 * Rstel
ax.set_xlim(-L, L)
ax.set_ylim(-L, L)
ax.set_zlim(-L, L)
ax.set_xticks([-L, 0, L])
ax.set_yticks([-L, 0, L])
ax.set_zticks([-L, 0, L])
ax.view_init(elev=15, azim=45)
ax.scatter(a,b,c, s=3)
scat = ax.scatter(a, b, c, s=5, color = 'tab:blue')

def update(frame):
    ax.view_init(elev=15, azim=frame)
    return scat,

# Animation
ani = FuncAnimation(
    fig,
    update,
    frames=360,     # full rotation
    interval=20     # ~60 fps on screen
)

# Save video
ani.save(
    r"C:/Users\adidu\Documents\Work stuff\Year 3\Computing Project\Old Python Files\Initial_Positions_Rotation.mp4",
    fps=60
)

plt.show()





u = 5*0 #np.random.rand(1000)
print(u)




