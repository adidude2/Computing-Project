# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 15:34:16 2025

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

"Define universal constants:"
G=6.6726e-11 #N-m2/kg2





def Posallocate(N, R):
    #Random Direction for vectors normalized
    r = np.random.normal(size=(int(N),3))
    rhat = r / np.linalg.norm(r, axis=1)[:, None]

    #Random Radius with uniform volume density
    u = np.random.rand(N)
    r = R * u**(1/3)

    #Combine Radius and Direction
    Points = rhat * r[:, None] 
    
    masses =  np.full(N, 10000)# 1e25/N)
    return Points, masses


a,mass = Posallocate(1000, 1000)
Xmatrix = a[:,0]
Ymatrix = a[:,1]
Zmatrix = a[:,2]

def HaloAccel(arrx, arry, arrz, rho0, rs):
    
    r = np.sqrt(arrx**2 + arry**2 + arrz**2)
    
    r_safe = np.where(r == 0, 1e-10, r)
    
    x = r_safe/rs
    
    Menc = 4*np.pi*rho0*rs**3*(
        
        np.log(1+x) - x/(1+x)
        
        )
    amag = -G * Menc / (r_safe**2)

    ax = amag * arrx / r_safe
    ay = amag * arry / r_safe
    az = amag * arrz / r_safe
    
    return ax, ay, az, Menc

a = HaloAccel(Xmatrix, Ymatrix, Zmatrix, 4, 2)


def force_vectorised(arrx, arry, arrz, mass):
    
    dx = arrx[:, None] - arrx[None, :]
    dy = arry[:, None] - arry[None, :]
    dz = arrz[:, None] - arrz[None, :]
    m = mass

    r2 = dx**2 + dy**2 + dz**2 + e**2
    r3 = r2**1.5

    Fx = -G * (m[:, None] * m[None, :]) * dx / r3
    Fy = -G * (m[:, None] * m[None, :]) * dy / r3
    Fz = -G * (m[:, None] * m[None, :]) * dz / r3

    # Remove self-interaction
    np.fill_diagonal(Fx, 0.0)
    np.fill_diagonal(Fy, 0.0)
    np.fill_diagonal(Fz, 0.0)

    # Net force on each particle
    #return arrx,arry
    return Fx.sum(axis=1), Fy.sum(axis=1), Fz.sum(axis=1)


#print(force_vectorised(a[:,0],a[:,1],a[:,2], mass))


def AccelCalc(arrx,arry, arrz, mass, t):
    Fx, Fy, Fz = force_vectorised(arrx, arry, arrz, mass)
    m = mass
    ax = Fx / m
    ay = Fy / m
    az = Fz / m
    
    
    hax, hay, haz = HaloAccel(
       arrx,
       arry,
       arrz,
       rho0 = 1e-16,   # tune later
       rs   = 5*R
    )

    ax += hax
    ay += hay
    az += haz
    
    
    dvx = ax * t
    dvy = ay * t
    dvz = az * t
    return dvx, dvy, dvz

#print(AccelCalc(Xmatrix, Ymatrix,odt))

def COMCalc(arrx,arry, arrz, mass):
    ms = mass
    COMx = (np.sum(ms * arrx))/(np.sum(ms))
    COMy = (np.sum(ms * arry))/(np.sum(ms))
    COMz = (np.sum(ms * arrz))/(np.sum(ms))
    return COMx, COMy, COMz
    
    
#x,y = COMCalc(Xmatrix, Ymatrix)
#print(x, y)


def RandomVels(a, Vmax):
    Xmatrix = a[:,0]
    Ymatrix = a[:,1]
    Zmatrix = a[:,2]
    L =len(Xmatrix)
    V = np.random.randn(L, 3)
    #a = np.array([[1,0,0],[0,1,0]])
    #V = np.array([[1,1,0],[0,1,1]])
    aa = np.sum(a * a, axis=1)
    dot = np.sum(V * a, axis=1)
    
    V_perp = V - (dot / aa)[:, None] * a
    norms = np.linalg.norm(V_perp, axis=1)
    Vhat = V_perp / norms[:, None]

    # Assign controlled speed
    speeds = np.random.uniform(0, Vmax, size=L)
    vels = speeds[:, None] * Vhat
    speeds = np.linalg.norm(vels, axis=1)
    print("speed min/median/max:",
      speeds.min(),
      np.median(speeds),
      speeds.max())
    return vels

def HaloVels(a, ):
    Xmatrix = a[:,0]
    Ymatrix = a[:,1]
    Zmatrix = a[:,2]
    r = np.sqrt(arrx**2 + arry**2 + Zmatrix**2)

#a,b,mass = Posallocate(1000, 1000)    
#vel =RandomVels(a, 5)
#print(vel)


def HugeFunc(T, t, ass, mass):
    a = ass
    Xmatrix = a[:,0]
    Ymatrix = a[:,1]
    Zmatrix = a[:,2]
    
    
    # Xmatrixl = Xmatrix - 10000
    # Xmatrixr = Xmatrix + 1000
    # Xmatrix = np.concatenate((Xmatrixl, Xmatrixr))
    # Ymatrix = np.concatenate((Ymatrix, Ymatrix))
    # Zmatrix = np.concatenate((Zmatrix, Zmatrix))  
    # mass =  np.concatenate((mass,mass))
    
    #vel = RandomVels(a, 1e-4)
    VXmatrix = np.zeros(len(Xmatrix))
    VYmatrix = np.zeros(len(Xmatrix))
    VZmatrix = np.zeros(len(Xmatrix))
    #VXmatrix = vel[:,0]
    #VYmatrix = vel[:,1]
    #VZmatrix = vel[:,2]
    Radius = np.zeros(len(Xmatrix))
    SavedX = []
    SavedY = []
    SavedZ = []
    SavedVX = []
    SavedVY = []
    SavedVZ = []
    SavedSteps = []
    SavedCOM = []
    for i in range (0,int(T-1)):
        A,B,C = COMCalc(Xmatrix,Ymatrix, Zmatrix, mass)
        if i % 1 == 0:
            SavedCOM.append((A, B, C))
            SavedX.append(Xmatrix.copy())
            SavedY.append(Ymatrix.copy())
            SavedZ.append(Zmatrix.copy())
            SavedVX.append(VXmatrix.copy())
            SavedVY.append(VYmatrix.copy())
            SavedVZ.append(VZmatrix.copy())
            SavedSteps.append(i)
            
        dvx, dvy, dvz = AccelCalc(Xmatrix, Ymatrix, Zmatrix, mass, t)
        Radius = np.sqrt((Xmatrix - A)**2 + (Ymatrix - B)**2 + (Zmatrix - C)**2)
        
        
        VXmatrix += dvx
        VYmatrix += dvy
        VZmatrix += dvz

        Xmatrix += VXmatrix * t
        Ymatrix += VYmatrix * t
        Zmatrix += VZmatrix * t
        
    return SavedCOM, SavedX, SavedY, SavedZ, SavedVX, SavedVY, SavedVZ, mass, SavedSteps 

def plotfig(R, a, b, SavedCOM, SavedX, SavedY, SavedZ, SavedSteps, t):
    
    #SavedCOM, SavedX, SavedY, SavedZ, SavedVX, SavedVY, SavedVZ, mass, SavedSteps = HugeFunc(T, t, a, b)
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect([1, 1, 1])

    scat = ax.scatter([], [], [])
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.set_zlim(-R, R)
    
    def update(frame):
        ax.cla()  # clear axes
        ax.set_box_aspect([1, 1, 1])
        
        A, B, C = SavedCOM[frame]

        X = SavedX[frame]
        Y = SavedY[frame]
        Z = SavedZ[frame]
        
        L = 3 * R
        ax.set_xlim(-L, L)
        ax.set_ylim(-L, L)
        ax.set_zlim(-L, L)
        ax.set_xticks([-L, 0, L])
        ax.set_yticks([-L, 0, L])
        ax.set_zticks([-L, 0, L])
        #ax.set_xticks([])
        #ax.set_yticks([])
        #ax.set_zticks([])
        #ax.set_axis_off()


        ax.scatter(
            X,
            Y,
            Z,
            s=5
        )
        ax.view_init(elev=15, azim=90)#+frame)  # elev is fixed, azim changes each frame
        
        ax.set_title(f"Time (days) = {SavedSteps[frame]*(t/(3600*24))}")
    
    ani = FuncAnimation(
        fig,
        update,
        frames=len(SavedX),
        interval=20  # ms between frames
        )
    
    ani.save(r"C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Old Python Files\nbody COM new.mp4", fps=60)

    plt.show()
    
    return ani, SavedX 


def EnergyPlot(a, b):
    #fig, axes = plt.subplots(2, 1, figsize=(9, 7))
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot()
    dt = datetime.timedelta(days=1).total_seconds()
    T = 365 * 60 * 60 * 24 * 0.5
    TotT = int(T / dt)
    SavedCOM, SavedX, SavedY, SavedZ, Vx, Vy, Vz, mass, SavedSteps = HugeFunc(TotT, dt, a, b)
    ani, _ = plotfig(R, a, b, SavedCOM, SavedX, SavedY, SavedZ, SavedSteps, dt)
    ms = mass
    X = SavedX
    Y = SavedY
    Z = SavedZ
    Np = len(X[0])
    halfposX = np.empty(Np)
    halfposY = np.empty(Np)
    halfposZ = np.empty(Np)
    #Uarr = np.empty(N*(N-1)/2)
    Utotarr = np.empty(len(X)-1)
    KEtotarr = np.empty(len(X)-1)
    for i in range (0,len(X)-1):
        #x = X[i]
        #y = Y[i]
        #z = Z[i]
        Utot = 0
        #dx = x[:,None] - x[None,:]
        #dy = y[:,None] - y[None,:]
        #dz = z[:,None] - z[None,:]
        #r = np.sqrt(dx**2 + dy**2 + dz**2 + e**2)
        Utot = 0
        midX = 0.5 * (X[i] + X[i+1])
        midY = 0.5 * (Y[i] + Y[i+1])
        midZ = 0.5 * (Z[i] + Z[i+1])
        positions = np.stack([midX, midY, midZ], axis=1)
        d = positions[:,None,:] - positions[None,:,:]
        r = np.sqrt(np.sum(d*d,axis=2)+e**2)
        U = -G * (ms[:,None]*ms[None,:]) / r
        np.fill_diagonal(U,0)
        Utot = np.sum(U)/2
        Utotarr[i] = Utot
        
        
        Ktot = 0
        vx = Vx[i+1]
        vy = Vy[i+1]
        vz = Vz[i+1]
        v2 = vx**2 + vy**2 + vz**2
        Ktot = 0.5 * np.sum(ms * v2)
        KEtotarr[i] = Ktot
    Etot = Utotarr + KEtotarr
    Uavg = np.mean(Utotarr)
    Kavg = np.mean(KEtotarr)
    Tavg = np.mean(Etot)
    x = np.linspace(0, T, len(Utotarr))/(60*60*24*365)
    ax.plot(x, Utotarr, label = 'Potential Energy', color = 'tab:red')
    ax.plot(x, KEtotarr, label = 'Kinetic Energy', color = 'tab:purple')
    ax.plot(x, Etot, label = 'Total Energy', color = 'tab:gray')
    ax.legend()
    fig.savefig(r'C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Old Python Files\Energy Deviation from meanNbody.png', transparent=True)
    return Utotarr, KEtotarr, Etot, ani, fig


N = 1000
R = 1000#100000
e = 0.1*R

odt = datetime.timedelta(days=1).total_seconds()
oT = 365 * 60 * 60 * 24 * 1
oTotT = oT / odt


a,b = Posallocate(N, R)
#ani, SavedX = plotfig(oTotT, odt, R, a, b)
U, KE, E, ani, fig = EnergyPlot(a,b)
#print(SavedX)





