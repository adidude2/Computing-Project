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


"ATTENTION ARRAYS GO ([ROW,COLUMN])"


def Posallocate(N, R, M):
    u = np.random.rand(N)
    radius = R * np.sqrt(u)
    #Random Direction for vectors normalized
    phi = np.random.rand(N) * 2 * np.pi
    x = radius * np.cos(phi)
    y = radius * np.sin(phi)
    #z = np.random.normal(scale=0.01, size=N)
    z= np.zeros(N)
    
    Points = np.column_stack((x,y,z))
    
    masses =  np.full(N, M/N)# 1e25/N)
    return Points, masses


#a,mass = Posallocate(1000, 1000, 100000)
#Xmatrix = a[:,0]
#Ymatrix = a[:,1]
#Zmatrix = a[:,2]

def HaloAccel(arrx, arry, arrz, Mhalo, Rvir, c):
    
    r = np.sqrt(arrx**2 + arry**2 + arrz**2)
    
    r_safe = np.where(r == 0, 1e-10, r)
    
    Rs = Rvir / c
    
    V = 4 * np.pi * Rs**3 * (
        
        np.log(1+c) - c/(1+c)
        
        )
    
    rho0 = Mhalo / V
    
    # rhor = rho0/( 
        
    #     (r_safe/Rs)*(1+(r_safe/Rs))**2
        
    #     )
    
    Menc = 4 * np.pi * rho0 * Rs**3 * (
        
        np.log((Rs+r)/Rs) - r/(Rs+r)
        
        )
    
    amag = -G * Menc / (r_safe**2)

    ax = amag * arrx / r_safe
    ay = amag * arry / r_safe
    az = amag * arrz / r_safe
    return ax, ay, az, Menc

    
#a = HaloAccel(Xmatrix, Ymatrix, Zmatrix, 4, 2)


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


def AccelCalc(arrx,arry, arrz, mass, t, Mhalo, Rvir, c):
    Fx, Fy, Fz = force_vectorised(arrx, arry, arrz, mass)
    m = mass
    ax = Fx / m
    ay = Fy / m
    az = Fz / m
    
    
    hax, hay, haz, M = HaloAccel(
       arrx,
       arry,
       arrz,
       Mhalo,
       Rvir,
       c
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


def HaloVels(a, Menc, mass):
    Xmatrix = a[:,0]
    Ymatrix = a[:,1]
    Zmatrix = a[:,2]
    L =len(Xmatrix)
    r = np.sqrt(Xmatrix**2 + Ymatrix**2)
    sort_idx = np.argsort(r)

    r_sorted = r[sort_idx]
    r_safe = np.where(r_sorted == 0, 1e-10, r_sorted)

    # Stellar enclosed mass
    m_sorted = mass[sort_idx]
    Mstar_enc = np.cumsum(m_sorted)

    # Halo enclosed mass
    Mhalo_sorted = Menc[sort_idx]

    # Total enclosed mass
    Mtot = Mhalo_sorted + Mstar_enc

    CircSpeed = np.sqrt(G * Mtot / r_safe)

    # Tangential direction
    a_sorted = a[sort_idx]
    tangent = np.column_stack((-a_sorted[:,1], a_sorted[:,0]))
    norms = np.linalg.norm(tangent, axis=1)
    Vhat = tangent / norms[:, None]

    vels_sorted = np.column_stack((CircSpeed[:,None] * Vhat, np.zeros(L)))

    # Unsort
    unsort_idx = np.argsort(sort_idx)
    return vels_sorted[unsort_idx]

def NonHaloVels(a, mass):
    Xmatrix = a[:,0]
    Ymatrix = a[:,1]
    Zmatrix = a[:,2]
    L =len(Xmatrix)

    
    
    r = np.sqrt(Xmatrix**2 + Ymatrix**2)
    sort_idx = np.argsort(r)
    a_sorted = a[sort_idx]
    r_sorted = np.sqrt(a_sorted[:,0]**2 + a_sorted[:,1]**2)
    r_safe = np.where(r_sorted == 0, 1e-10, r_sorted)
    m_sorted = mass[sort_idx]
    Menc = np.cumsum(m_sorted)
    CircSpeed = np.sqrt(G * Menc / r_safe)
    
    anew = np.column_stack((-a_sorted[:,1], a_sorted[:,0]))
    norms = np.linalg.norm(anew, axis=1)
    Vhat = anew / norms[:, None]
    velsxy = CircSpeed[:, None] * Vhat
    vels_sorted = np.column_stack((velsxy,np.zeros(L)))
    unsort_idx = np.argsort(sort_idx)
    vels = vels_sorted[unsort_idx]
    return vels


a,b = Posallocate(1000, 1000,1000)    
#vel =RandomVels(a, 5)
#print(vel)
#print(len(a[:,0]))


def HugeFunc(T, t, ass, mass, Mhalo, Rvir, c):
    a = ass
    Xmatrix = a[:,0]
    Ymatrix = a[:,1]
    Zmatrix = a[:,2]
    
    hax, hay, haz, Menc = HaloAccel(
       Xmatrix,
       Ymatrix,
       Zmatrix,
       Mhalo,
       Rvir,
       c
    )
    
    

    
    vel = HaloVels(a, Menc, mass)
    #vel = NonHaloVels(a, mass)
    
    VXmatrix = vel[:,0]
    VYmatrix = vel[:,1]
    VZmatrix = vel[:,2]
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
            
        dvx, dvy, dvz = AccelCalc(Xmatrix, Ymatrix, Zmatrix, mass, t, Mhalo, Rvir, c)
        
        
        VXmatrix += 0.5 * dvx
        VYmatrix += 0.5 * dvy
        VZmatrix += 0.5 * dvz

        Xmatrix += VXmatrix * t
        Ymatrix += VYmatrix * t
        Zmatrix += VZmatrix * t
        
        dvx, dvy, dvz = AccelCalc(Xmatrix, Ymatrix, Zmatrix, mass, t, Mhalo, Rvir, c)
        
        VXmatrix += 0.5 * dvx
        VYmatrix += 0.5 * dvy
        VZmatrix += 0.5 * dvz
        
    #print("v_dot_r:", VXmatrix*Xmatrix + VYmatrix*Ymatrix) 
    return SavedCOM, SavedX, SavedY, SavedZ, SavedVX, SavedVY, SavedVZ, mass, SavedSteps 

def plotfig(R, SavedCOM, SavedX, SavedY, SavedZ, SavedSteps, t):
    
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
        
        L = 2 * R
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
        ax.view_init(elev=90, azim=90)#+frame)  # elev is fixed, azim changes each frame
        time_gyr = (SavedSteps[frame] * t) / (1e9 * 3600*24*365)
        ax.set_title(f"Time (Gyr) = {time_gyr:.3f}")
    
    ani = FuncAnimation(
        fig,
        update,
        frames=len(SavedX),
        interval=20  # ms between frames
        )
    
    ani.save(r"C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Old Python Files\nbody COM new.mp4", fps=60)

    plt.show()
    
    return ani, SavedX 


def EnergyPlot(T, t, a, b, R, Mhalo, Rvir, c):
    #fig, axes = plt.subplots(2, 1, figsize=(9, 7))
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot()
    SavedCOM, SavedX, SavedY, SavedZ, Vx, Vy, Vz, mass, SavedSteps = HugeFunc(T, t, a, b, Mhalo, Rvir, c)
    ani, _ = plotfig(R, SavedCOM, SavedX, SavedY, SavedZ, SavedSteps, t)
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
    Virtotarr = np.empty(len(X)-1)
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
        # midX = 0.5 * (X[i] + X[i+1])
        # midY = 0.5 * (Y[i] + Y[i+1])
        # midZ = 0.5 * (Z[i] + Z[i+1])
        
        positions = np.stack([X, Y, Z], axis=1)
        #positions = np.stack([midX, midY, midZ], axis=1)
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
        
        Vir = 0
        Vir = (2*Ktot + Utot)
        Virtotarr[i] = Vir
    Etot = Utotarr + KEtotarr
    Uavg = np.mean(Utotarr)
    Kavg = np.mean(KEtotarr)
    Tavg = np.mean(Etot)
    seconds_per_year = 1e9 * 3600 * 24 * 365
    x = (np.array(SavedSteps[:-1]) * t) / seconds_per_year
    ax.plot(x, Utotarr, label = 'Potential Energy', color = 'tab:red')
    ax.plot(x, KEtotarr, label = 'Kinetic Energy', color = 'tab:purple')
    ax.plot(x, Etot, label = 'Total Energy', color = 'tab:gray')
    ax.plot(x, Virtotarr, label = 'Virial Energy', color = 'tab:blue')
    ax.set_xlabel("Time (Gyrs)")
    ax.set_ylabel("Energy (J)")
    ax.legend()
    fig.savefig(r'C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Old Python Files\Energy Deviation from meanNbody.png', transparent=True)
    print(len(Utotarr))
    print(len(KEtotarr))
    print((Etot - Etot[0]) / Etot[0])
    return Utotarr, KEtotarr, Etot, ani, fig


Msol = 1.989e+30
Kpc = 3.086e+19
Myr = 3600 * 24 * 365 * 1e6

N = 300
Rstel = 2.6 * Kpc
Rvir = 200 * Kpc
Mstel = 5.04e10 * Msol
Mhalo = 0 #0.97e12 * Msol
e = 0.1 * Kpc # softening
c = 9.4

odt = 0.05 * Myr
T_orbit = 212 * Myr
oT = T_orbit
oTotT = oT / odt
print(odt)

a,b = Posallocate(N, Rstel, Mstel)
Xmatrix = a[:,0]
Ymatrix = a[:,1]
Zmatrix = a[:,2]


#ani, SavedX = plotfig(oTotT,  odt, R, a, b)
U, KE, E, ani, fig = EnergyPlot(oTotT, odt , a, b, Rstel, Mhalo, Rvir, c)
#print(SavedX)
#rhor = HaloAccel(Xmatrix, Ymatrix, Zmatrix, Mhalo, Rvir, c)

# """

# SHOW THE PLOT FROM ABOVE SINCE YOURE NOT EVEN SHOWING THE Z DIRECTION ANYWAY
# ALSO REMEMBER TO CHANGE THE INITIAL VELS CALC IN THE BIG FUNC ANDDDD THE ACCEL CALC WHEN IM INCULDING THE HALO AND NOT

# """




