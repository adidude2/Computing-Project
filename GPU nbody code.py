# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 15:34:16 2025

@author: adidu
"""
import numpy as np
import cupy as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap
import scipy.integrate as integrate
import matplotlib.animation as animation
import pandas as pd
import datetime
import matplotlib.gridspec as gridspec

df = pd.read_excel(r"C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Stuff.xlsx")

#print(df.head())
array_all = df.to_numpy()
Particles = df['Particle Number (N)'].to_numpy()
#print(array_all)


odt = datetime.timedelta(days=1).total_seconds()
oT = 365 * 60 * 60 * 24 * 12
oTotT = oT / odt
#print(f'New one {TotT}' )

N = len(Particles)
Au = array_all[1,2]

#print(array_all[:,1])


#print(array_all[1,0])
#ATTENTION ARRAYS GO ([ROW,COLUMN])

#Define universal gravitation constant
G=6.6726e-11 #N-m2/kg2
e = 1

#cp.show_config()
#x = cp.arange(10)
#print(x)
print("GPU:", cp.cuda.runtime.getDeviceProperties(0)['computeCapability'])



def ForceCalc(m1,m2,x1,y1,x2,y2):
    dx = x2 - x1
    dy = y2 - y1
    r3 = (dx**2 + dy**2 + e**2)**(1.5)
    Fx = G * m1 * m2 * dx / r3
    Fy = G * m1 * m2 * dy / r3
    return Fx, Fy


def ForceIter(arrx,arry,x):
    array_Fx = np. zeros((N, N))
    array_Fy = np. zeros((N, N))
    for i in range (0,N):
        for j in range(i+1,N):
            mi = array_all[i,1]
            mj = array_all[j,1]
            xi = arrx[x,i]
            xj = arrx[x,j]
            yi = arry[x,i]
            yj = arry[x,j]
            a, b = ForceCalc(mi, mj, xi, yi, xj, yj)
            array_Fx[i,j] = a
            array_Fx[j,i] = -a
            array_Fy[i,j] = b
            array_Fy[j,i] = -b
    return array_Fx, array_Fy    

def force_vectorised(arra_x, arra_y,i):
    """
    x, y : (N,) arrays at a single timestep
    m    : (N,) masses
    """
    arrx = arra_x[i,:]
    arry = arra_y[i,:]
    
    dx = arrx[:, None] - arrx[None, :]
    dy = arry[:, None] - arry[None, :]
    m = array_all[:,1]

    r2 = dx**2 + dy**2 + e**2
    r3 = r2**1.5

    Fx = -G * (m[:, None] * m[None, :]) * dx / r3
    Fy = -G * (m[:, None] * m[None, :]) * dy / r3

    # Remove self-interaction
    np.fill_diagonal(Fx, 0.0)
    np.fill_diagonal(Fy, 0.0)

    # Net force on each particle
    return Fx.sum(axis=1), Fy.sum(axis=1)

#def AccelCalc(arrx,arry,x,t):
def AccelCalc(arra_x,arra_y,i,t):
    #Fx,Fy = ForceIter(arrx,arry,x)
    Fx,Fy = force_vectorised(arra_x, arra_y,i)

    #SumFx = np.sum(Fx, axis=1)
    #SumFy = np.sum(Fy, axis=1)
    m = array_all[:,1]
    ax = Fx / m
    ay = Fy / m
    dvx = ax * t
    dvy = ay * t
    #for i in range(0,N):
     #   mi = array_all[i,1]
      #  dvx[i] = (SumFx[i] * t) / mi
       # dvy[i] = (SumFy[i] * t) / mi
    return dvx, dvy


def COMCalc(arrx,arry,t):
    ms = array_all[:,1]
    #Sums1 = np. zeros(N)
    #Sums2 = np. zeros(N)
    #Sums1[i] = ms * arrx[t]
    COMx = (np.sum(ms * arrx[t]))/(np.sum(ms))
    #for i in range (0,N):
    #    Sums2[i] = ms[i] * arry[t,i]   
    COMy = (np.sum(ms * arry[t]))/(np.sum(ms))
    
    return COMx, COMy
    
    
#a,b = COMCalc()
#print(a, b)
    
def firstvel(T):
    Xmatrix = np.empty((int(T), N))
    Ymatrix = np.empty((int(T), N))
    for i in range (0,N):
        Xmatrix[0,i] = array_all[i,2]
        Ymatrix[0,i] = array_all[i,3]
    Fx, Fy = ForceIter(Xmatrix,Ymatrix,0)
    SumFx = np.sum(Fx, axis=1)
    SumFy = np.sum(Fy, axis=1)
    
    ms = array_all[:,1]
    xs = array_all[:,2]
    A,B = COMCalc(Xmatrix,Ymatrix,0)
    Fvel1 = np. zeros(N)
    Fvel2 = np. zeros(N)
    for i in range (0,N):
       Fvel2[i] = np.sqrt(abs((xs[i] - A) * SumFx[i]/ ms[i]))
    Fvel2[0] = -Fvel2[0]
    #for i in range (0,N):
       #Fvel1[i] = (xs[i] - A) * 2 * np.pi/(365*60*60*24)
    return Fvel2 #ONLY WORKS IN THE SUN AND EARTH SYSTEM AND MAYBE JUPITER
    

def InitVelCalc(T, t):
    Xmatrix = np.empty((int(T), N))
    Ymatrix = np.empty((int(T), N))
    for i in range (0,N):
        Xmatrix[0,i] = array_all[i,2]
        Ymatrix[0,i] = array_all[i,3]
    a = firstvel(T)
    initvx = np. zeros(N)
    initvy = np. zeros(N)
    xs = array_all[:,2]
    ys = array_all[:,3]
    A,B = COMCalc(Xmatrix,Ymatrix,0)
    delthetas = np. zeros(N)
    # This function is made assuming all objects start on the same plane
    for i in range (0,N):
        delthetas[i] = (a[i] * (t/2))/(xs[i] - A) 
    for i in range (0,N):
        initvx[i] = a[i] * np.sin(delthetas[i])
        initvy[i] = a[i] * np.cos(delthetas[i])
    return initvx, initvy
    #do the trig functions to be used in  later functions, gonna make this assuming every body starting with circular motion around a single centre of mass



def BigFunc(T, t):
    Xmatrix = np.zeros((int(T), N))
    Ymatrix = np.zeros((int(T), N))
    VXmatrix = np.zeros((int(T), N))
    VYmatrix = np.zeros((int(T), N))
    Eradius = np.zeros(int(T)-1)
    Jradius = np.zeros(int(T)-1)
    Sradius = np.zeros(int(T)-1)
    a,b = InitVelCalc(T, t)
    xs = array_all[:,2]
    for i in range (0,N):
        Xmatrix[0,i] = array_all[i,2]
        Ymatrix[0,i] = array_all[i,3]
        VXmatrix[0,i] = a[i]
        VYmatrix[0,i] = b[i]
    A1,B1 = COMCalc(Xmatrix,Ymatrix, 0)
    for i in range (0,int(T-1)):
        A,B = COMCalc(Xmatrix,Ymatrix, i)
        dvx, dvy = AccelCalc(Xmatrix, Ymatrix, i, t)
        #Eradius[i] = (xs[1] - A1) - np.sqrt((Xmatrix[i,1]- A)**2 + (Ymatrix[i,1]- B)**2)
        #Jradius[i] = (xs[2] - A1) - np.sqrt((Xmatrix[i,2]- A)**2 + (Ymatrix[i,2]- B)**2)
        #Sradius[i] = (xs[0] - A1) - np.sqrt((Xmatrix[i,0]- A)**2 + (Ymatrix[i,0]- B)**2)
        Eradius[i] = np.sqrt((Xmatrix[i,1]- A)**2 + (Ymatrix[i,1]- B)**2)
        Jradius[i] = np.sqrt((Xmatrix[i,2]- A)**2 + (Ymatrix[i,2]- B)**2)
        Sradius[i] = np.sqrt((Xmatrix[i,0]- A)**2 + (Ymatrix[i,0]- B)**2)
        VXmatrix[i+1] = VXmatrix[i] + dvx  
        VYmatrix[i+1] = VYmatrix[i] + dvy
        Xmatrix[i+1] = Xmatrix[i] + (VXmatrix[i+1] * t)
        Ymatrix[i+1] = Ymatrix[i] + (VYmatrix[i+1] * t)
    Ermean = (np.mean(Eradius)-(xs[1]-A1))/(xs[1]-A1)
    Jrmean = (np.mean(Jradius)/(xs[2]-A1))-1
    Srmean = (-np.mean(Sradius)/(xs[0]-A1))-1
    Eradius = (Eradius / (xs[1]-A1))-1
    Jradius = (Jradius / (xs[2]-A1))-1
    Sradius = -(Sradius / (xs[0]-A1))-1
    test = Sradius[0]
    #test = xs[1]
    return VXmatrix, VYmatrix, Xmatrix, Ymatrix, Eradius, Jradius, Sradius, Ermean, Jrmean, Srmean, test
    #return VXmatrix, VYmatrix, Xmatrix, Ymatrix, Eradius, Sradius, Ermean, Srmean, test

#Vx, Vy, X, Y, R1, R2, R3, Em, Jm, Sm, t = BigFunc(oTotT, odt)
#print(t)
#plt.plot(X[:,1]/Au, Y[:,1]/Au)
#plt.plot(X[:,0]/Au, Y[:,0]/Au)
#plt.plot(X[:,2]/Au, Y[:,2]/Au)