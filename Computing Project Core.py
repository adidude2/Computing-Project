# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 15:34:16 2025

@author: adidu
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
import pandas as pd
import datetime

df = pd.read_excel(r"C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Stuff.xlsx")

#print(df.head())
array_all = df.to_numpy()
Particles = df['Particle Number (N)'].to_numpy()
#print(array_all)


dt = datetime.timedelta(hours=1).total_seconds()
T = 365 * 60 * 60 * 24
TotT = T / dt
#print(f'New one {TotT}' )




N = len(Particles)
Vxs = df['VelocityX'].to_numpy()
Vys = df['VelocityY'].to_numpy()
#print(Vys)

#print(array_all[:,1])


#print(array_all[1,0])
#ATTENTION ARRAYS GO ([ROW,COLUMN])

#Define universal gravitation constant
G=6.6726e-11 #N-m2/kg2
e = 1

def ForceCalc(m1,m2,x1,y1,x2,y2):
    #Fx= G * (-m1 * m2) * (x1 - x2)/(abs(x1 - x2)**2 + e**2)**(3/2)
    #Fy= G * (-m1 * m2) * (y1 - y2)/(abs(y1 - y2)**2 + e**2)**(3/2)
    dx = x2 - x1
    dy = y2 - y1
    r3 = (dx**2 + dy**2 + e**2)**(1.5)
    Fx = G * m1 * m2 * dx / r3
    Fy = G * m1 * m2 * dy / r3
    return Fx, Fy

print()

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


#print(ForceIter())
def AccelCalc(arrx,arry,x,t):
    Fx,Fy = ForceIter(arrx,arry,x)

    SumFx = np.sum(Fx, axis=1)
    SumFy = np.sum(Fy, axis=1)
    dvx = np. zeros(N)
    dvy = np. zeros(N)
    for i in range(0,N):
        mi = array_all[i,1]
        dvx[i] = (SumFx[i] * t) / mi
        dvy[i] = (SumFy[i] * t) / mi
    return dvx, dvy


def COMCalc():  #Works out Centre of mass between two particles ONLY WORKS WITH 2 BODIES
    ms = array_all[:,1]
    xs = array_all[:,2]
    ys = array_all[:,3]
    Sums = np. zeros(N)
    for i in range (0,N):
        Sums[i] = ms[i] * xs[i]
    COMx = (np.sum(Sums))/(np.sum(ms))
    for i in range (0,N):
        Sums[i] = ms[i] * ys[i]   
    COMy = (np.sum(Sums))/(np.sum(ms))
    return COMx, COMy
    
a,b = COMCalc()
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
    A,B = COMCalc()
    Fvel1 = np. zeros(N)
    Fvel2 = np. zeros(N)
    for i in range (0,N):
       Fvel2[i] = np.sqrt(abs((xs[i] - A) * SumFx[i]/ ms[i]))
    Fvel2[0] = -Fvel2[0]
    #for i in range (0,N):
       #Fvel1[i] = (xs[i] - A) * 2 * np.pi/(365*60*60*24)
    return Fvel2 #ONLY WORKS IN THE SUN AND EARTH SYSTEM AND MAYBE JUPITER
    


#print(firstvel())

def InitVelCalc(T, t):
    a = firstvel(T)
    initvx = np. zeros(N)
    initvy = np. zeros(N)
    xs = array_all[:,2]
    ys = array_all[:,3]
    A,B = COMCalc()
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
    Xmatrix = np.empty((int(T), N))
    Ymatrix = np.empty((int(T), N))
    VXmatrix = np.zeros((int(T), N))
    VYmatrix = np.zeros((int(T), N))
    Eradius = np.zeros(int(T)-1)
    xs = array_all[:,2]
    a,b = InitVelCalc(T, t)
    A,B = COMCalc()
    for i in range (0,N):
        Xmatrix[0,i] = array_all[i,2]
        Ymatrix[0,i] = array_all[i,3]
        VXmatrix[0,i] = a[i]
        VYmatrix[0,i] = b[i]
    
    for i in range (0,int(T-1)):
        ax, ay = AccelCalc(Xmatrix, Ymatrix, i, t)
        Eradius[i] = (xs[1] - A) - np.sqrt((Xmatrix[i,1]- A)**2 + (Ymatrix[i,1]- B)**2)
        for j in range (0,N):
            VXmatrix[i+1,j] = VXmatrix[i,j] + ax[j]  #FIX FORCEITER BEFORE CONTINUING
            VYmatrix[i+1,j] = VYmatrix[i,j] + ay[j]
            Xmatrix[i+1,j] = Xmatrix[i,j] + (VXmatrix[i+1,j] * t)
            Ymatrix[i+1,j] = Ymatrix[i,j] + (VYmatrix[i+1,j] * t)
    return VXmatrix, VYmatrix, Xmatrix, Ymatrix, Eradius

Vx, Vy, X, Y, R = BigFunc(TotT, dt)



def dtvary():
    times = np.array([20, 30, 40, 50, 60])
    Y = len(times)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for i in range(0,Y):
        dt = datetime.timedelta(minutes=int(times[i])).total_seconds()
        T = 365 * 60 * 60 * 24
        TotT = T / dt
        Vx, Vy, X, Y, R = BigFunc(TotT, dt)
        x = np.linspace(0, T, len(R))
        p = times[i]
        axes[1].plot(x , R, label = f'dt = {p} mins')
    
    
    axes[0].plot(X[:,0], Y[:,0], label = 'Sun')
    axes[0].plot(X[:,1], Y[:,1], label = 'Earth')
    axes[0].scatter(*COMCalc(), color='red', marker='x', label='COM')
    axes[0].set_xlabel('x position (m)')
    axes[0].set_ylabel('y position (m)')
    axes[0].legend(loc='upper center', bbox_to_anchor=(0.8, 1.1))
    leg = axes[0].legend(loc='upper center', bbox_to_anchor=(0.8, 1.1))
    frame = leg.get_frame()
    frame.set_facecolor("none")   # transparent background
    frame.set_edgecolor("black")  # visible border
    frame.set_linewidth(1.5)      # normal border thickness
    frame.set_alpha(0.2) 
    axes[0].axis('equal')
    #axes[0].set_title(f'Earth-Sun Orbit')
    
    axes[1].set_ylabel('Deviation of Earth Orbit Radius from initial position (m)')
    axes[1].set_xlabel('Time (s)')
    leg = axes[1].legend(loc='upper center', bbox_to_anchor=(0.7, 1.1))
    frame = leg.get_frame()
    frame.set_facecolor("none")   # transparent background
    frame.set_edgecolor("black")  # visible border
    frame.set_linewidth(1.5)      # normal border thickness
    frame.set_alpha(0.2)
    plt.show()
    plt.savefig(r'C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Old Python Files\Earth Orbit Deviation.png', transparent=True)


dtvary()






#print(Vx, Vy, X, Y)
#print(X[0:5,0], Y[0:5,0])
#print(X[0:5,1], Y[0:5,1])
#plt.plot(X[0:5,0], Y[0:5,0])
#plt.plot(X[0:5,1], Y[0:5,1])
#plt.plot(X[:,2], Y[:,2])
#fig, axes = plt.subplots(1, 2, figsize=(12, 6))
#axes[0].plot(X[:,0], Y[:,0], label = 'Sun')
#axes[0].plot(X[:,1], Y[:,1], label = 'Earth')
#axes[0].scatter(*COMCalc(), color='red', marker='x', label='COM')
#axes[0].set_xlabel('x (m)')
#axes[0].set_ylabel('y (m)')
#axes[0].legend()
#axes[0].axis('equal')
#x = np.linspace(1,TotT,len(R))

#axes[1].plot(x , R, label = 'Idiot')
#axes[1].set_ylabel('Radius Of Earth From COM (m)')
#axes[1].set_xlabel('Timestep Number')
#axes[1].legend()
#plt.tight_layout()
#plt.show()
    
    