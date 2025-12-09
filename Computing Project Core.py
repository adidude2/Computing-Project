# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 15:34:16 2025

@author: adidu
"""
import numpy as np
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


odt = datetime.timedelta(hours=1).total_seconds()
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


def COMCalc(arrx,arry,t):
    ms = array_all[:,1]
    Sums1 = np. zeros(N)
    Sums2 = np. zeros(N)
    for i in range (0,N):
        Sums1[i] = ms[i] * arrx[t,i]
    COMx = (np.sum(Sums1))/(np.sum(ms))
    for i in range (0,N):
        Sums2[i] = ms[i] * arry[t,i]   
    COMy = (np.sum(Sums2))/(np.sum(ms))
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
        ax, ay = AccelCalc(Xmatrix, Ymatrix, i, t)
        #Eradius[i] = (xs[1] - A1) - np.sqrt((Xmatrix[i,1]- A)**2 + (Ymatrix[i,1]- B)**2)
        #Jradius[i] = (xs[2] - A1) - np.sqrt((Xmatrix[i,2]- A)**2 + (Ymatrix[i,2]- B)**2)
        #Sradius[i] = (xs[0] - A1) - np.sqrt((Xmatrix[i,0]- A)**2 + (Ymatrix[i,0]- B)**2)
        Eradius[i] = np.sqrt((Xmatrix[i,1]- A)**2 + (Ymatrix[i,1]- B)**2)
        Jradius[i] = np.sqrt((Xmatrix[i,2]- A)**2 + (Ymatrix[i,2]- B)**2)
        Sradius[i] = np.sqrt((Xmatrix[i,0]- A)**2 + (Ymatrix[i,0]- B)**2)
        for j in range (0,N):
            VXmatrix[i+1,j] = VXmatrix[i,j] + ax[j]  
            VYmatrix[i+1,j] = VYmatrix[i,j] + ay[j]
            Xmatrix[i+1,j] = Xmatrix[i,j] + (VXmatrix[i+1,j] * t)
            Ymatrix[i+1,j] = Ymatrix[i,j] + (VYmatrix[i+1,j] * t)
    Ermean = (np.mean(Eradius)/(xs[1]-A1))-1
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


def dtvary():
    times = np.array([12])
    Y = len(times)
    fig, axes = plt.subplots(1, 2, figsize=(21, 7))
    cmap = LinearSegmentedColormap.from_list("orange_black", ["orange", "black"])
    t = np.linspace(0, 0.5, Y)
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for i in range(0,Y):
        dt = datetime.timedelta(hours=int(times[i])).total_seconds()
        color = cmap(t[-1-i])
        T = 365 * 60 * 60 * 24 * 50
        TotT = T / dt
        Vx, Vy, X, Y, R1, R2, R3, Em, Jm, Sm, t = BigFunc(int(TotT), dt)
        #Vx, Vy, X, Y, R1, R3, Em, Sm, t = BigFunc(int(TotT), dt)
        x = np.linspace(0, T, len(R1)) / (60*60*24*365)
        p = times[i]
        axes[1].plot(x , R1,label = 'Earth') #f'dt = {p} mins')
        axes[1].plot(x , R3*4,label = 'Sun' ) #f'dt = {p} mins')
        axes[1].plot(x , R2*4,label = 'Jupiter') #f'dt = {p} mins')
    axes[0].plot(X[:,1]/Au, Y[:,1]/Au)
    axes[0].plot(X[:,0]/Au, Y[:,0]/Au)
    #axes[0].plot(X[:,0]/1000, Y[:,0]/1000, color = "#ff7f0e")
    axes[0].plot(X[:,2]/Au, Y[:,2]/Au)
    axes[0].plot(X[0,1]/Au, Y[0,1]/Au, marker='o', label = 'Earth', color = "#2178b4", markeredgecolor='black')
    axes[0].plot(X[0,0]/Au, Y[0,0]/Au, marker='o', label = 'Sun', color = "#ff7f0e", markeredgecolor='black')
    axes[0].plot(X[0,2]/Au, Y[0,2]/Au, marker='o', label = 'Jupiter', color = "#37a33c", markeredgecolor='black')
    axes[0].set_xlabel('x position (Au)')
    axes[0].set_ylabel('y position (Au)')
    axes[1].set_xlabel('Time (years)')
    axes[1].set_ylabel('Relative distance From COM')
    axes[0].legend(loc='upper center', bbox_to_anchor=(0.8, 1.1))
    leg = axes[0].legend(loc='upper center', bbox_to_anchor=(0.8, 1.1))
    frame = leg.get_frame()
    frame.set_facecolor("none")   # transparent background
    frame.set_edgecolor("black")  # visible border
    frame.set_linewidth(1.5)      # normal border thickness
    frame.set_alpha(0.2) 
    axes[0].axis('equal')
    leg = axes[1].legend( bbox_to_anchor=(0.7, 1.1))
    frame = leg.get_frame()
    frame.set_facecolor("none")   # transparent background
    frame.set_edgecolor("black")  # visible border
    frame.set_linewidth(1.5)      # normal border thickness
    frame.set_alpha(0.2)
    plt.savefig(r'C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Old Python Files\All 3 bodies LOOKING GOODDDD.png', transparent=True)
    avg = np.mean(R1)
    return avg
        
def dtvaryadv():
    #times = np.array([20, 30, 40, 50, 60])
    times = np.arange(720,100_000,500)
    Y = len(times)
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    Avgs = np.zeros(Y)
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for i in range(0,Y):
        dt = datetime.timedelta(minutes=int(times[i])).total_seconds()
        T = 365 * 60 * 60 * 24 * 100
        TotT = T / dt
        Vx, Vy, X, Y, R1, R2, R3, Em, Jm, Sm, t = BigFunc(TotT, dt)
        x = np.linspace(0, T, len(R1))
        p = times/(60*24)
        #axes[1].plot(x , R1,label = f'dt = {p} mins')
        #axes[1].plot(x , R2,label = f'dt = {p} mins')
        Avgs[i] = Em
    axes[1].loglog(p, Avgs)
    
    Vx, Vy, X, Y, R1, R2, R3, Em, Jm, Sm, t = BigFunc(oTotT, odt)
    axes[0].plot(X[:,0], Y[:,0], label = 'Sun')
    axes[0].plot(X[:,1], Y[:,1], label = 'Earth')
    axes[0].plot(X[:,2], Y[:,2], label = 'Jupiter')
    #axes[0].scatter(*COMCalc(), color='red', marker='x', label='COM')
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
    axes[1].set_xlabel('Timestep (Days)')
    #leg = axes[1].legend(loc='upper center', bbox_to_anchor=(0.7, 1.1))
    #frame = leg.get_frame()
    #frame.set_facecolor("none")   # transparent background
    #frame.set_edgecolor("black")  # visible border
    #frame.set_linewidth(1.5)      # normal border thickness
    #frame.set_alpha(0.2)
    plt.show()
    plt.savefig(r'C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Old Python Files\Earth Orbit Deviation.png', transparent=True)
    return Avgs

def EnergyPlot():
    fig, axes = plt.subplots(2, 1, figsize=(15, 5))
    ms = array_all[:,1]
    dt = datetime.timedelta(minutes=int(60)).total_seconds()
    T = 365 * 60 * 60 * 24
    TotT = T / dt
    Vx, Vy, X, Y, R1, R2, R3, Em, Jm, Sm, t = BigFunc(TotT, dt)
    halfposX = np.empty((len(X)-1, N))
    halfposY = np.empty((len(X)-1, N))
    #Uarr = np.empty(N*(N-1)/2)
    Utotarr = np.empty(len(halfposX))
    KEtotarr = np.empty(len(halfposX))
    
    for i in range (0,len(X)-1):
        Utot = 0
        for j in range (0,N):
            halfposX[i,j] = (X[i,j] + X[i+1,j]) / 2
            halfposY[i,j] = (Y[i,j] + Y[i+1,j]) / 2
        for l in range (0,N):
            for k in range(l+1,N):
                u = -G * (ms[l]*ms[k])/np.sqrt((halfposX[i,l] - halfposX[i,k])**2 + (halfposY[i,l] - halfposY[i,k])**2)
                Utot = Utot + u
        Utotarr[i] = Utot
        Ktot = 0
        for j in range (0,N):
            k = 1/2 * ms[j] * (Vx[i+1,j]**2 + Vy[i+1,j]**2)
            Ktot = Ktot + k
        KEtotarr[i] = Ktot
    Uavg = np.mean(Utotarr)
    Kavg = np.mean(KEtotarr)
    du = np.empty(len(halfposX))
    dk = np.empty(len(halfposX))
    for i in range (0, len(Utotarr)):
        du[i] = Utotarr[i] - Uavg
        dk[i] = KEtotarr[i] - Kavg
        
    x = np.linspace(0, T, len(Utotarr))
    #axes[0] = fig.add_axes([0.1, 0.35, 0.85, 0.6])
    fig = plt.figure(figsize=(12, 5))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])  # top 3x taller

    axes[0] = fig.add_subplot(gs[0])
    axes[1] = fig.add_subplot(gs[1])
    axes[0].plot(x, Utotarr)
    axes[0].plot(x, KEtotarr)
    axes[0] .xaxis.set_visible(False)
    axes[1].plot(x,du)
    axes[1].plot(x,dk)
    
    #axes[0].plot(X[:,0], Y[:,0], label = 'Sun')
    #axes[0].plot(X[:,1], Y[:,1], label = 'Earth')
    #axes[0].scatter(*COMCalc(), color='red', marker='x', label='COM')
    #axes[0].set_xlabel('x position (m)')
    #axes[0].set_ylabel('y position (m)')
    #axes[0].legend(loc='upper center', bbox_to_anchor=(0.8, 1.1))
    #leg = axes[0].legend(loc='upper center', bbox_to_anchor=(0.8, 1.1))
    #frame = leg.get_frame()
    #frame.set_facecolor("none")   # transparent background
    #frame.set_edgecolor("black")  # visible border
    #frame.set_linewidth(1.5)      # normal border thickness
    #frame.set_alpha(0.2) 
    #axes[0].axis('equal')
    #    for j in range (0,N):
    return Uavg, Kavg
            
            
    #return(halfposX,halfposY, TotT)
#U,K = EnergyPlot()
#print(U,K)

    
dtvary()
#A = dtvaryadv()
times = np.arange(720,100_000,500)
arr = np.array([7])
times = np.arange(720,100_000,500)
Y = len(times)
print(Y)





#DO DTS DIFFERENCES BY A FACTOR OF 2
#DO A PLOT OF EARTH COM DISTANCE WITH RESPECT TO JUPITER
#add a zoom in of the sun orbit on my earth jupiter plot
#
#
    
    