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
    return ax, ay, az, Menc, Rs, rho0

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

def AccelCalc(arrx,arry, arrz, mass, t, Mhalo, Rvir, c):
    Fx, Fy, Fz = force_vectorised(arrx, arry, arrz, mass)
    m = mass
    ax = Fx / m
    ay = Fy / m
    az = Fz / m
    
    
    hax, hay, haz, M, Rs, rho0 = HaloAccel(
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
    return dvx, dvy, dvz, M

def COMCalc(arrx,arry, arrz, mass):
    ms = mass
    COMx = (np.sum(ms * arrx))/(np.sum(ms))
    COMy = (np.sum(ms * arry))/(np.sum(ms))
    COMz = (np.sum(ms * arrz))/(np.sum(ms))
    return COMx, COMy, COMz

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

def BetterHaloVels(a, t, dvx, dvy, dvz, mass, f, fz):
    X = a[:,0]
    Y = a[:,1]
    Z = a[:,2]
    L =len(X)

    
    ax = dvx / t
    ay = dvy / t
    az = dvz / t
    
    r = np.sqrt(X**2 + Y**2)


    ar = (ax * X + ay * Y) / r
    
    vtrue2 = r * (-ar)
    vtrue2 = np.where(vtrue2 > 0, vtrue2, 0.0)
    vtan = np.sqrt(vtrue2)
    
    Rhat = np.column_stack((X,Y)) / r[:, None]
    
    vrad = np.random.normal(0, f * vtan)
    
    vz = np.random.normal(0, fz * vtan)
    
    anew = np.column_stack((-Y,X))
    norms = np.linalg.norm(anew, axis=1)
    Vhat = anew / norms[:, None]
    
    velsxy = vtan[:, None] * Vhat + vrad[:, None] * Rhat
    vels = np.column_stack((velsxy,vz))
    
    return vels

def HugeFunc(T, t, ass, mass, Mhalo, Rvir, c, f, fz):
    a = ass
    Xmatrix = a[:,0]
    Ymatrix = a[:,1]
    Zmatrix = a[:,2]
    
    hax, hay, haz, Menc, Rs, rho0 = HaloAccel(
       Xmatrix,
       Ymatrix,
       Zmatrix,
       Mhalo,
       Rvir,
       c
    )
    
    dvx, dvy, dvz, Menc = AccelCalc(
        Xmatrix, 
        Ymatrix, 
        Zmatrix, 
        mass, 
        t, 
        Mhalo, 
        Rvir, 
        c)

    
    #vel = HaloVels(a, Menc, mass)
    vel = BetterHaloVels(a, t, dvx, dvy, dvz, mass, f, fz)   
    
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
    SavedMenc = []
    for i in range (0,int(T-1)):
        A,B,C = COMCalc(Xmatrix,Ymatrix, Zmatrix, mass)
        if i % 10 == 0:
            SavedCOM.append((A, B, C))
            SavedX.append(Xmatrix.copy())
            SavedY.append(Ymatrix.copy())
            SavedZ.append(Zmatrix.copy())
            SavedVX.append(VXmatrix.copy())
            SavedVY.append(VYmatrix.copy())
            SavedVZ.append(VZmatrix.copy())
            SavedMenc.append(Menc.copy())
            SavedSteps.append(i)
            
        dvx, dvy, dvz, Menc = AccelCalc(Xmatrix, Ymatrix, Zmatrix, mass, t, Mhalo, Rvir, c)
        
        # Original Integrator
        # VXmatrix += dvx
        # VYmatrix += dvy
        # VZmatrix += dvz

        # Xmatrix += VXmatrix * t
        # Ymatrix += VYmatrix * t
        # Zmatrix += VZmatrix * t
        
        # Velocity Verlet (symplectic): keeps disk from collapsing
        # x_{n+1} = x_n + v_n*dt + 0.5*a_n*dt^2
        Xmatrix += VXmatrix * t + 0.5 * dvx * t
        Ymatrix += VYmatrix * t + 0.5 * dvy * t
        Zmatrix += VZmatrix * t + 0.5 * dvz * t

        # a_{n+1} at new position
        dvx_new, dvy_new, dvz_new, Menc = AccelCalc(Xmatrix, Ymatrix, Zmatrix, mass, t, Mhalo, Rvir, c)
        # v_{n+1} = v_n + 0.5*(a_n + a_{n+1})*dt
        VXmatrix += 0.5 * (dvx + dvx_new)
        VYmatrix += 0.5 * (dvy + dvy_new)
        VZmatrix += 0.5 * (dvz + dvz_new)
    return SavedCOM, SavedX, SavedY, SavedZ, SavedVX, SavedVY, SavedVZ, mass, SavedSteps, Rs, rho0

def plotfig(R, SavedCOM, SavedX, SavedY, SavedZ, SavedSteps, t, view, filename):
    
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
        
        Kpc = 3.086e+19
        
        L = 1 * R / Kpc
        ax.set_xlim(-L, L)
        ax.set_ylim(-L, L)
        ax.set_zlim(-L, L)
        ax.set_xticks([-L, 0, L])
        ax.set_yticks([-L, 0, L])
        ax.set_zticks([-L, 0, L])
        ax.set_xlabel("X (kpc)")
        ax.set_ylabel("Y (kpc)")
        ax.set_zlabel("Z (kpc)")
        #ax.set_xticks([])
        #ax.set_yticks([])
        #ax.set_zticks([])
        #ax.set_axis_off()
        

        ax.scatter(
            X / Kpc,
            Y / Kpc,
            Z / Kpc,
            s=2
        )
        if view == "face":
            ax.view_init(elev=90, azim=90)
        elif view == "edge":
            ax.view_init(elev=0, azim=90)
        elif view == "inclined":
            ax.view_init(elev=30, azim=45)
        time_gyr = (SavedSteps[frame] * t) / (1e6 * 3600*24*365)
        ax.set_title(f"Time (Myr) = {time_gyr:.3f}")
    
    ani = FuncAnimation(
        fig,
        update,
        frames=len(SavedX),
        interval=20  # ms between frames
        )
    
    
    ani.save(r"C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Old Python Files\\"
             + filename, fps=30)


    plt.show()
    
    return ani

def EnergyPlotfrac(t, R, Mhalo, Rvir, c, SavedCOM, SavedX, SavedY, SavedZ, Vx, Vy, Vz, mass, SavedSteps, Rs, rho0):
    fig, axes = plt.subplots(3, 1, figsize=(9, 10))
    #fig = plt.figure(figsize=(7, 5))
    #ax = fig.add_subplot()
    ani_face = plotfig(R, SavedCOM, SavedX, SavedY, SavedZ, SavedSteps, t,
                     view = "face",
                     filename="disk_faceon.mp4")
    ani_edge = plotfig(R, SavedCOM, SavedX, SavedY, SavedZ, SavedSteps, t,
                     view = "edge",
                     filename="disk_edgeon.mp4")
    #ani_inclined = plotfig(R, SavedCOM, SavedX, SavedY, SavedZ, SavedSteps, t,
    #                 view = "inclined",
    #                 filename="disk_incline.mp4")
    
    ms = mass
    X = SavedX
    Y = SavedY
    Z = SavedZ
    Np = len(X[0])
    #Uarr = np.empty(N*(N-1)/2)
    Utotarr = np.empty(len(X)-1)
    KEtotarr = np.empty(len(X)-1)
    Virtotarr = np.empty(len(X)-1)
    Lztot = np.empty(len(X)-1)
    for i in range (0,len(X)-1):
        x = X[i]
        y = Y[i]
        z = Z[i]
        #dx = x[:,None] - x[None,:]
        #dy = y[:,None] - y[None,:]
        #dz = z[:,None] - z[None,:]
        #r = np.sqrt(dx**2 + dy**2 + dz**2 + e**2)
        Utot = 0
        midX = 0.5 * (X[i] + X[i+1])
        midY = 0.5 * (Y[i] + Y[i+1])
        midZ = 0.5 * (Z[i] + Z[i+1])
        
        positions = np.stack([midX, midY, midZ], axis=1)
        #positions = np.stack([midX, midY, midZ], axis=1)
        d = positions[:,None,:] - positions[None,:,:]
        r = np.sqrt(np.sum(d*d,axis=2)+e**2)
        U = -G * (ms[:,None]*ms[None,:]) / r
        np.fill_diagonal(U,0)
        
        r = np.sqrt(x**2 + y**2 + z**2)
        r_safe = np.where(r==0, 1e-10, r)
        
        phihalo = (-4 * np.pi * G * rho0 * Rs**3 / r_safe) * np.log(
            
             1 + r_safe/Rs
            
             )
        
        Uhalo   = mass * phihalo
        
        U_stellar = np.sum(U) / 2
        
        Utot = U_stellar + np.sum(Uhalo)
        Utotarr[i] = Utot
        
        hax, hay, haz, _, _, _ = HaloAccel(x, y, z, Mhalo, Rvir, c)
        Omega_halo = np.sum(ms * (x * hax + y * hay + z * haz))
        
        vx = Vx[i+1]
        vy = Vy[i+1]
        vz = Vz[i+1]
        v2 = vx**2 + vy**2 + vz**2
        Ktot = 0.5 * np.sum(ms * v2)
        KEtotarr[i] = Ktot
        
        Vir = 0
        # Vir = 2K + Omega_stellar + Omega_halo; for gravity Omega_stellar = U_stellar
        Vir = 2 * Ktot + U_stellar + Omega_halo
        Virtotarr[i] = Vir
        
        Lz = 0
        Lz = np.sum(ms*(x*vy-y*vx))
        Lztot[i] = Lz

    Etot = Utotarr + KEtotarr
    Uavg = np.mean(Utotarr)
    Kavg = np.mean(KEtotarr)
    Tavg = np.mean(Etot)
    Vavg = np.mean(Virtotarr)
    du = np.empty(len(Utotarr))
    dk = np.empty(len(Utotarr))
    dTot = np.empty(len(Utotarr))
    dLz = np.empty(len(Utotarr))
    dVir = np.empty(len(Utotarr))
    for i in range (0, len(Utotarr)):
        du[i]   = (Utotarr[i] - Utotarr[0])/abs(Utotarr[0])
        dk[i]   = (KEtotarr[i] - KEtotarr[0])/abs(KEtotarr[0])
        dTot[i] = (Etot[i] - Etot[0])/abs(Etot[0])
        dVir[i] = np.abs(Virtotarr[i]) / np.abs(Utotarr[0])
        dLz[i]  = (Lztot[i] - Lztot[0])/abs(Lztot[0])
    
    max_du   = np.max(np.abs(du))
    max_dk   = np.max(np.abs(dk))
    max_de   = np.max(np.abs(dTot))
    max_dvir = np.max(np.abs(dVir))
    max_dlz  = np.max(np.abs(dLz))
        
    E0   = Etot[0]
    de   = 100 * (Etot - E0) / E0
    
    mean_de = np.mean(np.abs(de))
    
    
    

    #mean_dvir = 100 * np.mean(vir_ratio)
    
    seconds_per_year = 1e6 * 3600 * 24 * 365
    x = (np.array(SavedSteps[:-1]) * t) / seconds_per_year
    axes[0].plot(x, Utotarr, label = '$U$', color = 'green')
    axes[0].plot(x, KEtotarr, label = '$K$', color = 'orange')
    axes[0].plot(x, Etot, label = '$E$', color = 'purple')
    axes[0].plot(x, Virtotarr, label = '$2K + U$', color = 'tab:blue')
    axes[1].plot(x,du, label = r'$\Delta U/U_0$', color = 'green')
    axes[1].plot(x,dk, label = r'$\Delta K/K_0$', color = 'orange')
    axes[1].plot(x,dTot, label = r'$\Delta E/E_0$', color = 'purple')
    axes[1].plot(x,dVir, label = r'$(2K+U)/|U_0|$', color = 'tab:blue')
    axes[2].plot(x, dLz, label = r'$\Delta L_z/L_{z0}$', color = 'black')
    axes[0].tick_params(axis='x', labelbottom=False)
    axes[1].tick_params(axis='x', labelbottom=False)
    axes[2].set_xlabel("Time (Myrs)")
    axes[0].set_ylabel("Energy (J)")
    axes[1].set_ylabel(r'$\frac{\Delta E}{|E_0|}$', rotation='horizontal',  fontsize=14, labelpad=10)
    axes[2].set_ylabel(r'$\frac{\Delta L_z}{|L_{z0}|}$', rotation='horizontal',  fontsize=14, labelpad=14)
    axes[0].legend(fontsize=9)
    axes[1].legend(fontsize=9)
    #axes[2].legend()
    fig.savefig(r'C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Old Python Files\Energy Analysis.png', transparent=True)
    return Uavg, Kavg, Tavg, max_du, max_dk, max_de, max_dvir, max_dlz, ani_face, ani_edge, fig

def Toomre(X, Y, Vx, Vy, Rmax, mass, filename):
    x = X
    y = Y
    vx = Vx[0]
    vy = Vy[0]
    R = np.sqrt(x**2 + y**2)
    #phi = np.arctan2(Y[0], X[0])
    
    Vr = (x*vx + y*vy)/R
    vphi = (x*vy - y*vx)/R
    
    
    nbins = 15
    R_bins = np.linspace(0, Rmax, nbins+1)
    R_centers = 0.5*(R_bins[:-1] + R_bins[1:])
    
    Sigma = np.zeros(nbins)
    sigma_R = np.zeros(nbins)
    vphi_mean = np.zeros(nbins)
    
    for i in range(nbins):
        
        mask = (R >= R_bins[i]) & (R < R_bins[i+1])
        
        Mbin = np.sum(mass[mask])
        
        area = np.pi * (R_bins[i+1]**2 - R_bins[i]**2)
        
        Sigma[i] = Mbin / area

        if np.sum(mask) > 1:
            sigma_R[i] = np.std(Vr[mask])

        if np.sum(mask) > 0:
            vphi_mean[i] = np.mean(vphi[mask])
    
    Omega = vphi_mean / R_centers
    Omega2 = Omega**2
    dOmega2_dR = np.gradient(Omega2, R_centers)

    kappa = np.sqrt(4*Omega2 + R_centers*dOmega2_dR)
    
    Q = sigma_R * kappa / (3.36 * G * Sigma)
    
    plt.figure(figsize=(7, 5))
    plt.plot(R_centers/Kpc, Q)
    plt.axhline(1, linestyle='--')
    plt.axhline(2, linestyle=':')
    plt.xlabel('Radius (Kpc)')
    plt.ylabel('Toomre Q parameter')
    plt.title("Toomre Q against radius")
    plt.savefig(r'C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Energy diagnostic plots\Toomre Analysis\\'
                + filename, transparent=True)
    
    
def Toomreheatmap(X, Y, Vx, Vy, Rmax, mass, filename):
    x    = X
    y    = Y
    vx   = Vx[0]
    vy   = Vy[0]
    R    = np.sqrt(x**2 + y**2)
    phi  = np.arctan2(y, x)
    Vr   = (x*vx + y*vy)/R
    vphi = (x*vy - y*vx)/R
    
    
    Nr   = 12
    Nphi = 16
    
    R_edges   = np.linspace(0, Rmax, Nr+1)
    phi_edges = np.linspace(-np.pi, np.pi, Nphi+1)
    R_centers = 0.5*(R_edges[:-1] + R_edges[1:])
    
    Sigma     = np.zeros((Nr,Nphi))
    sigma_R   = np.zeros((Nr,Nphi))
    vphi_mean = np.zeros(Nr)
    Q         = np.zeros((Nr,Nphi))
    
    
    for i in range(Nr):

        R_inner = R_edges[i]
        R_outer = R_edges[i+1]
        
        mask1 = (R >= R_inner) & (R < R_outer)

        for j in range(Nphi):

            phi_inner = phi_edges[j]
            phi_outer = phi_edges[j+1]

            mask2 = (
                (R >= R_inner) & (R < R_outer) &
                (phi >= phi_inner) & (phi < phi_outer)
                )
    
            A = 0.5 * (R_outer**2 - R_inner**2)*(phi_outer - phi_inner)
            Mbin = np.sum(mass[mask2])
            Sigma[i,j] = Mbin / A
            
            if np.sum(mask2) > 1:
                sigma_R[i,j] = np.std(Vr[mask2])
                
        if np.sum(mask1) > 0:
            vphi_mean[i] = np.mean(vphi[mask1])
                
    Omega = vphi_mean / R_centers
    Omega2 = Omega**2
    dOmega2_dR = np.gradient(Omega2, R_centers)

    kappa = np.sqrt(4*Omega2 + R_centers*dOmega2_dR)
    
    for i in range(Nr):
        for j in range(Nphi):
            
            if Sigma[i,j] > 0 and sigma_R[i,j] > 0:
                Q[i,j] = sigma_R[i,j] * kappa[i] / (3.36 * G * Sigma[i,j])
            else:
                Q[i,j] = np.nan
    
    Phi, R = np.meshgrid(phi_edges, R_edges)        
    fig, ax = plt.subplots(subplot_kw={'projection':'polar'}, figsize=(7,7))

    pcm = ax.pcolormesh(Phi, R/Kpc, Q, cmap='inferno', shading='auto', vmin=0, vmax=3)

    fig.colorbar(pcm, ax=ax, label='Toomre Q')
    
    ax.set_ylabel("Radius (kpc)")
    
    # ticks_kpc = [2,4,6,8,10,12]

    # ax.set_yticks(np.array(ticks_kpc) * Kpc)
    # ax.set_yticklabels(ticks_kpc)

    # ax.set_rlabel_position(135)
    ax.set_title("Toomre Q Disk Map")
    fig.savefig(r'C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Energy diagnostic plots\Toomre Analysis\\'
                + filename, transparent=True)
            



Msol            = 1.989e+30
Kpc             = 3.086e+19
Myr             = 3600 * 24 * 365 * 1e6

N               = 4000              # Number of bodies
Rd              = 1.3 * 2.15 * Kpc  # Scale length
Rmax            = 5 * Rd            # Truncation Radius
z0              = 0.3 * Kpc         # Thickness of the disk
Rvir            = 200 * Kpc         # Halo virial radius
Mstel           = 5.04e10 * Msol    # Total mass of Stars
Mhalo           = 0.97e12 * Msol    # Mass of Whole Halo

Mhalo_effective = 1.5 * Mhalo       # Effective halo mass in simulation
Mstel_effective = 0.9 * Mstel       # Effective stellar mass in simulation

e               = 0.3 * Kpc         # softening
c               = 9.4               # halo concentration
sigma_R_factor  = 0.30              # radial velocity dispersion
sigma_z_factor  = 0.1               # vertical velocity dispersion


odt = 0.1 * Myr                     # Timestep
T_orbit = 212 * Myr * 0.05             # Simulation time

oT = T_orbit
oTotT = oT / odt
print(oTotT)

a,b = Posallocate(N, Rd, Rmax, Mstel_effective, z0)
X = a[:,0]
Y = a[:,1]
Z = a[:,2]


#ani, SavedX = plotfig(oTotT,  odt, R, a, b)
SavedCOM, SavedX, SavedY, SavedZ, Vx, Vy, Vz, mass, SavedSteps, Rs, rho0 = HugeFunc(
   
    oTotT,
    odt,
    a,
    b,
    Mhalo_effective,
    Rvir,
    c,
    sigma_R_factor,
    sigma_z_factor
    
    )

# Uavg, Kavg, Tavg, max_du, max_dk, max_de, max_dvir, max_dlz, ani_face, ani_edge, fig = EnergyPlotfrac(
     
#     odt,
#     Rmax,
#     Mhalo_effective,
#     Rvir,
#     c,
#     SavedCOM,
#     SavedX,
#     SavedY,
#     SavedZ,
#     Vx,
#     Vy,
#     Vz,
#     mass,
#     SavedSteps,
#     Rs,
#     rho0
    
#     )

# print("Mean Kinetic:", Kavg)
# print("Mean Potential:", Uavg)
# print("Mean Total:", Tavg)
# print("Max Energy Fractional Difference:", max_de)
# print("Max Kinetic Fractional Difference:", max_dk)
# print("Max Potential Fractional Difference:", max_du)
# print("Max Virial Residue:", max_dvir)
# print("Max Angular Momentum Fractional Difference:", max_dlz)


#print(SavedX)
#rhor = HaloAccel(Xmatrix, Ymatrix, Zmatrix, Mhalo, Rvir, c)

filename1 = "Toomre radial plot (4000p).png"
filename2 = "Toomre Heatmap (4000p).png"

Toomre(X, Y, Vx, Vy, Rmax, b, filename1)
Toomreheatmap(X, Y, Vx, Vy, Rmax, b, filename2)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')
ax.set_box_aspect([1, 1, 1])
L = 1 * Rmax / Kpc
ax.set_xlim(-L, L)
ax.set_ylim(-L, L)
ax.set_zlim(-L, L)
ax.set_xticks([-L, 0, L])
ax.set_yticks([-L, 0, L])
ax.set_zticks([-L, 0, L])
ax.view_init(elev=90, azim=90)
ax.scatter(X/Kpc,Y/Kpc,Z/Kpc, s=3)
plt.show()

fig.savefig(r'C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Energy diagnostic plots\Toomre Analysis\\'
            + "4000p positions.png", transparent=True)


# """

# SHOW THE PLOT FROM ABOVE SINCE YOURE NOT EVEN SHOWING THE Z DIRECTION ANYWAY
# ALSO REMEMBER TO CHANGE THE INITIAL VELS CALC IN THE BIG FUNC ANDDDD THE ACCEL CALC WHEN IM INCULDING THE HALO AND NOT

# """