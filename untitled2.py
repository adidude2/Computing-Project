# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 20:36:38 2026

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


def Toomreheatmap(X, Y, Vx, Vy, SavedSteps, Rmax, mass, filename, t):
    fig, ax = plt.subplots(subplot_kw={'projection':'polar'}, figsize=(7,7))
    
    Nr   = 12
    Nphi = 16
    
    R_edges   = np.linspace(0, Rmax, Nr+1)
    phi_edges = np.linspace(-np.pi, np.pi, Nphi+1)
    R_centers = 0.5*(R_edges[:-1] + R_edges[1:])
    
    def update(frame):
        ax.cla()  # clear axes

        x = X[frame]
        y = Y[frame]
        vx   = Vx[frame]
        vy   = Vy[frame]
        R    = np.sqrt(x**2 + y**2)
        phi  = np.arctan2(y, x)
        R_safe = np.where(R == 0, 1e-10, R)

        Vr   = (x*vx + y*vy)/R_safe
        vphi = (x*vy - y*vx)/R_safe
        
        
        
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
        
        Phi, RR = np.meshgrid(phi_edges, R_edges)        

        pcm = ax.pcolormesh(Phi, RR/Kpc, Q, cmap='inferno', shading='auto', vmin=0, vmax=3)

        if frame == 0:
            fig.colorbar(pcm, ax=ax, label='Toomre Q')
        
        ax.set_ylabel("Radius (kpc)")

        # ax.set_rlabel_position(135)
        ax.set_title("Toomre Q Disk Map")
        time_gyr = (SavedSteps[frame] * t) / (1e6 * 3600*24*365)
        ax.set_title(f"Toomre Q Disk Map - Time (Myr) = {time_gyr:.3f}")
        
    ani = FuncAnimation(
        fig,
        update,
        frames = range(0, len(X), 5),
        interval=20  # ms between frames
        )
    
    
    
    ani.save(r"C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Old Python Files\\"
             + filename, fps=30)
    
    return ani
    
    
    