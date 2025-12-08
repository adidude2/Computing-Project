# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 11:07:00 2025

@author: adidu
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import datetime

# --- Load data (unchanged, assuming same columns as original) ---
df = pd.read_excel(r"C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Stuff.xlsx")
array_all = df.to_numpy()
Particles = df['Particle Number (N)'].to_numpy()
N = len(Particles)

# --- Constants ---
G = 6.6726e-11            # gravitational constant, N m^2 / kg^2
EPS = 1.0                 # softening length (m) - keep small, tune to avoid singularities

# --- Utility functions ---

def ForceCalc(m1, m2, x1, y1, x2, y2, eps=EPS):
    """
    Pairwise gravitational force from body 2 on body 1.
    Returns Fx, Fy (signed forces on body1).
    """
    dx = x2 - x1
    dy = y2 - y1
    r2 = dx*dx + dy*dy + eps*eps
    r3 = r2 * np.sqrt(r2)
    Fx = G * m1 * m2 * dx / r3
    Fy = G * m1 * m2 * dy / r3
    return Fx, Fy

def ForceIter(arrx, arry, idx, eps=EPS):
    """
    Compute NxN matrix of pairwise forces at timestep index `idx`.
    arrx, arry are 2D arrays (time, particle).
    Returns Fx_matrix, Fy_matrix where Fx[i,j] is force on i due to j.
    """
    Fx = np.zeros((N, N), dtype=float)
    Fy = np.zeros((N, N), dtype=float)
    for i in range(N):
        xi = arrx[idx, i]
        yi = arry[idx, i]
        for j in range(i+1, N):
            xj = arrx[idx, j]
            yj = arry[idx, j]
            mi = array_all[i, 1]
            mj = array_all[j, 1]
            fxi, fyi = ForceCalc(mi, mj, xi, yi, xj, yj, eps=eps)
            # force on i by j is +, on j by i is - (Newton III)
            Fx[i, j] = fxi
            Fy[i, j] = fyi
            Fx[j, i] = -fxi
            Fy[j, i] = -fyi
    return Fx, Fy

def AccelCalc(arrx, arry, idx, eps=EPS):
    """
    Compute accelerations (ax, ay) for all particles at time index `idx`.
    Returns arrays ax (m/s^2), ay (m/s^2)
    """
    Fx, Fy = ForceIter(arrx, arry, idx, eps=eps)
    SumFx = np.sum(Fx, axis=1)
    SumFy = np.sum(Fy, axis=1)
    ax = SumFx / array_all[:, 1]  # a = F / m
    ay = SumFy / array_all[:, 1]
    return ax, ay

def COMCalc(arrx, arry, idx):
    """
    Centre of mass at time index `idx`.
    """
    ms = array_all[:, 1]
    x = arrx[idx, :]
    y = arry[idx, :]
    M = np.sum(ms)
    COMx = np.sum(ms * x) / M
    COMy = np.sum(ms * y) / M
    return COMx, COMy

# --- Initial velocity helper(s) ---

def circular_velocity(central_index, positions_x, positions_y):
    """
    Compute circular orbital speed around `central_index` mass for each other body,
    using v = sqrt(G*M_center/r). Returns vx, vy across all bodies such that the
    velocity is perpendicular to radius vector (right-hand rotated).
    central_index: index of central mass (e.g. Sun = 0).
    positions_x, positions_y: 1D arrays of initial positions.
    """
    M_center = array_all[central_index, 1]
    vx = np.zeros(N)
    vy = np.zeros(N)
    cx = positions_x[central_index]
    cy = positions_y[central_index]
    for i in range(N):
        if i == central_index:
            continue
        rx = positions_x[i] - cx
        ry = positions_y[i] - cy
        r = np.hypot(rx, ry)
        if r == 0:
            continue
        speed = np.sqrt(G * M_center / r)
        # choose perpendicular direction (rotate +90°)
        # unit radial vector (rx/r, ry/r); perpendicular is (-ry/r, rx/r)
        ux = -ry / r
        uy = rx / r
        vx[i] = speed * ux
        vy[i] = speed * uy
    return vx, vy

# Your original empirical 'firstvel' approach left intact as alternative:
def firstvel_empirical(arrx0, arry0):
    """
    Original empirical approach (kept for compatibility).
    arrx0, arry0: arrays with initial positions (1D)
    Returns speeds (1D array).
    """
    ms = array_all[:, 1]
    xs = arrx0
    comx = np.sum(ms * xs) / np.sum(ms)
    # get forces from initial positions:
    # build temporary 2D arrays with only row 0 set
    tmpX = np.zeros((1, N))
    tmpY = np.zeros((1, N))
    tmpX[0, :] = arrx0
    tmpY[0, :] = arry0
    Fx, Fy = ForceIter(tmpX, tmpY, 0)
    SumFx = np.sum(Fx, axis=1)
    speeds = np.zeros(N)
    for i in range(N):
        val = (xs[i] - comx) * SumFx[i] / ms[i]
        speeds[i] = np.sqrt(abs(val)) if val != 0 else 0.0
    # preserve the odd sign choice from original if you want:
    speeds[0] = -speeds[0]
    return speeds

# --- Main integrator (BigFunc) ---

def BigFunc(n_steps, dt_seconds, init_vel_method='circular', central_index=0, eps=EPS):
    """
    Integrate the system for n_steps (int) with timestep dt_seconds (float).
    init_vel_method: 'circular' or 'empirical'
    central_index: index used by circular velocity initializer (default 0 -> Sun)
    Returns: VX, VY, X, Y arrays of shape (n_steps, N) and Eradius, Jradius, Sradius arrays (n_steps-1)
    """
    n_steps = int(n_steps)
    X = np.zeros((n_steps, N), dtype=float)
    Y = np.zeros((n_steps, N), dtype=float)
    VX = np.zeros((n_steps, N), dtype=float)
    VY = np.zeros((n_steps, N), dtype=float)

    # initial positions from the input file
    X[0, :] = array_all[:, 2]
    Y[0, :] = array_all[:, 3]

    # initial velocities
    if init_vel_method == 'circular':
        init_vx, init_vy = circular_velocity(central_index, X[0, :], Y[0, :])
    else:
        speeds = firstvel_empirical(X[0, :], Y[0, :])
        # orient speeds tangentially like original InitVelCalc tried to do
        # small-angle assumption preserved: put all motion perpendicular to radius vector
        init_vx = np.zeros(N)
        init_vy = np.zeros(N)
        cx, cy = COMCalc(X, Y, 0)
        for i in range(N):
            rx = X[0, i] - cx
            ry = Y[0, i] - cy
            r = np.hypot(rx, ry)
            if r == 0:
                continue
            ux = -ry / r
            uy = rx / r
            init_vx[i] = speeds[i] * ux
            init_vy[i] = speeds[i] * uy

    VX[0, :] = init_vx
    VY[0, :] = init_vy

    # preallocate diagnostic arrays
    Eradius = np.zeros(n_steps - 1) if n_steps > 1 else np.array([])
    Jradius = np.zeros(n_steps - 1) if n_steps > 1 else np.array([])
    Sradius = np.zeros(n_steps - 1) if n_steps > 1 else np.array([])

    # initial centre-of-mass for reference
    A1, B1 = COMCalc(X, Y, 0)
    xs_init = array_all[:, 2]

    for i in range(n_steps - 1):
        # accelerations at step i
        ax, ay = AccelCalc(X, Y, i, eps=eps)   # m/s^2

        # diagnostics (example: radius relative to initial)
        # NOTE: this keeps the same indexing used in your original code:
        Eradius[i] = np.hypot(X[i, 1] - COMCalc(X, Y, i)[0], Y[i, 1] - COMCalc(X, Y, i)[1])  # Earth
        Jradius[i] = (xs_init[2] - A1) - np.hypot(X[i, 2] - COMCalc(X, Y, i)[0], Y[i, 2] - COMCalc(X, Y, i)[1])  # Jupiter-ish
        Sradius[i] = (xs_init[0] - A1) - np.hypot(X[i, 0] - COMCalc(X, Y, i)[0], Y[i, 0] - COMCalc(X, Y, i)[1])  # Sun-ish

        # update velocities and positions (explicit)
        VX[i+1, :] = VX[i, :] + ax * dt_seconds
        VY[i+1, :] = VY[i, :] + ay * dt_seconds
        X[i+1, :] = X[i, :] + VX[i+1, :] * dt_seconds
        Y[i+1, :] = Y[i, :] + VY[i+1, :] * dt_seconds

    # normalize Eradius by initial Earth x (preserve your original return format)
    if len(Eradius) > 0:
        Eradius = Eradius / xs_init[1]

    return VX, VY, X, Y, Eradius, Jradius, Sradius

# --- Example helper to convert human times and run ---
def hours_to_seconds(h):
    return datetime.timedelta(hours=h).total_seconds()

def run_example_one_year(dt_hours=1.0):
    dt = hours_to_seconds(dt_hours)
    seconds_per_year = 365 * 24 * 60 * 60
    n_steps = int(seconds_per_year / dt)
    return BigFunc(n_steps, dt, init_vel_method='circular', central_index=0, eps=EPS)

# --- Small plotting helper (reworked dtvary) ---
def dtvary():
    times_hours = np.array([12], dtype=float)   # example list of dt in hours
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    cmap = LinearSegmentedColormap.from_list("orange_black", ["orange", "black"])
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for i, th in enumerate(times_hours):
        dt = hours_to_seconds(th)
        T_seconds = 365 * 24 * 60 * 60 * 50  # 20 years in seconds
        n_steps = int(T_seconds // dt)
        VX, VY, X, Y, R1, R2, R3 = BigFunc(n_steps, dt, init_vel_method='circular', central_index=0)
        x = np.linspace(0, T_seconds, len(R1)) / (60*60*24*365)  # years
        axes[1].plot(x, R1, label=f'dt = {int(th)} hours')
        axes[0].plot(X[:, 0], Y[:, 0], label='Sun')
        axes[0].plot(X[:, 1], Y[:, 1], label='Earth')
        axes[0].plot(X[:, 2], Y[:, 2], label='Jupiter')
    axes[0].set_xlabel('x position (m)')
    axes[0].set_ylabel('y position (m)')
    axes[0].legend()
    axes[0].axis('equal')
    plt.show()

dtvary()
