# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 20:48:36 2025

@author: adidu
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
import pandas as pd
import datetime

def drawing():
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8, 8))

# Positions of the two bodies
    body1 = (0, 0)
    body2 = (2, 2)

# Plot the two bodies
    ax.scatter(*body1, s=800, color='black', label='Mass m₁')
    ax.scatter(*body2, s=800, color='orange', label='Mass m₂')

# Draw the line connecting them (gravitational/spring line)
    ax.plot([body1[0], body2[0]], [body1[1], body2[1]], 'k:', linewidth=2)

# Annotate the masses
    ax.text(body1[0] - 0.5, body1[1], r"$m_1$", ha='center', fontsize=14,  color='black')
    ax.text(body2[0] + 0.5, body2[1], r"$m_2$", ha='center', fontsize=14,  color='orange')

# Annotate the distance
    mid_x = (body1[0] + body2[0]) / 2
    mid_y = (body1[1] + body2[1]) / 2
    ax.text(mid_x, mid_y + 0.3, r"$r$", ha='center', fontsize=14)

# Compute center of mass (for equal masses it's at x=0)
    com_x = (body1[0] + body2[0]) / 2
    com_y = (body1[1] + body2[1]) / 2

# Plot center of mass
    ax.scatter(com_x, com_y, color='black', s=100, marker='x', label='Center of Mass')
    ax.text(com_x, com_y - 0.4, r"COM", ha='center', fontsize=12)

# ---- Add vectors ----
    arrow_scale = 0.8

# Velocity vectors (tangential)
    #ax.arrow(body1[0], body1[1], 0, -arrow_scale, head_width=0.2, color='green', length_includes_head=True)
    #ax.arrow(body2[0], body2[1], 0, arrow_scale, head_width=0.2, color='green', length_includes_head=True)
    #ax.text(body1[0] - 0.4, body1[1] - 0.6, r"$\vec{v}_1$", fontsize=12, color='green')
    #ax.text(body2[0] + 0.4, body2[1] + 0.6, r"$\vec{v}_2$", fontsize=12, color='green')

# Force vectors (toward each other)
    ax.arrow(body1[0], body1[1], (body1[0] + body2[0]) / 2 - 0.25, (body1[1] + body2[1]) / 2 - 0.25, head_width=0.1, color='white', length_includes_head=True)
    ax.arrow(body2[0], body2[1], -((body1[0] + body2[0])/ 2 - 0.25),-((body1[1] + body2[1]) / 2 - 0.25), head_width=0.1, color='white', length_includes_head=True)
    ax.text(body1[0] + 0.5, body1[1] + 0.2, r"$\vec{F}_{12}$", fontsize=12, color='white')
    ax.text(body2[0] - 0.65, body2[1] - 0.3, r"$\vec{F}_{21}$", fontsize=12, color='white')
    
# ---- Formatting ----
    #ax.axhline(0, color='gray', linewidth=0.5)
    #ax.axvline(0, color='gray', linewidth=0.5)
    ax.set_xlim(-1, 4)
    ax.set_ylim(-3, 3)
    ax.set_aspect('equal', 'box')
    ax.axis('off')

    #ax.set_title("Two-Body System: Velocities and Forces", fontsize=15)
    #ax.legend(loc='upper right', labelspacing=2)
    plt.savefig(r'C:\Users\adidu\Documents\Work stuff\Year 3\Computing Project\Old Python Files\Force Diagram.png', transparent=True)
    plt.show()
    
drawing()