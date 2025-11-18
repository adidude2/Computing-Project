# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 21:06:03 2024

@author: adidu
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

#Define universal gravitation constant
G=6.6726e-11 #N-m2/kg2

#Reference quantities
mE = 5.9742e24
mM = 7.35e22
d = 3.844e8


#L1 point determining code
mew = (mM / (mM + mE))

roots = np.roots([1, mew-3, 3-(2*mew), -mew, 2*mew, -mew])
print(f'Roots = {roots[4]}')
print(type(roots[4]))
Y = roots[4]
A = Y.real
lr = (d * A)-15280797.926
print(f'Radius from the Moon = {lr}')
print(f'Radius from the Earth = {d-lr}')
ir = 7000e3 #d*(mM/(3 * mE))**(1/3) #+ 2991932.6395

alpha = np.pi*(1-(1/(2 * np.sqrt(2))) * np.sqrt(((ir/(d))+1)**3))
print(f'Alpha = {alpha}')


rE = d * (mM / (mM + mE))
rM = d * (mE / (mM + mE))
T = (((4*np.pi**2)/(G*(mE + mM))) * d**3)**0.5
#T = 2419200
print(T)
omega=2*np.pi/T

# create a time array up to t_max sampled at steps of dt
dt = 10
t_min = 0.
t_max =  3*T/4#(1 - ((alpha+0.2)/(2*np.pi))) * 5829.24657 #3/4 * 5829.24657 T/5 
t = np.arange(t_min,t_max,dt)
nstep = int((t_max-t_min)/dt)

# create a positions array for both massive bodies
x1 = np.empty([len(t)])
y1 = np.empty([len(t)])
x2 = np.empty([len(t)])
y2 = np.empty([len(t)])

# Model the movement of the 2 bodies using simple harmonic motion
for i in range(len(t)):
    x1[i]=np.cos(omega*t[i])*rE  # set the massive particle positions
    y1[i]=np.sin(omega*t[i])*rE 
    x2[i]=-np.cos(omega*t[i])*rM
    y2[i]=-np.sin(omega*t[i])*rM
    
    


# first postionas and speeds for the rocket
#multiplier for taylor ex = 1.04860411065608946
dr = 7000e3 #d*(mM/(3 * mE))**(1/3) #+ 2991932.6395
x= x1[0] - dr
y= 0
r=np.sqrt(x*x+y*y)
w = (2 * np.pi)/ T
vx = 0
vy= -np.sqrt((G*mE)/dr)      #2 * np.pi * dr/T # velocity needed to stay in the L2 point
#print(dr + 2991932.6395)
#print(dr + 3664097.11465)





xs = np.empty([len(t)+1])
ys = np.empty([len(t)+1])
dxs = np.empty([len(t)+1])
dys = np.empty([len(t)+1])
drs = np.empty([len(t)+1])
xs[0] = x
ys[0] = y
dxs[0] = vx
dys[0] = vy
plt.figure(figsize=(10,11))
def iniRungeK():
    dt = 10
    for i in range(len(t)):
        dE = np.sqrt((xs[i] - x1[i])**2 + (ys[i] - y1[i])**2)
        dM = np.sqrt((xs[i] - x2[i])**2 + (ys[i] - y2[i])**2)
        d2x = (-G * mE * (xs[i] - x1[i]) / (dE**3)) - (G * mM * (xs[i] - x2[i]) / (dM**3))
        d2y = (-G * mE * (ys[i] - y1[i]) / (dE**3)) - (G * mM * (ys[i] - y2[i]) / (dM**3))
        
        zx1 = xs[i] + (dt * dxs[i] / 2)
        dzx1 = dxs[i] + (dt * d2x / 2)
        zy1 = ys[i] + (dt * dys[i] / 2)
        dzy1 = dys[i] + (dt * d2y / 2)
        
        dE = np.sqrt((zx1 - x1[i])**2 + (zy1 - y1[i])**2)
        dM = np.sqrt((zx1 - x2[i])**2 + (zy1 - y2[i])**2)
        zx2 = xs[i] + (dt * dzx1 / 2)
        dzx2 = dxs[i] + ((-G * mE * (zx1 - x1[i]) / (dE**3)) - (G * mM * (zx1 - x2[i]) / (dM**3))) * (dt/2)
        zy2 = ys[i] + (dt * dzy1 / 2)
        dzy2 = dys[i] + ((-G * mE * (zy1 - y1[i]) / (dE**3)) - (G * mM * (zy1 - y2[i]) / (dM**3))) * (dt/2)
        
        dE = np.sqrt((zx2 - x1[i])**2 + (zy2 - y1[i])**2)
        dM = np.sqrt((zx2 - x2[i])**2 + (zy2 - y2[i])**2)
        zx3 = xs[i] + (dt * dzx2)
        dzx3 = dxs[i] + ((-G * mE * (zx2 - x1[i]) / (dE**3)) - (G * mM * (zx2 - x2[i]) / (dM**3))) * dt
        zy3 = ys[i] + (dt * dzy2)
        dzy3 = dys[i] + ((-G * mE * (zy2 - y1[i]) / (dE**3)) - (G * mM * (zy2 - y2[i]) / (dM**3))) * dt
        
        dE = np.sqrt((zx3 - x1[i])**2 + (zy3 - y1[i])**2)
        dM = np.sqrt((zx3 - x2[i])**2 + (zy3 - y2[i])**2)
        xs[i+1] = xs[i] + (dt/6) * (dxs[i] + (2 * dzx1) + (2 * dzx2) + dzx3)
        dxs[i+1] = dxs[i] + (dt/6) * (d2x + (2 * (dzx2 - dxs[i]) * (2/dt)) + (2 * (dzx3 - dxs[i]) * (1/dt)) + ((-G * mE * (zx3 - x1[i]) / (dE**3)) - (G * mM * (zx3 - x2[i]) / (dM**3))))
        ys[i+1] = ys[i] + (dt/6) * (dys[i] + (2 * dzy1) + (2 * dzy2) + dzy3)
        dys[i+1] = dys[i] + (dt/6) * (d2y + (2 * (dzy2 - dys[i]) * (2/dt)) + (2 * (dzy3 - dys[i]) * (1/dt)) + ((-G * mE * (zy3 - y1[i]) / (dE**3)) - (G * mM * (zy3 - y2[i]) / (dM**3))))
    
    #dist = np.sqrt((xs[-1] - x2[-1])**2 + (ys[-1] - y2[-1])**2)
    edist = np.sqrt((xs[-1] - x1[-1])**2 + (ys[-1] - y1[-1])**2)
    
    #return dist #use this line for optidr

    #print(f'Initial Distance From Earth = {dr}')
    #print(f'Distance From Earth after analysis = {edist}')
    #print(f'Difference in distance = {abs(edist-dr)}')
    #print(dxs[-1])
    #print(dys[-1])
    plt.plot(x1,y1, label = "Earth", color = 'blue')
    plt.plot(x2,y2, label = "Moon", color = 'orange')
    plt.plot(xs,ys, label = "Rocket Initial", color = 'red')
    plt.legend()
    
iniRungeK()   
 
r = np.sqrt((xs[-1] - x1[-1])**2 + (ys[-1] - y1[-1])**2)


Th = np.pi * np.sqrt(((ir+d-lr)**3)/(8 * G * mE))
# create a time array up to t_max sampled at steps of dt
dt = 10
t_min = t_max
t_max = t_max + Th #-15000 #T/500
t = np.arange(t_min,t_max,dt)
nstep = int((t_max-t_min)/dt)


x1 = np.empty([len(t)])
y1 = np.empty([len(t)])
x2 = np.empty([len(t)])
y2 = np.empty([len(t)])

# Model the movement of the 2 bodies using simple harmonic motion
for i in range(len(t)):
    x1[i]=np.cos(omega*t[i])*rE  # set the massive particle positions
    y1[i]=np.sin(omega*t[i])*rE 
    x2[i]=-np.cos(omega*t[i])*rM
    y2[i]=-np.sin(omega*t[i])*rM

dv = np.sqrt(G * mE/r) * (np.sqrt(2*(d-lr)/(r+(d-lr)))-1) 
theta = np.arctan(dys[-1]/dxs[-1])
if dys[-1] > 0 and dxs[-1] < 0:
    theta = np.pi + theta
if dys[-1] < 0 and dxs[-1] < 0:
    theta = np.pi + theta
xcs = np.empty([len(t)+1])
ycs = np.empty([len(t)+1])
dxcs = np.empty([len(t)+1])
dycs = np.empty([len(t)+1])
vs = np.empty([len(t)+1])
    #drcs = np.empty([len(t)+1])
       
xcs[0] = xs[-1]
ycs[0] = ys[-1] #positions stay the same
dxcs[0] = dxs[-1] + dv * np.cos(theta)
dycs[0] = dys[-1] + dv * np.sin(theta) #velocities need to experience the boost

print(f'Boost to start elliptical = {dv}')
def contRungeK():
    
    for i in range(len(t)):
        dE = np.sqrt((xcs[i] - x1[i])**2 + (ycs[i] - y1[i])**2)
        dM = np.sqrt((xcs[i] - x2[i])**2 + (ycs[i] - y2[i])**2)
        d2x = (-G * mE * (xcs[i] - x1[i]) / (dE**3)) - (G * mM * (xcs[i] - x2[i]) / (dM**3))
        d2y = (-G * mE * (ycs[i] - y1[i]) / (dE**3)) - (G * mM * (ycs[i] - y2[i]) / (dM**3))
        
        zx1 = xcs[i] + (dt * dxcs[i] / 2)
        dzx1 = dxcs[i] + (dt * d2x / 2)
        zy1 = ycs[i] + (dt * dycs[i] / 2)
        dzy1 = dycs[i] + (dt * d2y / 2)
        
        dE = np.sqrt((zx1 - x1[i])**2 + (zy1 - y1[i])**2)
        dM = np.sqrt((zx1 - x2[i])**2 + (zy1 - y2[i])**2)
        zx2 = xcs[i] + (dt * dzx1 / 2)
        dzx2 = dxcs[i] + ((-G * mE * (zx1 - x1[i]) / (dE**3)) - (G * mM * (zx1 - x2[i]) / (dM**3))) * (dt/2)
        zy2 = ycs[i] + (dt * dzy1 / 2)
        dzy2 = dycs[i] + ((-G * mE * (zy1 - y1[i]) / (dE**3)) - (G * mM * (zy1 - y2[i]) / (dM**3))) * (dt/2)
        
        dE = np.sqrt((zx2 - x1[i])**2 + (zy2 - y1[i])**2)
        dM = np.sqrt((zx2 - x2[i])**2 + (zy2 - y2[i])**2)
        zx3 = xcs[i] + (dt * dzx2)
        dzx3 = dxcs[i] + ((-G * mE * (zx2 - x1[i]) / (dE**3)) - (G * mM * (zx2 - x2[i]) / (dM**3))) * dt
        zy3 = ycs[i] + (dt * dzy2)
        dzy3 = dycs[i] + ((-G * mE * (zy2 - y1[i]) / (dE**3)) - (G * mM * (zy2 - y2[i]) / (dM**3))) * dt
        
        dE = np.sqrt((zx3 - x1[i])**2 + (zy3 - y1[i])**2)
        dM = np.sqrt((zx3 - x2[i])**2 + (zy3 - y2[i])**2)
        xcs[i+1] = xcs[i] + (dt/6) * (dxcs[i] + (2 * dzx1) + (2 * dzx2) + dzx3)
        dxcs[i+1] = dxcs[i] + (dt/6) * (d2x + (2 * (dzx2 - dxcs[i]) * (2/dt)) + (2 * (dzx3 - dxcs[i]) * (1/dt)) + ((-G * mE * (zx3 - x1[i]) / (dE**3)) - (G * mM * (zx3 - x2[i]) / (dM**3))))
        ycs[i+1] = ycs[i] + (dt/6) * (dycs[i] + (2 * dzy1) + (2 * dzy2) + dzy3)
        dycs[i+1] = dycs[i] + (dt/6) * (d2y + (2 * (dzy2 - dycs[i]) * (2/dt)) + (2 * (dzy3 - dycs[i]) * (1/dt)) + ((-G * mE * (zy3 - y1[i]) / (dE**3)) - (G * mM * (zy3 - y2[i]) / (dM**3))))
        v = np.sqrt(dxcs[i]**2 + dycs[i]**2)
        vs[i] = v
        #if v < 100 :
            #break
    #dist = np.sqrt((xs[-1] - x2[-1])**2 + (ys[-1] - y2[-1])**2)
    #edist = np.sqrt((xcs[-1] - x1[-1])**2 + (ycs[-1] - y1[-1])**2)
    
    #return dist #use this line for optidr

    #print(f'Initial Distance From Earth = {dr}')
    #print(f'Distance From Earth after analysis = {edist}')
    #print(f'Difference in distance = {abs(edist-dr)}')
    #print(dxs[-1],dxcs[0])
    #print(dys[-1],dycs[0])
    #print(alpha)
    #plt.plot(x1,y1, color = 'blue')
    #plt.plot(x2,y2, color = 'orange')
    #plt.plot(xcs,ycs, label = "Rocket Transfer", color = 'green')
    plt.legend()  
    
    
    

#contRungeK()


# create a time array up to t_max sampled at steps of dt
dt = 10
t_min = t_max
t_max = t_max + 4*T/2 #T/500
t = np.arange(t_min,t_max,dt)
nstep = int((t_max-t_min)/dt)


x1 = np.empty([len(t)])
y1 = np.empty([len(t)])
x2 = np.empty([len(t)])
y2 = np.empty([len(t)])

# Model the movement of the 2 bodies using simple harmonic motion
for i in range(len(t)):
    x1[i]=np.cos(omega*t[i])*rE  # set the massive particle positions
    y1[i]=np.sin(omega*t[i])*rE 
    x2[i]=-np.cos(omega*t[i])*rM
    y2[i]=-np.sin(omega*t[i])*rM


v = np.sqrt(dxcs[-1]**2 + dycs[-1]**2)
#print(f'Final velocity = {v}') 
    
vl=2 * np.pi * (rM-lr)/T

dv = vl - v
theta = np.arctan(dycs[-1]/dxcs[-1])
if dycs[-1] > 0 and dxcs[-1] < 0:
    theta = np.pi + theta
if dycs[-1] < 0 and dxcs[-1] < 0:
    theta = np.pi + theta
xccs = np.empty([len(t)+1])
yccs = np.empty([len(t)+1])
dxccs = np.empty([len(t)+1])
dyccs = np.empty([len(t)+1])
    #drcs = np.empty([len(t)+1])
       
xccs[0] = xcs[-1]
yccs[0] = ycs[-1] #positions stay the same
dxccs[0] = vl * np.cos(theta)#dxcs[-1] + dv * np.cos(theta)
dyccs[0] = vl * np.sin(theta)#dycs[-1] + dv * np.sin(theta) #velocities need to experience the boost

print(f'Boost to orbit the moon = {dv}')
print(f'Orbital Speed = {vl}')
def bcontRungeK():
    
    
    for i in range(len(t)):
        dE = np.sqrt((xccs[i] - x1[i])**2 + (yccs[i] - y1[i])**2)
        dM = np.sqrt((xccs[i] - x2[i])**2 + (yccs[i] - y2[i])**2)
        d2x = (-G * mE * (xccs[i] - x1[i]) / (dE**3)) - (G * mM * (xccs[i] - x2[i]) / (dM**3))
        d2y = (-G * mE * (yccs[i] - y1[i]) / (dE**3)) - (G * mM * (yccs[i] - y2[i]) / (dM**3))
        
        zx1 = xccs[i] + (dt * dxccs[i] / 2)
        dzx1 = dxccs[i] + (dt * d2x / 2)
        zy1 = yccs[i] + (dt * dyccs[i] / 2)
        dzy1 = dyccs[i] + (dt * d2y / 2)
        
        dE = np.sqrt((zx1 - x1[i])**2 + (zy1 - y1[i])**2)
        dM = np.sqrt((zx1 - x2[i])**2 + (zy1 - y2[i])**2)
        zx2 = xccs[i] + (dt * dzx1 / 2)
        dzx2 = dxccs[i] + ((-G * mE * (zx1 - x1[i]) / (dE**3)) - (G * mM * (zx1 - x2[i]) / (dM**3))) * (dt/2)
        zy2 = yccs[i] + (dt * dzy1 / 2)
        dzy2 = dyccs[i] + ((-G * mE * (zy1 - y1[i]) / (dE**3)) - (G * mM * (zy1 - y2[i]) / (dM**3))) * (dt/2)
        
        dE = np.sqrt((zx2 - x1[i])**2 + (zy2 - y1[i])**2)
        dM = np.sqrt((zx2 - x2[i])**2 + (zy2 - y2[i])**2)
        zx3 = xccs[i] + (dt * dzx2)
        dzx3 = dxccs[i] + ((-G * mE * (zx2 - x1[i]) / (dE**3)) - (G * mM * (zx2 - x2[i]) / (dM**3))) * dt
        zy3 = yccs[i] + (dt * dzy2)
        dzy3 = dyccs[i] + ((-G * mE * (zy2 - y1[i]) / (dE**3)) - (G * mM * (zy2 - y2[i]) / (dM**3))) * dt
        
        dE = np.sqrt((zx3 - x1[i])**2 + (zy3 - y1[i])**2)
        dM = np.sqrt((zx3 - x2[i])**2 + (zy3 - y2[i])**2)
        xccs[i+1] = xccs[i] + (dt/6) * (dxccs[i] + (2 * dzx1) + (2 * dzx2) + dzx3)
        dxccs[i+1] = dxccs[i] + (dt/6) * (d2x + (2 * (dzx2 - dxccs[i]) * (2/dt)) + (2 * (dzx3 - dxccs[i]) * (1/dt)) + ((-G * mE * (zx3 - x1[i]) / (dE**3)) - (G * mM * (zx3 - x2[i]) / (dM**3))))
        yccs[i+1] = yccs[i] + (dt/6) * (dyccs[i] + (2 * dzy1) + (2 * dzy2) + dzy3)
        dyccs[i+1] = dyccs[i] + (dt/6) * (d2y + (2 * (dzy2 - dyccs[i]) * (2/dt)) + (2 * (dzy3 - dyccs[i]) * (1/dt)) + ((-G * mE * (zy3 - y1[i]) / (dE**3)) - (G * mM * (zy3 - y2[i]) / (dM**3))))
    
    #dist = np.sqrt((xs[-1] - x2[-1])**2 + (ys[-1] - y2[-1])**2)
    #edist = np.sqrt((xcs[-1] - x1[-1])**2 + (ycs[-1] - y1[-1])**2)
    
    #return dist #use this line for optidr

    #print(f'Initial Distance From Earth = {dr}')
    #print(f'Distance From Earth after analysis = {edist}')
    #print(f'Difference in distance = {abs(edist-dr)}')
    #print(dxs[-1],dxcs[0])
    #print(dys[-1],dycs[0])
    #print(alpha)
    #plt.plot(x1,y1 , color = 'blue',)
    #plt.plot(x2,y2, color = 'orange')#, marker = 'o')
    #plt.plot(xccs,yccs, label = "Rocket Final", color = 'brown')#, marker = '>')
    plt.legend()


#bcontRungeK()

