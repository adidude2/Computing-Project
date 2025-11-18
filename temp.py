import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

#Define universal gravitation constant
G=6.6726e-11 #N-m2/kg2

#Reference quantities
mE = 5.97219e24
mM = 7.35e22
d = 3.844e8
rE = d * (mM / (mM + mE))
rM = d * (mE / (mM + mE))
T = (((4*np.pi**2)/(G*(mE + mM))) * d**3)**0.5
#T = 2419200
omega=2*np.pi/T

# create a time array up to t_max sampled at steps of dt
dt = 10
t_min = 0.
t_max =   T*2
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
    
    
    
mew = (mM / (mM + mE))

roots = np.roots([1, mew-3, 3-(2*mew), -mew, 2*mew, -mew])
print(f'Roots = {roots[4]}')
print(type(roots[4]))
Y = roots[4]
A = Y.real
lr = d * A
print(lr)

# first postionas and speeds for the rocket
#multiplier for taylor ex = 1.04860411065608946
dr = d*(mM/(3 * mE))**(1/3) #+ 2991932.6395
x= x2[0] + lr
y= 0
r=np.sqrt(x*x+y*y)
w = (2 * np.pi)/ T
vx = 0
vy=-2 * np.pi * (rM-lr)/T # velocity needed to stay in the L2 point
#print(dr + 2991932.6395)
#print(dr + 3664097.11465)
print (x)

xs = np.empty([len(t)+1])
ys = np.empty([len(t)+1])
dxs = np.empty([len(t)+1])
dys = np.empty([len(t)+1])
drs = np.empty([len(t)+1])
xs[0] = x
ys[0] = y
dxs[0] = vx
dys[0] = vy
    
def TaylorEx():
    for i in range(len(t)):
        dE = np.sqrt((xs[i] - x1[i])**2 + (ys[i] - y1[i])**2)
        dM = np.sqrt((xs[i] - x2[i])**2 + (ys[i] - y2[i])**2)
        d2x = (-G * mE * (xs[i] - x1[i]) / (dE**3)) - (G * mM * (xs[i] - x2[i]) / (dM**3))
        d2y = (-G * mE * (ys[i] - y1[i]) / (dE**3)) - (G * mM * (ys[i] - y2[i]) / (dM**3))
        
        xs[i+1] = xs[i] + dt * dxs[i] + ((dt**2)/2) * d2x
        dxs[i+1] = dxs[i] + dt * d2x
        ys[i+1] = ys[i] + dt * dys[i] + ((dt**2)/2) * d2y
        dys[i+1] = dys[i] + dt * d2y
        
    dist = np.sqrt((xs[-1] - x2[-1])**2 + (ys[-1] - y2[-1])**2)

    print(f'Initial Distance From Moon = {dr}')
    print(f'Distance From Moon after one orbit = {dist}')
    print(f'Difference in distance = {abs(dist-dr)}')
    plt.plot(x1,y1, label = "Earth")
    plt.plot(x2,y2, label = "Moon")
    plt.plot(xs,ys, label = "Rocket")
    plt.legend()
    
def RungeK():
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
    
    dist = np.sqrt((xs[-1] - x2[-1])**2 + (ys[-1] - y2[-1])**2)
    
    #return dist #use this line for optidr
    print(d2x)    
    print(f'Initial Distance From Moon = {dr}')
    print(f'Distance From Moon after one orbit = {dist}')
    print(f'Difference in distance = {abs(dist-dr)}')
    plt.plot(x1,y1, label = "Earth")
    plt.plot(x2,y2, label = "Moon")
    plt.plot(xs,ys, label = "Rocket")
    plt.legend()
    

def optidr():
    a = 3664097.11465
    b = 3664098.114755
    c = (0.000001)
    r = np.arange(a, b, c)
    drs = np.empty([len(r)])
    dins = np.empty([len(r)])
    for i in range(len(r)):
        ds = dr + (a + i * c)
        dins[i] = i
        x= x2[0] - ds
        xs[0] = x
        dl = RungeK()
        drs[i] = abs(dl-ds)
        print(abs(dl-ds))
        
        
    #print(f'Initial Distance From Moon = {ds}')
    #print(f'Distance From Moon after one orbit = {dl}')
    #print(f'Difference in distance = {abs(dl-ds)}')
    #plt.plot(x1,y1, label = "Earth")
    #plt.plot(x2,y2, label = "Moon")
    #plt.plot(xs,ys, label = "Rocket")
    #plt.legend()
    plt.plot(dins, drs,)
        
        
        
RungeK()
#optidr()
#TaylorEx()

#IF YOU HAVE TIME TO FIX IT FIX THE OPTIDRS VY VELOCITY BCS IT'S AFFECTED BY THE PERIOD.





