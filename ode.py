import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

G = 6.67408*10**-11
M = 1.989*10**30
AU = 1.4959787*10**11
omega = 2*np.pi/(365.263629*24*60*60)
mu = G*M

t = np.linspace(0,365*24*60*60,365*40)
T = t/(24*60*60)

# ax = -mu/r^3 * x
# ay = -mu/r^3 * y
# r = np.sqrt(x**2 + y**2)

# vx = x'
# vx' = (-mu/r**3)*x

# vx(0) = 0
# vy(0) = r*omega

alpha0 = -19.55309963
alpha0R = alpha0*(np.pi/180)

def Solve_M(val):
    E = 2*val
    i = 1
    while i < 10:

        der = .5*(1 - e*np.cos(E))
        Fun = .5*(E - e*np.sin(E)) - val

        h = -1*Fun/der
        E = E + h

        i = i + 1
    return E    

def orbit(z,t,mu):
    x,vx,y,vy = z
    r = np.sqrt(x**2 + y**2)
    r3 = r**3
    dzdt = [vx,(-mu/r3)*x,vy,(-mu/r3)*y]
    return dzdt

z0 = [0.9832*AU,0,0,30290]

sol = odeint(orbit,z0,t,args=(mu,))

Xau = sol[:,0]/AU
Yau = sol[:,2]/AU

Vxkm = sol[:,1]/1000
Vykm = sol[:,3]/1000

thetaR = np.arctan(Yau/Xau)
thetaD = thetaR*(180/np.pi)

b = (max(Yau) - min(Yau))/2
a = (max(Xau)-min(Xau))/2
e = np.sqrt((a**2 - b**2)/a**2)

print(e)
print(a)
#plt.plot(Xau,Yau,'b')
plt.plot(T,Xau,'r')
plt.grid()
plt.show()