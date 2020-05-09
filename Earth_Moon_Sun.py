import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

G = 6.67408*10**-11
M = 1.989*10**30
m = 5.972*10**24
AU = 1.4959787*10**11
omega = 2*np.pi/(365.263629*24*60*60)
mu = G*M
a = .999771*AU
e = 0.0167

factor = np.sqrt((1-e)/(1+e))
SiderialDay = 0.99726956    # Solar Days (24 hr)
SiderialYear = 365.2563629  # Solar Days (24 hr)

AreaTotal = np.pi 

Area_Time_Ratio = AreaTotal/SiderialYear
Time_Area_Ratio = SiderialYear/AreaTotal

t = np.linspace(0,365*24*60*60,365*40)
T = t/(24*60*60)

thetaD = []
thetaR = []
R = []

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

for i in t:
    AreaSwept = Area_Time_Ratio*(i/(24*60*60))
    while AreaSwept > np.pi:
        AreaSwept = AreaSwept - np.pi

    E = Solve_M(AreaSwept)
    theta = 2*np.arctan(np.tan(E/2)/factor)
    thetad = theta*(180/np.pi)
    Re = a*(1-e**2)/(1 + e*np.cos(theta))

    thetaD.append(thetad)
    thetaR.append(theta)
    R.append(Re)

# ax = -mu/r^3 * x
# ay = -mu/r^3 * y
# r = np.sqrt(x**2 + y**2)

# vx = x'
# vx' = (-mu/r**3)*x

# vx(0) = 0
# vy(0) = r*omega

def orbit(z,t,G,M,m,thetaR):
    x,y,vx,vy = z
    R = a*(1-e**2)/(1+e*np.cos(thetaR))
    X = x - (R*np.cos(thetaR)+a*e)
    Y = y - (R*np.sin(thetaR))
    alpha = np.arctan(Y/X)
    r = X/np.cos(alpha)
    

    ax = -G*M*np.cos(thetaR)/R**2 - G*m*np.cos(alpha)/r**2
    ay = -G*M*np.sin(thetaR)/R**2 - G*m*np.sin(alpha)/r**2

    dzdt = [vx,vy,ax,ay]
    return dzdt


VxE0 = 0
VyE0 = 30290

VxM0 = 0
VyM0 = 1000

Vx = VxE0 + VxM0
Vy = VyE0 + VyM0

z0 = [(0.9832+.002562)*AU,-0.000834,Vx,Vy]

sol = odeint(orbit,z0,t,args=(G,M,m,thetaR))

Xau = sol[:,0]/AU
Yau = sol[:,1]/AU

Vxkm = sol[:,2]/1000
Vykm = sol[:,3]/1000





#plt.plot(Xau,Yau,'b')
plt.plot(Xau,Yau,'r')
plt.grid()
plt.show()