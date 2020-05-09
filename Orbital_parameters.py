import numpy as np 
import scipy as sp 

rtest = [-3670,-3870,4400]
vtest = [4.7,-7.4,1]

mu = 398600
J2earth = 0.00108263
Rearth = 6371

def mod(r):
    s = 0
    for i in r:
        s = s + i**2
    return s**.5
def cap(r):
    v = []
    for i in r:
        directionCosine = i/mod(r)
        v.append(directionCosine)
    return v
def Orbit_parameters_r_v(r,v):
    h = np.cross(r,v)
    H = mod(h)
    kcap = [0,0,1]

    rcap = cap(r)
    vcap = cap(v)

    vm = mod(v)
    rm = mod(r)

    vradial = np.dot(v,r)/rm

    N = np.cross(kcap,h)
    Ncap = cap(N)
    
    c1 = (1/mu)*(vm**2 - mu/rm)
    c2 = -(1/mu)*np.dot(r,v)
    
    a1 = np.multiply(c1,r)  
    a2 = np.multiply(c2,v)  
    
    e = np.add(a1,a2)   # eccentricity vector
    ecap = cap(e)


    em = mod(e) # eccentricity
    i = (np.arccos(h[2]/mod(h)))*(180/np.pi)    # inclination in degrees
    ascensionNode = (np.arccos(N[0]/mod(N)))*(180/np.pi)
    argumentPerigee = np.arccos(np.dot(Ncap,ecap))*(180/np.pi)
    trueAnomoly = np.arccos(np.dot(rcap,ecap))*(180/np.pi)

    rp = (H**2/mu)*(1/(1+em))
    ra = (H**2/mu)*(1/(1-em))

    a = .5*(ra + rp)
    T = (2*np.pi)*(a**1.5)/np.sqrt(mu)

    if N[1] < 0 :
        if ascensionNode < 180 :
            ascensionNode = 360 - ascensionNode
    if N[1] > 0:
        if ascensionNode > 180:
            ascensionNode = 360 -ascensionNode
    
    if vradial > 0: # flying away from peregee theta within 180
        if trueAnomoly > 180 :
            trueAnomoly = 360-trueAnomoly

    if vradial < 0: # flying tow peregee theta greater than 180
        if trueAnomoly < 180 :
            trueAnomoly = 360-trueAnomoly

    if e[2] > 0 :
        if argumentPerigee > 180 :
            argumentPerigee = 360 - argumentPerigee
    if e[2] < 0:
        if argumentPerigee < 180 :
            argumentPerigee = 360 - argumentPerigee
    
    
    JtN =   -(1.5*np.sqrt(mu)*J2earth*(Rearth**2))
    JtD =   (a**3.5)*(1-em**2)**2
    
    Jterm = JtN/JtD
    
    ut = Jterm*np.cos(i*(np.pi)/180)       # Rad per sec
    wt = Jterm*(2.5*(np.sin(i*(np.pi/180)))**2 - 2)        # Rad per sec
    
    dot = [ut,wt]
    O = [em,i,ascensionNode,argumentPerigee,trueAnomoly,a]
    return O,dot,H,T,rm
def printPar(O):
    print("eccentricity of orbit e  ",O[0])
    print("inclination of orbit i  ",O[1])
    print("AscNode of orbit U ",O[2])
    print("arg of perigee of orbit w  ",O[3])
    print("true anomoly of orbit theta  ",O[4])
    print("semi major axis of orbit a  ",O[5])
def solve_M_E(M0,e):
    i = 1
    E = M0
    while i < 10:
        Fun = E - e*np.sin(E) - M0
        der = 1 - e*np.cos(E)
        h = -Fun/der
        E = E + h
        i = i + 1
    return E
def Matrix_periFocal_geocentric(i,u,w):
    I = [[1,0,0],[0,np.cos(i),np.sin(i)],[0,-np.sin(i),np.cos(i)]]
    W = [[np.cos(w),np.sin(w),0],[-np.sin(w),np.cos(w),0],[0,0,1]]
    U = [[np.cos(u),np.sin(u),0],[-np.sin(u),np.cos(u),0],[0,0,1]]

    Q1 = np.matmul(W,I)
    Q = np.matmul(Q1,U)
    Qt = Q.T
    return Qt
def perifocal_coordinates(r,v,deltaT):
    orbit = Orbit_parameters_r_v(r,v)[0]

    e = orbit[0]
    i = orbit[1]*(np.pi/180)
    u = orbit[2]*(np.pi/180)
    w = orbit[3]*(np.pi/180)
    theta = orbit[4]*(np.pi/180)

    a = orbit[5]
    h = Orbit_parameters_r_v(r,v)[2]
    T = Orbit_parameters_r_v(r,v)[3]
    rr = Orbit_parameters_r_v(r,v)[4]
    n = 2*np.pi/T   # Mean Rotation speed
    factor = np.sqrt(1-e)/np.sqrt(1+e)

    wt = Orbit_parameters_r_v(r,v)[1][1]
    ut = Orbit_parameters_r_v(r,v)[1][0]

    E0 = 2*np.arctan(np.tan(theta/2)*factor)
    t0 = (E0 - e*np.sin(E0))/n

    tf = t0 + deltaT

    nf = tf/T
    deltaT_perigee  = (nf - int(nf))*T
    Mean_Anomaly_tf = n*deltaT_perigee
    Ef = solve_M_E(Mean_Anomaly_tf,e)
    thetaf = 2*np.arctan(np.tan(Ef/2)/factor)

    rtf = ((h**2)/mu)/(1 + e*np.cos(thetaf))

    r_pq = [[rtf*np.cos(thetaf)],[rtf*np.sin(thetaf)],[0]]
    v_pq = [[-(mu/h)*np.sin(thetaf)],[(mu/h)*(e+np.cos(thetaf))],[0]]


    wf = w + deltaT*wt
    uf = u + deltaT*ut

    rgeo_tf = np.matmul(Matrix_periFocal_geocentric(i,uf,wf),r_pq)
    vgeo_tf = np.matmul(Matrix_periFocal_geocentric(i,uf,wf),v_pq)

    return rgeo_tf,vgeo_tf
def orbitTransfer_Hohmann_P(rp1,ra1,rp2,ra2):   # firing at peregee
    a1 = (rp1+ra1)/2
    a2 = (rp2+ra2)/2
    vp1 = np.sqrt(mu*(2/rp1 -1/a1))
    vp2 = np.sqrt(mu*(2/rp1 -1/a2))
    dv = vp2 - vp1
    return dv
def orbitTransfer_Hohmann_A(rp1,ra1,rp2,ra2):   # firing at peregee
    a1 = (rp1+ra1)/2
    a2 = (rp2+ra2)/2
    va1 = np.sqrt(mu*(2/ra1 -1/a1))
    va2 = np.sqrt(mu*(2/ra1 -1/a2))
    dv = va2 - va1
    return dv
def orbitTransfer_Hohmann_Circularization_A(rp1,ra1):   # firing at peregee
    a1 = (rp1+ra1)/2
    r2 = ra1
    va1 = np.sqrt(mu*(2/ra1 -1/a1))
    vr2 = np.sqrt(mu*(1/r2))
    dv = vr2 - va1
    return dv
def orbitTransfer_Hohmann_Circularization_P(rp1,ra1):   # firing at peregee
    a1 = (rp1+ra1)/2
    r2 = rp1
    vp1 = np.sqrt(mu*(2/rp1 -1/a1))
    vr2 = np.sqrt(mu*(1/r2))
    dv = vr2 - vp1
    return dv

print(orbitTransfer_Hohmann_P(6858,7178,6858,22378))
