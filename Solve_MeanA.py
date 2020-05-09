import numpy as np
import scipy as sp
from sympy import *

e = 0.0167 
val = 213*np.pi/180
E = val
i = 1
while i < 10:
    area = .5*(E - e*np.sin(E))
    
    der = .5*(1 - e*np.cos(E))
    Fun = area - val

    h = -1*Fun/der
    
    print(h)
    E = E + h

    i = i + 1

print(E*180/(np.pi) - 360)