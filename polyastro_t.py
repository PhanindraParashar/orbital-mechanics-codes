import numpy as np
import matplotlib.pyplot as plt
from poliastro.plotting import plot
from astropy import units as u
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from poliastro.examples import iss
from poliastro import iod
plt.ion()  # To immediately show plots

plt.style.use("seaborn")  # Recommended

r = [-12045, -3490, 2500] * u.km
v = [-3.457, 6.618, 2.533] * u.km / u.s

ss = Orbit.from_vectors(Earth, r, v)


plot(ss)

date_launch = time.Time('2011-11-26 15:02', scale='utc')
date_arrival = time.Time('2012-08-06 05:17', scale='utc')
tof = date_arrival - date_launch

ss0 = Orbit.from_body_ephem(Earth, date_launch)
ssf = Orbit.from_body_ephem(Mars, date_arrival)


(v0, v), = iod.lambert(Sun.k, ss0.r, ssf.r, tof)

iss.propogate(30*u.min)