import numpy as np
import matplotlib.pyplot as plt
plt.ion()  # To immediately show plots

from astropy import units as u

from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit

plt.style.use("seaborn")