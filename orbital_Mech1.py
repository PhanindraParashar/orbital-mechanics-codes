from orbital import earth, KeplerianElements, Maneuver, plot

from scipy.constants import kilo
import matplotlib.pyplot as plt

orbit = KeplerianElements.with_altitude(1000 * kilo, body=earth)
man = Maneuver.hohmann_transfer_to_altitude(10000 * kilo)
plot(orbit, title='Maneuver 1', maneuver=man)
plt.show()