import matplotlib.pyplot as plt
from parfive import Downloader

import astropy.units as u
from astropy.time import Time

from sunkit_spex.extern.stix import STIXLoader
from sunkit_spex.legacy.fitting.fitter import Fitter, load


#---------------Input parameters----------------

spec_file="../data/stix/stx_spectrum_2411138743.fits"
srm_file="../data/stix/stx_srm_2411138743.fits"


#-----------------------------------------------

time_profile_size = (9, 6)
spec_plot_size = (10, 10)
joint_spec_plot_size = (25, 10)
tol = 1e-20
spec_font_size = 18
xlims, ylims = [3, 100], [1e-1, 1e6]

plt.rcParams["font.size"] = spec_font_size


stix_spec = STIXLoader(spectrum_file=spec_file, srm_file=srm_file)
plt.figure(layout="tight")

# the line that actually plots
stix_spec.lightcurve(energy_ranges=[[4, 15], [15, 30], [30, 60]])
plt.xticks(rotation=45)
#plt.xlabel("Time (UT)")
plt.savefig("stix_lightcurve.png", dpi=300)
plt.show()