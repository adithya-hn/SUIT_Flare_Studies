
''''
Created on 17 oct 2025
@author: adithya-hn

A simple example to load and plot STIX spectrum and lightcurve using SunPy and SunKit-Spex

'''


import matplotlib.pyplot as plt
from parfive import Downloader
import astropy.units as u
from astropy.time import Time
from sunkit_spex.extern.stix import STIXLoader
from sunkit_spex.legacy.fitting.fitter import Fitter, load
from datetime import datetime , timedelta

#---------------Input parameters----------------

spec_file="/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/data/stix/stx_spectrum_2411012243.fits"
srm_file="/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/data/stix/stx_srm_2411012243.fits"
Start_t="2024-11-01T12:30:30"
End_t="2024-11-01T14:45:30"
case='8_Nov01'

#-----------------------------------------------

time_profile_size = (9, 6)
spec_plot_size = (10, 10)
joint_spec_plot_size = (25, 10)
tol = 1e-20
spec_font_size = 18
xlims, ylims = [3, 100], [1e-1, 1e6]

plt.rcParams["font.size"] = spec_font_size
plt.figure(layout="tight",figsize=(12,8))
stix_spec = STIXLoader(spectrum_file=spec_file, srm_file=srm_file)
stix_spec.update_event_times(start=Time(Start_t), end=Time(End_t))
plot_range=[datetime.fromisoformat(Start_t)- timedelta(minutes=10),datetime.fromisoformat(End_t)+ timedelta(minutes=10)]
stix_spec.lightcurve(energy_ranges=[[4, 8],[9,12],[13,15], [16, 30]])
plt.xlim(plot_range[0],plot_range[1])
plt.xticks(rotation=45)
plt.savefig(f"stix_{case}_{Start_t}_{End_t}.png", dpi=300)
plt.show()