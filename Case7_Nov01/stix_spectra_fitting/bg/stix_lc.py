
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

spec_file="../../data/stix/stx_spectrum_2410315184.fits"
srm_file="../../data/stix/stx_srm_2410315184.fits"
Start_t="2024-10-31T09:16:30"
End_t="2024-11-01T02:20:00"
case='7_Nov01'


start_background_time = "2024-11-01T01:07:30"
end_background_time   = "2024-11-01T01:10:00"
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
stix_spec.update_background_times(start=Time(start_background_time), end=Time(end_background_time))

plot_range=[datetime.fromisoformat(Start_t)- timedelta(minutes=10),datetime.fromisoformat(End_t)+ timedelta(minutes=10)]
stix_spec.lightcurve(energy_ranges=[[4, 8],[9,12],[13,25], [22, 30]])
plt.xlim(plot_range[0],plot_range[1])
plt.xticks(rotation=45)
plt.savefig(f"stix_{case}_{Start_t}_{End_t}.png", dpi=300)
plt.show()