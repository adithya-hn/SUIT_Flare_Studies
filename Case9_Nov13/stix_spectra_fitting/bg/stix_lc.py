
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
# spec_file="../../data/stix/stx_spectrum_2411127496.fits"
# srm_file="../../data/stix/stx_srm_2411127496.fits"
# Start_t="2024-11-12T22:22:00" #no date before this time
# End_t="2024-11-12T23:12:00"

spec_file="../../data/stix/stx_spectrum_2411125825.fits"
srm_file="../../data/stix/stx_srm_2411125825.fits"
Start_t="2024-11-12T23:08:00" 
End_t="2024-11-13T00:20:30"
case='9_Nov13'

start_background_time = "2024-11-12T23:08:00"
end_background_time   = "2024-11-12T23:11:00"
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
# print(stix_spec.header)
stix_spec.update_event_times(start=Time(Start_t), end=Time(End_t))
stix_spec.update_background_times(start=Time(start_background_time), end=Time(end_background_time))

plot_range=[datetime.fromisoformat(Start_t)- timedelta(minutes=10),datetime.fromisoformat(End_t)+ timedelta(minutes=10)]
stix_spec.lightcurve(energy_ranges=[[4, 8],[9,12],[13,25], [22, 30]])
plt.xlim(plot_range[0],plot_range[1])
plt.xticks(rotation=45)
plt.savefig(f"stix_{case}_{Start_t}_{End_t}.png", dpi=300)
plt.show()