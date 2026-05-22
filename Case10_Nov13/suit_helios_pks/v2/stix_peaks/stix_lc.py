
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
from datetime import datetime , timedelta
import matplotlib.dates as mdates
from sys import path as sys_path
import pandas as pd
import scienceplots
plt.style.use('science')

#---------------Input parameters----------------

Start_t = "2024-11-13T14:57:00"
End_t   = "2024-11-13T17:12:00"


spec_file="stx_spectrum_2411132589.fits"
srm_file="stx_srm_2411132589.fits" 
case='10_Nov13'


#-----------------------------------------------

time_profile_size = (9, 6)
spec_plot_size = (10, 10)
joint_spec_plot_size = (25, 10)
tol = 1e-20
spec_font_size = 18
# xlims, ylims = [3, 100], [1, 1e6]

plt.rcParams["font.size"] = spec_font_size
plt.figure(layout="tight",figsize=(12,8))
stix_spec = STIXLoader(spectrum_file=spec_file, srm_file=srm_file)
stix_spec.update_event_times(start=Time("2024-11-13T12:10:00"), end=Time("2024-11-13T12:11:01"))
# stix_spec.update_background_times(start=Time(start_background_time), end=Time(end_background_time))

plot_range=[datetime.fromisoformat(Start_t)- timedelta(minutes=0),datetime.fromisoformat(End_t)+ timedelta(minutes=1)]
lc=stix_spec.lightcurve(energy_ranges=[[4, 9], [10, 30]])
plt.axvline(datetime.fromisoformat("2024-11-01T02:13:25"),ls='--',lw=2,color='b',label='Impuslive phase start time',alpha=0.7)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=5))
plt.xlim(plot_range[0],plot_range[1])
plt.title ('STIX Lightcurve (2024-Nov-01)')
plt.xlabel('Time (UT)')
# plt.xticks(rotation=45)
plt.legend(loc='upper left')
plt.ylim(1,1e5)
plt.savefig(f"stix_{case}_onest_lc.png", dpi=300)
plt.close()


energy_ranges = [[4, 9],[10, 30]]
band_labels   = ["4-9 keV", "10-30 keV"]



# lines = ax.get_lines()   # list of all curves
all_lines  = lc.get_lines()
data_lines = all_lines[1::2]   # indices 1, 3, 5, 7  ← real cadence (4879 pts)

lc_data = {}
for line, label in zip(data_lines, band_labels):
    times  = mdates.num2date(line.get_xdata())
    counts = line.get_ydata()
    lc_data[label] = {"times": times, "counts": counts}
    print(f"{label}: {len(times)} pts | [{counts.min():.2f}, {counts.max():.2f}]")

# ── DataFrame ─────────────────────────────────────────────────────────────────
df = pd.DataFrame({"time": lc_data[band_labels[0]]["times"]})
for label in band_labels:
    df[label] = lc_data[label]["counts"]

df.to_csv("stix_lightcurves.csv", index=False)