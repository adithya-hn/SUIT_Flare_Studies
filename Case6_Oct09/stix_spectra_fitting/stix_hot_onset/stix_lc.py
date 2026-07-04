
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

spec_file="/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/stix/with_bg/stx_spectrum_2410088145.fits"
srm_file="/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/data/stix/stx_srm_2411012243.fits"
Start_t="2024-10-09T01:00:00"
End_t="2024-10-09T01:50:00"
case='6_Oct9'


start_background_time = "2024-10-09T01:14:30"
end_background_time   = "2024-10-09T01:19:00"

thermal_range    = [datetime.fromisoformat("2024-10-09T01:20:00"),
                    datetime.fromisoformat("2024-10-09T01:30:40")]

nonthermal_range = [datetime.fromisoformat("2024-10-09T01:30:40"),
                    datetime.fromisoformat("2024-10-09T01:33:00")]

hot_onset_time   =  datetime.fromisoformat("2024-10-09T01:24:40")



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
stix_spec.update_event_times(start=Time("2024-10-08T23:29:00"), end=Time("2024-10-08T23:29:30"))
stix_spec.update_background_times(start=Time(start_background_time), end=Time(end_background_time))
plt.axvspan(*thermal_range,    alpha=0.1, color="tomato",      zorder=0)#, label="Thermal fit range")
plt.axvspan(*nonthermal_range, alpha=0.1, color="mediumpurple", zorder=0)#, label="Non-thermal fit range")
plot_range=[datetime.fromisoformat(Start_t)- timedelta(minutes=0),datetime.fromisoformat(End_t)+ timedelta(minutes=1)]
lc=stix_spec.lightcurve(energy_ranges=[[4, 8],[9,12],[13,25], [22, 30]])
plt.axvline(datetime.fromisoformat("2024-10-09T01:31:10"),ls='-',lw=1,color='b',label='Impuslive phase start time')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=5))
plt.xlim(plot_range[0],plot_range[1])
plt.title ('STIX Lightcurve (2024-Oct-09)')
plt.xlabel('Time (UT)')
plt.xticks(rotation=45)
plt.legend(loc='upper right')
plt.ylim(5,1e6)
plt.savefig(f"stix_{case}_onest_lc.png", dpi=300)
plt.show()


energy_ranges = [[4, 8], [9, 12], [13, 22], [22, 30]]
band_labels   = ["4-8 keV", "9-12 keV", "13-22 keV", "22-30 keV"]



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