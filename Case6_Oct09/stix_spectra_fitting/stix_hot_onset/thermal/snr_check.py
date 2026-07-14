from astropy.time import Time, TimeDelta
import pandas as pd
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from parfive import Downloader
import astropy.units as u
from astropy.time import Time
from sunkit_spex.extern.stix import STIXLoader
from sunkit_spex.legacy.fitting.fitter import Fitter, load
from datetime import datetime, timedelta
import astropy.time as atime
import numpy as np


flare_start = Time("2024-10-09T01:24:40")
flare_end   = Time("2024-10-09T01:25:00")
#---------------Input parameters----------------
spec_file="../stx_spectrum_2410088145.fits"
srm_file="../stx_srm_2411012243.fits"

case='6_Oct9'

start_background_time = "2024-10-08T23:33:30"
end_background_time   = "2024-10-08T23:37:00"

#-----------------------------------------------



bin_size_20 = TimeDelta(20, format='sec')
bin_size_30 = TimeDelta(30, format='sec')

outfile = f"{case}_stix_timeResolved_nth_fit.csv"

columns = ["time_start","time_end",
"T","T_er1","T_er2",
"EM","EM_er1","EM_er2","L"]

df = pd.DataFrame(columns=columns)
log = {"20s": [], "30s": [], "skipped": []}

t = flare_start
i = 0
stix_spec = STIXLoader(spectrum_file=spec_file, srm_file=srm_file)
t1 = t + bin_size_20
stix_spec.update_event_times(start=t, end=t1)
stix_spec.update_background_times(atime.Time(start_background_time), atime.Time(end_background_time))
# print([a for a in dir(stix_spec) if any(k in a.lower() for k in ('count', 'bg', 'background', 'error'))])
counts_total = stix_spec._count_rate_perspec         # or stix_spec.data.counts, etc.
counts_total_err = stix_spec._count_rate_error_perspec   # if available; else sqrt(N)
# print(f" Total counts: {len(counts_total)} ± {counts_total_err}")
# print(len(counts_total))
spec_data = stix_spec._loaded_spec_data

# print(f" Counts: {len(counts)} ± {np.sqrt(counts)}")
# print(spec_data.keys())
bg_counts = spec_data["extras"]["background_counts"]           # per-channel background counts, over the background window
bg_count_error = spec_data["extras"]["background_count_error"]
bg_rate = spec_data["extras"]["background_rate"]                # counts/s per channel
bg_rate_error = spec_data["extras"]["background_rate_error"]
bg_eff_exposure = spec_data["extras"]["background_effective_exposure"]

count_rate = spec_data["count_rate"]  # counts/s per channel
count_rate_er = spec_data["count_rate_error"]  # counts er/s per channel
energy_bins = spec_data["count_channel_bins"]   # energy edges for each channel
energy_mids = spec_data["count_channel_mids"]

fig, ax = plt.subplots()

# step/staircase plot, respecting bin widths (standard for spectra)
edges = np.append(energy_bins[:, 0], energy_bins[-1, 1])
# print(f"Edges: {edges}")
ax.stairs(count_rate, edges, label="Event count rate")

# add error bars at bin centers
ax.errorbar(energy_mids, count_rate, yerr=count_rate_er, fmt="-", ecolor="black", alpha=0.6)
# ax.errorbar(energy_mids, bg_rate, yerr=bg_rate_error, fmt="none", ecolor="orange", alpha=0.6)
# plot background rate as a dashed line
ax.step(energy_mids, bg_rate, where="mid", label="Background rate", color="gray", linestyle="--")

ax.set_xlabel("Energy [keV]")
ax.set_ylabel("Count rate [counts s$^{-1}$ keV$^{-1}$]")
# ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()
ax.set_xlim(6,30)
plt.close()

# print(f" Channel bins: {energy_bins}")

# net (background-subtracted) signal
net_rate = count_rate - bg_rate

# Poisson errors add in quadrature — both source and bkg contribute noise
net_rate_error = np.sqrt(count_rate_er**2 + bg_rate_error**2)

# avoid divide-by-zero for empty bins
with np.errstate(divide="ignore", invalid="ignore"):
    snr = np.where(net_rate_error > 0, net_rate / net_rate_error, 0)

energy_mids = spec_data["count_channel_mids"]

fig, ax = plt.subplots()
ax.plot(energy_mids, snr, drawstyle="steps-mid")
ax.axhline(3, color="red", linestyle="--", label="SNR = 3 threshold")
ax.set_xlabel("Energy [keV]")
ax.set_ylabel("SNR")
# ax.set_xscale("log")
ax.set_xlim(6,30)
plt.legend()
plt.show()