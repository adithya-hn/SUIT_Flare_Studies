import subprocess
import os
from astropy.io import fits
import glob
from datetime import datetime,timedelta
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter
from scipy.signal import find_peaks

command = [
    "python3",
    "LcFromEvtFile_V06.py",
    "--fits_file", "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/helios/2024/11/01/HLS_20241031_235953_43193sec_lev1_V111/events/evt.fits",
    "--mode", "CDTE_TOTAL",
    "--start_time", "2024-11-01T02:00:00.000",
    "--end_time",   "2024-11-01T02:31:00.000",
    "--time_bin_size", "30",
    "--output", "LightCurve.csv"
]

# Inputs that the script expects after running
user_inputs = "1\n22 30\nn \n"

result = subprocess.run(
    command,
    input=user_inputs,
    text=True,
    capture_output=True
)

print(result.stdout)
print(result.stderr)
h_lc=pd.read_csv('LightCurve.csv')
ht=np.array(h_lc["Time (UTC)"], dtype='datetime64[s]');
h_count=np.array(h_lc["10.00-30.00 keV"],dtype=float)/60
h_er=np.sqrt(h_count)/60

t_cut = ht
h_count_cut = h_count
h_er_cut = h_er
plt.figure()
plt.errorbar(t_cut,h_count_cut,yerr=h_er_cut)
# plt.errorbar(lc_t,lc_count/60,yerr=np.sqrt(lc_count)/60)
plt.xlabel('Time-UTC',fontsize=14);plt.ylabel('Counts/s ',fontsize=14);
plt.tick_params(labelsize=14);
plt.yscale('log')
plt.show()

counts=h_count_cut
count_er=h_er_cut
dt=t_cut
window = 10  # choose based on timescale over which you expect NO flares
bkg = median_filter(counts, size=window)

# Avoid zeros in background to prevent division by zero
bkg_safe = np.clip(bkg, 1e-6, None)
sigma = 1.253*np.sqrt(bkg_safe*60)/60
z = (counts - bkg_safe) / (sigma/np.sqrt(window))    # "Gaussian-equivalent" significance per bin
sigma_thresh = 3  # or 4, 5, etc.
min_separation_bins = 2  # minimum distance between distinct peaks

peaks, props = find_peaks(
    z,
    height=sigma_thresh,    # only bins with z >= sigma_thresh
    distance=min_separation_bins
)

peak_times_ = (dt[peaks].astype('datetime64[s]')).astype(str)
peak_times = dt[peaks]
peak_counts = counts[peaks]
peak_bkg   = bkg_safe[peaks]
peak_sig   = z[peaks]  # or z_poisson[peaks]
#print((peak_times).astype('datetime64[s]'))
np.savetxt('helios_peaks.csv',np.c_[peak_times_,peak_counts],header='date_time,helio_count',comments='',delimiter=',',fmt='%s')



plt.figure(figsize=(14,8))
plt.errorbar(dt,bkg,yerr=np.sqrt(bkg_safe*60)/(60*np.sqrt(window)))
for pk in peak_times:
    plt.axvline(pk,alpha=0.4)
plt.errorbar(dt,counts,yerr=count_er)
plt.yscale('log')
plt.savefig('helios_peak.png',dpi=300)
plt.show()
