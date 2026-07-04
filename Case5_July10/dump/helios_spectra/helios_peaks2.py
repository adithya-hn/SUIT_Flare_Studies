
'''
CREATE ON 7 DEC 2025 
AUTHOR: ADITHYA 
PURPOSE: FIND THE PEAKS IN THE HELIOS LIGHT CURVE DATA
CORRECTION HISTORY:
02 JAN 2026 : MEDIAN OF 15 WITH ERROR IS APPLIED.


'''

import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import find_peaks
from scipy.stats import poisson
import matplotlib.pyplot as plt
import datetime

data=np.loadtxt('helios_CdTe_c5.csv',delimiter=',',skiprows=1,dtype='str')
dt=np.array(data[:,0],dtype='datetime64')
counts=np.array(data[:,1],dtype=float)
count_er=np.array(data[:,2],dtype=float)
# c is 1D array of counts
window = 15  # choose based on timescale over which you expect NO flares
bkg = median_filter(counts, size=window)

# window = 15
# n_mc = 5000

# mc_bkg = np.zeros((n_mc, len(counts)))

# for i in range(n_mc):
#     x = np.random.poisson(counts)
#     mc_bkg[i] = median_filter(x, size=window)

# bkg = np.median(mc_bkg, axis=0)
# bkg_err = np.std(mc_bkg, axis=0)

# Avoid zeros in background to prevent division by zero
bkg_safe = np.clip(bkg, 1e-6, None)
sigma = np.sqrt(bkg_safe)
z = (counts - bkg_safe) / (sigma/np.sqrt(window))    # "Gaussian-equivalent" significance per bin
sigma_thresh = 2  # or 4, 5, etc.
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
plt.errorbar(dt,bkg,yerr=np.sqrt(bkg_safe)/np.sqrt(window))
for pk in peak_times:
    plt.axvline(pk,alpha=0.4)
plt.errorbar(dt,counts,yerr=count_er)
plt.yscale('log')
plt.savefig('helios_peak.png',dpi=300)
plt.show()