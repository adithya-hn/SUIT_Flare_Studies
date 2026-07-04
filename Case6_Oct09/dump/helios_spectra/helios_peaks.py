
'''
CREATE ON 7 DEC 2025 
AUTHOR: ADITHYA 
PURPOSE: FIND THE PEAKS IN THE HELIOS LIGHT CURVE DATA
'''

import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import find_peaks
from scipy.stats import poisson
import matplotlib.pyplot as plt

data=np.loadtxt('binned_counts_12_30kev_lc.csv',delimiter=',',skiprows=1,dtype='str')
dt=np.array(data[:,0],dtype='datetime64')
counts=np.array(data[:,1],dtype=float)
# c is 1D array of counts
window = 10  # choose based on timescale over which you expect NO flares
bkg = median_filter(counts, size=window)
# Avoid zeros in background to prevent division by zero
bkg_safe = np.clip(bkg, 1e-6, None)
sigma = np.sqrt(bkg_safe)

z = (counts - bkg_safe) / sigma    # "Gaussian-equivalent" significance per bin

sigma_thresh = 3.0  # or 4, 5, etc.
min_separation_bins = 3  # minimum distance between distinct peaks

peaks, props = find_peaks(
    z,
    height=sigma_thresh,    # only bins with z >= sigma_thresh
    distance=min_separation_bins
)

peak_times = dt[peaks]
peak_counts = counts[peaks]
peak_bkg   = bkg_safe[peaks]
peak_sig   = z[peaks]  # or z_poisson[peaks]


plt.errorbar(dt,bkg,yerr=np.sqrt(bkg_safe))
for pk in peak_times:
    plt.axvline(pk,alpha=0.4)
plt.errorbar(dt,counts,yerr=np.sqrt(counts))
plt.yscale('log')
plt.show()