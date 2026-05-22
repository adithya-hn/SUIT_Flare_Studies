

'''
Created on 7 Feb 2026
@author: adithya-hn
purpose: find strong peaks in th suit peaks

'''



from astropy.convolution import convolve, Gaussian1DKernel
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.dates as mdates
import pandas as pd
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
import numpy.ma as ma
import seaborn as sns
import glob
import os
from scipy.optimize import curve_fit
from scipy.ndimage import median_filter
from scipy.signal import savgol_filter
import scienceplots
plt.style.use('science')

scol =sns.color_palette("colorblind")

suit_fl="../Diff_img_data_NB04.csv"
stix_df=pd.read_csv("stix_lightcurves.csv")
h_lc=pd.read_csv('Helios_LightCurve.csv')

ht=np.array(h_lc["Time (UTC)"], dtype='datetime64[s]');
h_count=np.array(h_lc["10.00-30.00 keV"],dtype=float)/60
h_er=np.sqrt(h_count)/60

# stix_df=np.array(stix_df["time"], dtype='datetime64[s]');
stix_df["time"] = pd.to_datetime(stix_df["time"])
stix_df = stix_df.set_index("time")
stix_binned = stix_df.resample("60s").sum().reset_index()

# print(stix_binned)
stix_t=np.array(stix_binned["time"], dtype='datetime64[s]');
stix_count=np.array(stix_binned["10-30 keV"],dtype=float)/60
stix_er=np.sqrt(stix_count)/60

all_times = []
all_ints = []

data = np.genfromtxt(suit_fl, delimiter=',', names=True, dtype=None)
flux=np.array(data['diff_count'])   
t = np.array(data['Date'], dtype='datetime64[ms]')
area=np.array(data['area'],dtype=float)

x = np.array(data['diff_count'], dtype=float)#/area

def find_peaks_min_rise(flux, min_rise_points=3):
    peaks = []
    n = len(flux)

    for i in range(min_rise_points, n-1):
        rising = True
        
        # Check consecutive rise
        for k in range(min_rise_points):
            if flux[i-k-1] >= flux[i-k]:
                rising = False
                break
        
        # Check slope change to negative
        if rising and flux[i] > flux[i+1]:
            peaks.append(i)

    return np.array(peaks)


kernel = Gaussian1DKernel(1)
flux_smooth = convolve(flux, kernel, boundary='extend')
peaks,p = find_peaks(flux)

window = 10  # choose based on timescale over which you expect NO flares
bkg = median_filter(h_count, size=window)
# Avoid zeros in background to prevent division by zero
bkg_safe = np.clip(bkg, 1e-6, None)
sigma = 1.253*np.sqrt(bkg_safe*60)/60
z = (h_count - bkg_safe) / (sigma/np.sqrt(window))    # "Gaussian-equivalent" significance per bin
sigma_thresh = 3  # or 4, 5, etc.
min_separation_bins = 2  # minimum distance between distinct peaks

h_peaks, props = find_peaks( z,
    height=sigma_thresh,    # only bins with z >= sigma_thresh
    distance=min_separation_bins)

peak_times_ = (ht[h_peaks].astype('datetime64[s]')).astype(str)
peak_counts = h_count[h_peaks]
peak_bkg   = bkg_safe[h_peaks]
peak_sig   = z[h_peaks]  # or z_poisson[h_peaks]

#---------------------

bkg_stix = median_filter(stix_count, size=window)
bkg_safe_ = np.clip(bkg_stix, 1e-6, None)
sigma_ = 1.253*np.sqrt(bkg_safe_*60)/60
z_ = (stix_count - bkg_safe_) / (sigma_/np.sqrt(window))    # "Gaussian-equivalent" significance per bin
sigma_thresh = 3  # or 4, 5, etc.
min_separation_bins = 2  # minimum distance between distinct peaks

stix_peaks, props = find_peaks( z_,
    height=sigma_thresh,    # only bins with z >= sigma_thresh
    distance=min_separation_bins)

stix_pt = (stix_t[stix_peaks].astype('datetime64[s]')).astype(str)
stix_pk_cts = stix_count[stix_peaks]
stix_pk_bkg   = bkg_safe_[stix_peaks]
stix_pk_sig   = z_[stix_peaks]  # or z_poisson[h_peaks]
#---------------------

print('X-ray peaks: ',len(peak_times_))

fig,ax=plt.subplots(figsize=(12,6))
ax2=plt.twinx()

dt=np.array(t[peaks].astype('datetime64[s]')).astype(str)
pk=np.array(x[peaks])


ax.errorbar(ht,h_count,yerr=h_er,label='HEL1OS')
ax.errorbar(stix_t,stix_count,yerr=stix_er,label='STIX')
# ax2.plot(t,flux_smooth,'m')
# ax.errorbar(ht,bkg,yerr=np.sqrt(bkg_safe*60)/(60*np.sqrt(window)),fmt='g--',label='Background')
ax2.plot(t,x,'ro-',markersize=3,label='SUIT') 
ax2.set_ylim(1e4,1e12)
ax.set_ylim(1e-4,1e4)

for i in range(len(h_peaks)):
    plt.axvline(ht[h_peaks[i]], color='b',alpha=0.2)

for i in range(len(stix_peaks)):
    plt.axvline(stix_t[stix_peaks[i]], color='g',alpha=0.2)

for i in range(len(peaks)):
    plt.axvline(t[peaks[i]], color='r',alpha=0.2)

np.savetxt('helios_peaks.csv',np.c_[peak_times_,peak_counts],header='date_time,helio_count',comments='',delimiter=',',fmt='%s')
np.savetxt('suit_diff_peaks.csv',np.c_[dt,pk],header='date_time,suit_count',comments='',delimiter=',',fmt='%s')
np.savetxt('stix_peaks.csv',np.c_[peak_times_,peak_counts],header='date_time,helio_count',comments='',delimiter=',',fmt='%s')

# plt.title(f'{case_id}')
ax.set_yscale('log')
ax2.set_yscale('log')
plt.title('Light curves')
ax2.set_ylabel('SUIT Excess intensity')
ax.set_ylabel('X-ray count rate')
ax.set_xlabel('Time')
ax.legend()
ax2.legend()
time_formatter = mdates.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(f'suit_strong_local_peaks.png',dpi=300)
plt.show()

