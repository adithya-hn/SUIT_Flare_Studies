# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 11:24:41 2026

@author: srika
"""

import os
from astropy.io import fits
import glob
from datetime import datetime,timedelta
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter
from scipy.signal import find_peaks

e1 = glob.glob('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/helios/2024/07/10/HLS_20240710_000051_43141sec_lev1_V211/events/*.fits');

aa = fits.open(e1[0]);
c1 = aa[1].data;        # Cdte1 event data
c2 = aa[2].data;    

t_start = pd.to_datetime(c1['utc-isot'][0]); t_end =  pd.to_datetime(c1['utc-isot'][-1]);


#norm_v=norm[0]


tvec = pd.date_range(t_start,t_end,freq='1min');     # make a time vector with 1 minute binning
# tvec = pd.date_range(t_start,t_end,freq='20s');
obtvec = np.linspace(c1.hlsobt[0],c1.hlsobt[-1],len(tvec) );






# select events that are in a particular energy: for bandwise lightcurve
emin = 10; emax = 30;
# emin = 5; emax = 90;

fc1 = c1[np.logical_and(c1.ener>emin,c1.ener<emax)];  # cdte1 events
fc2 = c2[np.logical_and(c2.ener>emin,c2.ener<emax)];  # cdte2 events

# convert the utc datestrings to datetime:
#evt_time1 = np.array([pd.to_datetime(x) for x in fc1['utc-isot']]);
# the above operation takes a lot of time to convert each datestring to datetime
# plot to check that the time-correlation from UTC to HLSOBT is done properly:
# plt.figure();
# plt.plot(evt_time1,fc1.hlsobt,'.');
# plt.plot(tvec,obtvec,'.');
# plt.show()

# # thus use the obtvec that maps to the UTC time-bins:
hc1 = np.histogram(fc1.hlsobt,obtvec)[0];
hc2 = np.histogram(fc2.hlsobt,obtvec)[0];
cad = np.mean(np.diff(obtvec));
tmid = tvec[0:-1] + timedelta(seconds=cad*0.5); # generate the mid-point of the time intervals for the lightcurve

# plt.figure();
# # plt.plot(tmid,hc1);
# plt.errorbar(tmid,hc1,np.sqrt(hc1));
# plt.errorbar(tmid,hc2,np.sqrt(hc2));
# plt.xlabel('Time-UTC',fontsize=14);plt.ylabel(f'Counts in {np.floor(cad)}s',fontsize=14);
# plt.tick_params(labelsize=14);
# plt.show()

pd.DataFrame({'time':tmid,'cdte_counts':(hc1+hc2)/cad,'cdte_count_er':np.sqrt(hc1+hc2)/cad}).to_csv('helios_cdte_bandlc_1min_cad.csv',index=False)

plt.figure();
# plt.plot(tmid,hc1);
plt.errorbar(tmid,(hc1+hc2)/cad,np.sqrt(hc1+hc2)/cad);

plt.xlabel('Time-UTC',fontsize=14);plt.ylabel('Counts/s ',fontsize=14);
plt.tick_params(labelsize=14);
plt.show()


helios_lc=pd.read_csv('helios_cdte_bandlc_1min_cad.csv');
ht = np.array(helios_lc['time'], dtype='datetime64[s]');
h_count = np.array(helios_lc['cdte_counts'],dtype=float);
h_er = np.array(helios_lc['cdte_count_er'],dtype=float);

# h_lc=pd.read_csv('LightCurve.csv')
# lc_t=np.array(h_lc["Time (UTC)"], dtype='datetime64[s]');
# lc_count=np.array(h_lc["10.00-30.00 keV"],dtype=float)

t_st =  np.datetime64("2024-07-10T03:59:00")
t_ed =  np.datetime64("2024-07-10T06:20:00")
mask = (ht >= t_st) & (ht <= t_ed)
t_cut = ht[mask]
h_count_cut = h_count[mask]
h_er_cut = h_er[mask]
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
window = 15  # choose based on timescale over which you expect NO flares
bkg = median_filter(counts, size=window)

# Avoid zeros in background to prevent division by zero
bkg_safe = np.clip(bkg, 1e-6, None)
sigma = np.sqrt(bkg_safe)
z = (counts - bkg_safe) / (sigma/np.sqrt(window))    # "Gaussian-equivalent" significance per bin
sigma_thresh = 1  # or 4, 5, etc.
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



# # read cdte1 spectra, spectrogram files

# s1 = glob.glob('./cdte/hel1os_cdte_spectra_cdte1.fits');
# s2 = glob.glob('./cdte/hel1os_cdte_spectra_cdte2.fits');

# ss1 = fits.open(s1[0]);  ss2 = fits.open(s2[0]);

# sp1= ss1[1].data; sp2 = ss2[1].data;
# lc1 = np.array([np.sum(x) for x in sp1.COUNTS]); lc2 = np.array([np.sum(x) for x in sp2.COUNTS]); 
# exp = np.mean(sp1.EXPOSURE);

# plt.figure(); plt.errorbar(x=sp1.TSTART,y=lc1/exp,yerr=np.sqrt(lc1)/exp); plt.errorbar(x= sp2.TSTART,y=lc2/exp,yerr=np.sqrt(lc2)/exp);
# plt.ylabel('Counts/s ',fontsize=14);
# # plt.figure();















