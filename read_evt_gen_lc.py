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

e1 = glob.glob('./events/*.fits');

aa = fits.open(e1[0]);
c1 = aa[1].data;        # Cdte1 event data
c2 = aa[2].data;    

t_start = pd.to_datetime(c1['utc-isot'][0]); t_end =  pd.to_datetime(c1['utc-isot'][-1]);

# tvec = pd.date_range(t_start,t_end,freq='1min');     # make a time vector with 1 minute binning
tvec = pd.date_range(t_start,t_end,freq='20s');
obtvec = np.linspace(c1.hlsobt[0],c1.hlsobt[-1],len(tvec) );



# select events that are in a particular energy: for bandwise lightcurve
# emin = 10; emax = 30;
emin = 5; emax = 90;

fc1 = c1[np.logical_and(c1.ener>emin,c1.ener<emax)];  # cdte1 events
fc2 = c2[np.logical_and(c2.ener>emin,c2.ener<emax)];  # cdte2 events

# convert the utc datestrings to datetime:
# evt_time1 = np.array([pd.to_datetime(x) for x in fc1['utc-isot']]);
# the above operation takes a lot of time to convert each datestring to datetime
# plot to check that the time-correlation from UTC to HLSOBT is done properly:
# plt.figure();
# plt.plot(evt_time1,fc1.hlsobt,'.');
# plt.plot(tvec,obtvec,'.');

# thus use the obtvec that maps to the UTC time-bins:
hc1 = np.histogram(fc1.hlsobt,obtvec)[0];
hc2 = np.histogram(fc2.hlsobt,obtvec)[0];
cad = np.mean(np.diff(obtvec));
tmid = tvec[0:-1] + timedelta(seconds=cad*0.5); # generate the mid-point of the time intervals for the lightcurve

plt.figure();
# plt.plot(tmid,hc1);
plt.errorbar(tmid,hc1,np.sqrt(hc1));
plt.errorbar(tmid,hc2,np.sqrt(hc2));
plt.xlabel('Time-UTC',fontsize=14);plt.ylabel(f'Counts in {np.floor(cad)}s',fontsize=14);
plt.tick_params(labelsize=14);


plt.figure();
# plt.plot(tmid,hc1);
plt.errorbar(tmid,hc1/cad,np.sqrt(hc1)/cad);
plt.errorbar(tmid,hc2/cad,np.sqrt(hc2)/cad);
plt.xlabel('Time-UTC',fontsize=14);plt.ylabel('Counts/s ',fontsize=14);
plt.tick_params(labelsize=14);


# read cdte1 spectra, spectrogram files

s1 = glob.glob('./cdte/hel1os_cdte_spectra_cdte1.fits');
s2 = glob.glob('./cdte/hel1os_cdte_spectra_cdte2.fits');

ss1 = fits.open(s1[0]);  ss2 = fits.open(s2[0]);

sp1= ss1[1].data; sp2 = ss2[1].data;
lc1 = np.array([np.sum(x) for x in sp1.COUNTS]); lc2 = np.array([np.sum(x) for x in sp2.COUNTS]); 
exp = np.mean(sp1.EXPOSURE);

plt.figure(); plt.errorbar(x=sp1.TSTART,y=lc1/exp,yerr=np.sqrt(lc1)/exp); plt.errorbar(x= sp2.TSTART,y=lc2/exp,yerr=np.sqrt(lc2)/exp);
plt.ylabel('Counts/s ',fontsize=14);
# plt.figure();















