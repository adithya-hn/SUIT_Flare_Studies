import warnings
warnings.simplefilter('ignore')
import numpy as np
from astropy.io import fits

import datetime

import os
import glob

import matplotlib.pyplot as plt
%matplotlib inline
# plt.style.use(['science','notebook','grid'])

import matplotlib.dates as mdates
import matplotlib.ticker as mticker

plt.rcParams.update({'font.size': 18})


lc_file = 'AL1_SLX_L1_20241025_v1.1/SDD2/AL1_SOLEXS_20241025_SDD2_L1.lc.gz'

hdu_lc = fits.open(lc_file)
lc_data = hdu_lc[1].data

def timestamp2datetime(timestamp_arr):
    datetime_arr = []
    for ta in timestamp_arr:
        datetime_arr.append(datetime.datetime.fromtimestamp(ta, tz=datetime.timezone.utc))

    return datetime_arr

lc_datetime = timestamp2datetime(lc_data['TIME'])

#--------------------------------------Light curve--------------------
fig, ax = plt.subplots()
ax.plot(lc_datetime,lc_data['COUNTS'],ms=0,ls='-',lw=3)
ax.set_yscale('log')
ax.set_xlabel('Time (s)')
ax.set_ylabel(r'SoLEXS counts/s')
ax.set_title('SoLEXS Lightcurve Data');

ax.xaxis.set(major_formatter=mdates.DateFormatter(r"%H:%M"))
fig.autofmt_xdate()
plt.show()
#------------------------------rebinnig------------------------------

def rebin_lc(lc, time_arr ,rebin_sec): #lc: counts per second
    datetime_arr = time_arr
    extra_bins = len(lc) % rebin_sec
    if extra_bins != 0:
        lc = lc[:-extra_bins]
    new_bins = int(len(lc)/rebin_sec)
    new_lc = lc.reshape((new_bins, rebin_sec)).sum(axis=1)
    new_tm = np.arange(new_bins)*rebin_sec


    new_datetime_arr = []
    for ii in new_tm:
        new_datetime_arr.append(datetime_arr[int(ii)])

    new_lc = new_lc/rebin_sec 

    return new_lc, new_datetime_arr
rebin_sec = 60
lc_60, time_arr_60 = rebin_lc(lc_data['COUNTS'],lc_data['TIME'],rebin_sec)
datetime_arr_60 = timestamp2datetime(time_arr_60)

fig, ax = plt.subplots()
ax.plot(lc_datetime,lc_data['COUNTS'],ms=0,ls='-',lw=3)
ax.plot(datetime_arr_60,lc_60,ms=0,ls='-',lw=3,color='r')

ax.set_yscale('log')
ax.set_xlabel('Time (s)')
ax.set_ylabel(r'SoLEXS counts/s')
ax.set_title('SoLEXS Lightcurve Data');
ax.xaxis.set(major_formatter=mdates.DateFormatter(r"%H:%M"))
fig.autofmt_xdate()
plt.show()

#-------------------------- spectral fit prep----------------------
pi_file = 'AL1_SLX_L1_20241025_v1.1/SDD2/AL1_SOLEXS_20241025_SDD2_L1.pi.gz'
hdu_pi = fits.open(pi_file)
pi_data = hdu_pi[1].data

fig, ax = plt.subplots()

ax.step(pi_data['CHANNEL'][300],pi_data['COUNTS'][300])
ax.step(pi_data['CHANNEL'][27200],pi_data['COUNTS'][27200])

ax.set_yscale('log')
plt.show()