import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.plot(); plt.close();

matplotlib.rcParams['figure.figsize'] = (12,8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 22})

import warnings
import pandas as pd
from astropy.io import fits
from astropy.time import Time
import datetime
import sunkit_spex
from sunkit_spex.legacy.fitting.fitter import Fitter



lc_filename = './HLS_20240710_121028_42567sec_lev1_V211/cdte/lightcurve_cdte1.fits' #[i for i in os.listdir() if i.endswith('.lc.gz')][0]
gti_filename = './HLS_20240710_121028_42567sec_lev1_V211/aux/gticdte1.fits'

# Loading light curve (LC) file
lc_hdul = fits.open(lc_filename)

# Loading good time intervals (GTI) file
if(gti_filename is not None):
    gti_hdul = fits.open(gti_filename)
    gti_start_arr_UTC = Time(gti_hdul[1].data['tstart'], format='mjd').to_datetime()
    gti_stop_arr_UTC = Time(gti_hdul[1].data['tstop'], format='mjd').to_datetime()

lc_hdul.info()

lcurve_idx = 3 # can vary between 1 to 5
plot_times_UTC = np.array([datetime.datetime.strptime(i[2:], '%y-%m-%dT%H:%M:%S.%f') for i in lc_hdul[lcurve_idx].data['ISOT']]).astype(datetime.datetime)
plot_lcurve = lc_hdul[lcurve_idx].data['CTR']

fig = plt.figure()
plt.errorbar(x=plot_times_UTC, y=plot_lcurve, yerr=np.sqrt(plot_lcurve), fmt='.', capsize=5, c='b', label='Light curve');
if(gti_filename is not None):
    for n,i,j in zip(np.arange(len(gti_start_arr_UTC)), gti_start_arr_UTC, gti_stop_arr_UTC):
        plt.axvspan(i,j, color='green', alpha=0.5, label='_'*n + 'Good Time Intervals')

plt.legend(loc='upper left')
plt.grid(); plt.yscale('log')
plt.xlim(plot_times_UTC[0], plot_times_UTC[-1])
plt.xlabel('Time [UT]')
plt.ylabel('Counts')
plt.xticks(rotation=90)
plt.title(f"{lc_hdul[lcurve_idx].header['EXTNAME']}")
ax = plt.gca()
plt.show()