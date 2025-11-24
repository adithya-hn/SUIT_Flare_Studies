import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import warnings
import pandas as pd
from astropy.io import fits
from astropy.time import Time
import datetime
import sunkit_spex
from sunkit_spex.legacy.fitting.fitter import Fitter

path = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/helios/'

lc_filename = f'{path}lightcurve_cdte1.fits' #[i for i in os.listdir() if i.endswith('.lc.gz')][0]
gti_filename = f'{path}gticdte1.fits'

# Loading light curve (LC) file
lc_hdul = fits.open(lc_filename)

# Loading good time intervals (GTI) file
if(gti_filename is not None):
    gti_hdul = fits.open(gti_filename)
    gti_start_arr_UTC = Time(gti_hdul[1].data['tstart'], format='mjd').to_datetime()
    gti_stop_arr_UTC = Time(gti_hdul[1].data['tstop'], format='mjd').to_datetime()

lc_hdul.info()

plot_range=[datetime.datetime(2024,7,10,3,30),datetime.datetime(2024,7,10,6,30)]


lcurve_idx = 1 # can vary between 1 to 5

plot_times_UTC1 = np.array([datetime.datetime.strptime(i[2:], '%y-%m-%dT%H:%M:%S.%f') for i in lc_hdul[1].data['ISOT']]).astype(datetime.datetime)
plot_lcurve1 = lc_hdul[1].data['CTR']
plot_times_UTC2 = np.array([datetime.datetime.strptime(i[2:], '%y-%m-%dT%H:%M:%S.%f') for i in lc_hdul[2].data['ISOT']]).astype(datetime.datetime)
plot_lcurve2 = lc_hdul[2].data['CTR']
plot_times_UTC3 = np.array([datetime.datetime.strptime(i[2:], '%y-%m-%dT%H:%M:%S.%f') for i in lc_hdul[3].data['ISOT']]).astype(datetime.datetime)
plot_lcurve3 = lc_hdul[3].data['CTR']
plot_times_UTC4 = np.array([datetime.datetime.strptime(i[2:], '%y-%m-%dT%H:%M:%S.%f') for i in lc_hdul[4].data['ISOT']]).astype(datetime.datetime)
plot_lcurve4 = lc_hdul[4].data['CTR']
plot_times_UTC5 = np.array([datetime.datetime.strptime(i[2:], '%y-%m-%dT%H:%M:%S.%f') for i in lc_hdul[5].data['ISOT']]).astype(datetime.datetime)
plot_lcurve5 = lc_hdul[5].data['CTR']
tol = 1e-20


idx1 = (plot_times_UTC1>plot_range[0]) & (plot_times_UTC1<plot_range[1])
idx2 = (plot_times_UTC2>plot_range[0]) & (plot_times_UTC2<plot_range[1])
idx3 = (plot_times_UTC3>plot_range[0]) & (plot_times_UTC3<plot_range[1])
idx4 = (plot_times_UTC4>plot_range[0]) & (plot_times_UTC4<plot_range[1])
idx5 = (plot_times_UTC5>plot_range[0]) & (plot_times_UTC5<plot_range[1])
plt.errorbar(x=plot_times_UTC1[idx1], y=plot_lcurve1[idx1], yerr=np.sqrt(plot_lcurve1[idx1]), fmt='.', capsize=5, c='b', label='Light curve (5- 20 keV)')
#plt.errorbar(x=plot_times_UTC2[idx2], y=plot_lcurve2[idx2], yerr=np.sqrt(plot_lcurve2[idx2]), fmt='.', capsize=5, c='g', label='Light curve (20- 30 keV)')
#plt.errorbar(x=plot_times_UTC3[idx3], y=plot_lcurve3[idx3], yerr=np.sqrt(plot_lcurve3[idx3]), fmt='.', capsize=5, c='r', label='Light curve (30- 40 keV)')
#plt.errorbar(x=plot_times_UTC4[idx4], y=plot_lcurve4[idx4], yerr=np.sqrt(plot_lcurve4[idx4]), fmt='.', capsize=5, c='c', label='Light curve (40- 60 keV)')
#plt.errorbar(x=plot_times_UTC5[idx5], y=plot_lcurve5[idx5], yerr=np.sqrt(plot_lcurve5[idx5]), fmt='.', capsize=5, c='m', label='Light curve (1.8- 90 keV)',alpha=0.5)
#plt.errorbar(x=plot_times_UTC5[idx5], y=plot_lcurve2[idx2]+plot_lcurve1[idx1], yerr=np.sqrt(plot_lcurve5[idx5]), fmt='.', capsize=5, c='m', label='Light curve (5- 30 keV)')
if(gti_filename is not None):
    for n,i,j in zip(np.arange(len(gti_start_arr_UTC)), gti_start_arr_UTC, gti_stop_arr_UTC):
        plt.axvspan(i,j, color='green', alpha=0.5, label='_'*n + 'Good Time Intervals')
integration_interval = (datetime.datetime(2024,11,1,1,55,00), datetime.datetime(2024,11,1,1,58,00))
# integration_interval_2 = (datetime.datetime(2024,10,3,12,6,0), datetime.datetime(2024,10,3,12,8,40))
#plt.axvspan(*integration_interval, color='r', alpha=0.5, label='Spectrum Integration Period')
# plt.axvspan(*integration_interval_2, color='b', alpha=0.5, label='Spectrum Background')
plt.legend(loc='best')
plt.grid(); plt.yscale('log')
plt.xlim(plot_times_UTC1[idx1][0], plot_times_UTC1[idx1][-1])
plt.xlabel('Time [UT]')
plt.ylabel('Counts')
plt.xticks(rotation=90)
plt.savefig('c4_helios_lc.png',dpi=300)
plt.show()
 