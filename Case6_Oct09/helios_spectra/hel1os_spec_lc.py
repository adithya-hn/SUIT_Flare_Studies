import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.time import Time
import warnings
import pandas as pd
from astropy.io import fits
from astropy.time import Time
import datetime
import sunkit_spex
from sunkit_spex.legacy.fitting.fitter import Fitter

#------------------------------------------------------------------------------------------------

path = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case9_Nov13/data/helios/'
start_dt,end_dt=datetime.datetime(2024,11,13,0,22),datetime.datetime(2024,11,13,0,37)
pi_filename1 = f'{path}hel1os_cdte_spectra_cdte1.fits'
gti_filename1 = f'{path}gticdte1.fits'

pi_filename2 = f'{path}hel1os_cdte_spectra_cdte2.fits'
gti_filename2 = f'{path}gticdte2.fits'

cdte_srf_filename = f'{path}CdTeResponseReader/hel1os_cdte_srf_v03.fits'

#------------------------------------------------------------------------------------------------

pi_hdul1 = fits.open(pi_filename1)
pi_hdul2 = fits.open(pi_filename2)
#pi_hdul.info()


cdte_srf_hdul = fits.open(cdte_srf_filename)
enes = np.mean(np.array([cdte_srf_hdul[2].data['E_MIN'],cdte_srf_hdul[2].data['E_MAX']]).T, axis=1)
out_chan_bins = np.array([cdte_srf_hdul[2].data['E_MIN'],cdte_srf_hdul[2].data['E_MAX']]).T
inp_chan_bins = np.array([cdte_srf_hdul[1].data['ENERG_LO'],cdte_srf_hdul[1].data['ENERG_HI']]).T
out_ene_binwidth = np.mean(cdte_srf_hdul[2].data['E_MAX'] - cdte_srf_hdul[2].data['E_MIN'])
inp_ene_binwidth = np.mean(cdte_srf_hdul[1].data['ENERG_HI'] - cdte_srf_hdul[1].data['ENERG_LO'])
out_chan_widths = np.ones(cdte_srf_hdul[1].data['MATRIX'].shape[1])*out_ene_binwidth
inp_chan_widths = np.ones(cdte_srf_hdul[1].data['MATRIX'].shape[0])*inp_ene_binwidth

photon_energy = np.mean(out_chan_bins, axis=1) 
idx_gt12 =  (photon_energy >=12.0) & (photon_energy <=30.0 )

cdte1_counts = pi_hdul1[1].data['COUNTS']       # shape (Ntimes, Nchan)
cdte2_counts = pi_hdul2[1].data['COUNTS']  
lc_gt12 =np.sum((cdte1_counts[2:, idx_gt12]+cdte2_counts[1:, idx_gt12]), axis=1)
lc_gt12_err = np.sqrt(lc_gt12)


tstart = pi_hdul1[1].data['TSTART']   # seconds
tstop  = pi_hdul1[1].data['TSTOP']    # seconds
tmid = 0.5 * (tstart + tstop)        # mid timestamp in second
date_obs = pi_hdul1[1].header['TSTART']    # string
t = Time((date_obs+tmid/86400), format='mjd')
times_utc = t.to_datetime() 


df = pd.DataFrame({'counts': lc_gt12.astype(float)}, index=pd.to_datetime(times_utc[2:]))
df_1min = df.resample('1min').sum()
df_cut = df_1min.loc[start_dt:end_dt]
df_cut.to_csv('binned_counts_12_30kev_lc.csv')


plt.figure(figsize=(14,8))
#plt.xlim(start_dt, end_dt)
plt.errorbar(df_cut.index,df_cut['counts'],yerr=np.sqrt(df_cut['counts']),color='k',fmt='.-', capsize=5)
plt.xlabel('Time(UT)')
plt.ylabel('Counts')
plt.title('HEL1OS counts')
plt.yscale('log')
plt.savefig('Helios_lc_from_spec_12_30.png',dpi=300)
plt.show()
