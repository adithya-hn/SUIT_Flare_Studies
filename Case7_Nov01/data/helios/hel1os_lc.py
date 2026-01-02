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

path = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/helios/'

lc_filename1 = f'{path}lightcurve_cdte2.fits' #[i for i in os.listdir() if i.endswith('.lc.gz')][0]
gti_filename1 = f'{path}gticdte2.fits'

lc_filename2 = f'{path}lightcurve_cdte2.fits' #[i for i in os.listdir() if i.endswith('.lc.gz')][0]
gti_filename2 = f'{path}gticdte2.fits'


# Loading light curve (LC) file
lc_hdul1 = fits.open(lc_filename1)
lc_hdul2 = fits.open(lc_filename2)

# Loading good time intervals (GTI) file
if(gti_filename1 is not None):
    gti_hdul = fits.open(gti_filename1)
    gti_start_arr_UTC = Time(gti_hdul[1].data['tstart'], format='mjd').to_datetime()
    gti_stop_arr_UTC = Time(gti_hdul[1].data['tstop'], format='mjd').to_datetime()

#lc_hdul1.info()

plt.figure(figsize=(14,8))
plot_range=[datetime.datetime(2024,11,1,0,16),datetime.datetime(2024,11,1,2,31)]


plot_times_UTC1 = np.array([datetime.datetime.strptime(i[2:], '%y-%m-%dT%H:%M:%S.%f') for i in lc_hdul1[1].data['ISOT']]).astype(datetime.datetime)
plot_lcurve1 = lc_hdul1[1].data['CTR']
plot_times_UTC2 = np.array([datetime.datetime.strptime(i[2:], '%y-%m-%dT%H:%M:%S.%f') for i in lc_hdul1[2].data['ISOT']]).astype(datetime.datetime)
plot_lcurve2 = lc_hdul1[2].data['CTR']
plot_times_UTC3 = np.array([datetime.datetime.strptime(i[2:], '%y-%m-%dT%H:%M:%S.%f') for i in lc_hdul1[3].data['ISOT']]).astype(datetime.datetime)
plot_lcurve3 = lc_hdul1[3].data['CTR']
plot_times_UTC4 = np.array([datetime.datetime.strptime(i[2:], '%y-%m-%dT%H:%M:%S.%f') for i in lc_hdul1[4].data['ISOT']]).astype(datetime.datetime)
plot_lcurve4 = lc_hdul1[4].data['CTR']
plot_times_UTC5 = np.array([datetime.datetime.strptime(i[2:], '%y-%m-%dT%H:%M:%S.%f') for i in lc_hdul1[5].data['ISOT']]).astype(datetime.datetime)
plot_lcurve5 = lc_hdul1[5].data['CTR']



idx1 = (plot_times_UTC1>plot_range[0]) & (plot_times_UTC1<plot_range[1])
idx2 = (plot_times_UTC2>plot_range[0]) & (plot_times_UTC2<plot_range[1])
idx3 = (plot_times_UTC3>plot_range[0]) & (plot_times_UTC3<plot_range[1])
idx4 = (plot_times_UTC4>plot_range[0]) & (plot_times_UTC4<plot_range[1])
idx5 = (plot_times_UTC5>plot_range[0]) & (plot_times_UTC5<plot_range[1])

plot_times_UTC2_ = np.array([datetime.datetime.strptime(i[2:], '%y-%m-%dT%H:%M:%S.%f') for i in lc_hdul2[1].data['ISOT']]).astype(datetime.datetime)
plot_lcurve2_ = lc_hdul2[1].data['CTR']

times1  = pd.to_datetime(plot_times_UTC1)
counts1 = np.asarray(plot_lcurve1, dtype=float)
df1 = pd.DataFrame({'c1': counts1}, index=times1)

# Second file
times2  = pd.to_datetime(plot_times_UTC2_)
counts2 = np.asarray(plot_lcurve2_, dtype=float)
df2 = pd.DataFrame({'c2': counts2}, index=times2)

df1_1min = df1.resample('1min').sum()
df2_1min = df2.resample('1min').sum()
df_join = df1_1min.join(df2_1min, how='outer')

df_join['total'] = df_join['c1'].fillna(0) + df_join['c2'].fillna(0)
df_join = df_join.loc[plot_range[0]:plot_range[1]]


#plt.errorbar(x=plot_times_UTC1[idx1], y=plot_lcurve1[idx1], yerr=np.sqrt(plot_lcurve1[idx1]), fmt='.', capsize=5, c='b', label='Light curve (5- 20 keV)')
#plt.errorbar(x=plot_times_UTC2[idx2], y=plot_lcurve2[idx2], yerr=np.sqrt(plot_lcurve2[idx2]), fmt='.', capsize=5, c='g', label='Light curve (20- 30 keV)')
#plt.errorbar(x=plot_times_UTC3[idx3], y=plot_lcurve3[idx3], yerr=np.sqrt(plot_lcurve3[idx3]), fmt='.', capsize=5, c='r', label='Light curve (30- 40 keV)')
#plt.errorbar(x=plot_times_UTC4[idx4], y=plot_lcurve4[idx4], yerr=np.sqrt(plot_lcurve4[idx4]), fmt='.', capsize=5, c='c', label='Light curve (40- 60 keV)')
#plt.errorbar(x=plot_times_UTC5[idx5], y=plot_lcurve5[idx5], yerr=np.sqrt(plot_lcurve5[idx5]), fmt='.', capsize=5, c='m', label='Light curve (1.8- 90 keV)',alpha=0.5)


plt.errorbar(df_join.index,df_join['total'],yerr=np.sqrt(df_join['total']),fmt='.-',capsize=5)
plt.yscale('log')
plt.close()

pi_filename = f'{path}hel1os_cdte_spectra_cdte1.fits'
gti_filename = f'{path}gticdte1.fits'
pi_hdul = fits.open(pi_filename)
#pi_hdul.info()
integration_interval = (datetime.datetime(2024,11,1,1,55,00), datetime.datetime(2024,11,1,1,58,00))
pi_times_MJD = (pi_hdul[1].data['TSTART']/86400) + pi_hdul[1].header['TSTART']

cadence = np.mean(pi_hdul[1].data['TSTOP'] - pi_hdul[1].data['TSTART'])

pi_idx = (pi_times_MJD>Time(integration_interval[0].isoformat(), format='isot').mjd) & (pi_times_MJD<Time(integration_interval[1].isoformat(), format='isot').mjd)

int_time = len(np.where(pi_idx)[0])*cadence

spec = np.sum(pi_hdul[1].data['COUNTS'][pi_idx,:], axis=0)
spec_err = np.sqrt(spec)

channels = np.arange(len(spec))

cdte_arf_filename = f'{path}CdTeResponseReader/hel1os_cdte_arf_v03.fits'
cdte_arf_hdul = fits.open(cdte_arf_filename)
inp_enes = np.mean(np.array([cdte_arf_hdul[1].data['ENERG_LO'], cdte_arf_hdul[1].data['ENERG_HI']]).T, axis=1)
cdte_srf_filename = f'{path}CdTeResponseReader/hel1os_cdte_srf_v03.fits'
cdte_srf_hdul = fits.open(cdte_srf_filename)
enes = np.mean(np.array([cdte_srf_hdul[2].data['E_MIN'],cdte_srf_hdul[2].data['E_MAX']]).T, axis=1)


out_chan_bins = np.array([cdte_srf_hdul[2].data['E_MIN'],cdte_srf_hdul[2].data['E_MAX']]).T;
inp_chan_bins = np.array([cdte_srf_hdul[1].data['ENERG_LO'],cdte_srf_hdul[1].data['ENERG_HI']]).T;
out_ene_binwidth = np.mean(cdte_srf_hdul[2].data['E_MAX'] - cdte_srf_hdul[2].data['E_MIN'])
inp_ene_binwidth = np.mean(cdte_srf_hdul[1].data['ENERG_HI'] - cdte_srf_hdul[1].data['ENERG_LO'])
out_chan_widths = np.ones(cdte_srf_hdul[1].data['MATRIX'].shape[1])*out_ene_binwidth
inp_chan_widths = np.ones(cdte_srf_hdul[1].data['MATRIX'].shape[0])*inp_ene_binwidth

energy_edges = out_chan_bins
all_counts = pi_hdul[1].data['COUNTS']    
maskE = energy_edges[:,0] >= 30.0    # example: >12 keV
counts_gt12 = np.sum(all_counts[:, maskE], axis=1)
photon_energy = np.mean(out_chan_bins, axis=1) 
idx_gt12 =  (photon_energy >=12.0) # & (photon_energy <=55.0 )

all_counts = pi_hdul[1].data['COUNTS']       # shape (Ntimes, Nchan)
lc_gt12 =np.sum(all_counts[:, idx_gt12], axis=1)
lc_gt12_err = np.sqrt(lc_gt12)


tstart = pi_hdul[1].data['TSTART']   # seconds
tstop  = pi_hdul[1].data['TSTOP']    # seconds
tmid = 0.5 * (tstart + tstop)        # mid timestamp in second
date_obs = pi_hdul[1].header['TSTART']    # string
t = Time((date_obs+tmid/86400), format='mjd')
times_utc = t.to_datetime() 
#times_utc = [t0 datetime.timedelta(seconds=float(s)) for s in tmid]
print(times_utc) 
print(pi_times_MJD)


 

all_counts = pi_hdul[1].data['COUNTS']         # shape (Nspec, 511)
lc_gt12 = np.sum(all_counts[:, idx_gt12], axis=1)
lc_gt12_err = np.sqrt(lc_gt12)

df = pd.DataFrame({'counts': lc_gt12.astype(float)}, index=pd.to_datetime(times_utc))
df_1min = df.resample('1min').sum()

plt.errorbar(df_1min.index,df_1min['counts'],yerr=np.sqrt(df_1min['counts']),fmt='.-', capsize=5)
plt.yscale('log')
plt.show()
