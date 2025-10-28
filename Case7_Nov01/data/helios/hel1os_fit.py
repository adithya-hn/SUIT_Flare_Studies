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

path = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/helios/'

lc_filename = f'{path}lightcurve_cdte2.fits' #[i for i in os.listdir() if i.endswith('.lc.gz')][0]
gti_filename = f'{path}gticdte2.fits'

# Loading light curve (LC) file
lc_hdul = fits.open(lc_filename)

# Loading good time intervals (GTI) file
if(gti_filename is not None):
    gti_hdul = fits.open(gti_filename)
    gti_start_arr_UTC = Time(gti_hdul[1].data['tstart'], format='mjd').to_datetime()
    gti_stop_arr_UTC = Time(gti_hdul[1].data['tstop'], format='mjd').to_datetime()

lc_hdul.info()

plot_range=[datetime.datetime(2024,11,1,1,30),datetime.datetime(2024,11,1,1,59)]

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
plt.errorbar(x=plot_times_UTC2[idx2], y=plot_lcurve2[idx2], yerr=np.sqrt(plot_lcurve2[idx2]), fmt='.', capsize=5, c='g', label='Light curve (20- 30 keV)')
plt.errorbar(x=plot_times_UTC3[idx3], y=plot_lcurve3[idx3], yerr=np.sqrt(plot_lcurve3[idx3]), fmt='.', capsize=5, c='r', label='Light curve (30- 40 keV)')
plt.errorbar(x=plot_times_UTC4[idx4], y=plot_lcurve4[idx4], yerr=np.sqrt(plot_lcurve4[idx4]), fmt='.', capsize=5, c='c', label='Light curve (40- 60 keV)')
#plt.errorbar(x=plot_times_UTC5[idx5], y=plot_lcurve5[idx5], yerr=np.sqrt(plot_lcurve5[idx5]), fmt='.', capsize=5, c='m', label='Light curve (80- 90 keV)',alpha=0.5)
if(gti_filename is not None):
    for n,i,j in zip(np.arange(len(gti_start_arr_UTC)), gti_start_arr_UTC, gti_stop_arr_UTC):
        plt.axvspan(i,j, color='green', alpha=0.5, label='_'*n + 'Good Time Intervals')
integration_interval = (datetime.datetime(2024,11,1,1,55,00), datetime.datetime(2024,11,1,1,58,00))
# integration_interval_2 = (datetime.datetime(2024,10,3,12,6,0), datetime.datetime(2024,10,3,12,8,40))
plt.axvspan(*integration_interval, color='r', alpha=0.5, label='Spectrum Integration Period')
# plt.axvspan(*integration_interval_2, color='b', alpha=0.5, label='Spectrum Background')
plt.legend(loc='best')
plt.grid(); plt.yscale('log')
plt.xlim(plot_times_UTC1[idx1][0], plot_times_UTC1[idx1][-1])
plt.xlabel('Time [UT]')
plt.ylabel('Counts')
plt.xticks(rotation=90)
plt.savefig('c7_helios_lc.png',dpi=300)
plt.show()

pi_filename = f'{path}hel1os_cdte_spectra_cdte1.fits'
gti_filename = f'{path}gticdte1.fits'
pi_hdul = fits.open(pi_filename)
pi_hdul.info()

pi_times_MJD = (pi_hdul[1].data['TSTART']/86400) + pi_hdul[1].header['TSTART']
cadence = np.mean(pi_hdul[1].data['TSTOP'] - pi_hdul[1].data['TSTART'])

pi_idx = (pi_times_MJD>Time(integration_interval[0].isoformat(), format='isot').mjd) & (pi_times_MJD<Time(integration_interval[1].isoformat(), format='isot').mjd)
int_time = len(np.where(pi_idx)[0])*cadence

# pi_idx_2 = (pi_times_MJD>Time(integration_interval_2[0].isoformat(), format='isot').mjd) & (pi_times_MJD<Time(integration_interval_2[1].isoformat(), format='isot').mjd)
# int_time_2 = len(np.where(pi_idx_2)[0])*cadence

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

srm = np.zeros_like(cdte_srf_hdul[1].data['MATRIX']);    # spectral response matrix that is combination of ARF and SRF
for k in range(len(cdte_srf_hdul[1].data['MATRIX'])):
    srm[k,:] = cdte_arf_hdul[1].data['SPECRESP'][k]*cdte_srf_hdul[1].data['MATRIX'][k,:]

custom_dict = {"count_channel_bins":out_chan_bins, #
               "count_channel_binning":out_chan_widths, #
               "counts":spec, # in counts
               "count_error":spec_err, # in counts
               "count_rate":spec/(int_time*out_ene_binwidth), # in counts/s/keV
               "count_rate_error":np.sqrt(spec)/(int_time*out_ene_binwidth), # in counts/s/keV
               "photon_channel_bins":inp_chan_bins, #
               "photon_channel_binning":inp_chan_widths, #
               "photon_channel_mids":np.mean(inp_chan_bins,axis=1), #
               "effective_exposure":int_time, #
               "srm":srm,}; #

fitter = Fitter(custom_dict)
fitter.model = "(f_vth)"#+thick_fn)"
fitter.loglikelihood = "cstat"
# fitter.loglikelihood = "loglike"
# fitter.energy_fitting_range = [8,15]

fitter.params["T1_spectrum1"] = {"Status":"free", "Value":10, "Bounds":(1, 100)}
fitter.params["EM1_spectrum1"] = {"Status":"free", "Value":28100, "Bounds":(10, 100000)}

# fitter.params["total_eflux1_spectrum1"] = {"Status": "fix", "Value": 0.9, "Bounds": (1e-1, 1e1)}
# fitter.params["index1_spectrum1"] = {"Status": "fix", "Value": 6, "Bounds": (3, 1e2)}
# fitter.params["e_c1_spectrum1"] = {"Status": "fix", "Value": 1.1e1, "Bounds": (1.1e1, 2e2)}

# hel1_spec_fit = fitter.fit(tol=tol, options={"maxiter": 1000})

# fitter.energy_fitting_range = [10, 30]

# # sort model parameters
# fitter.params["T1_spectrum1"] = "fix"
# fitter.params["EM1_spectrum1"] = "fix"
# fitter.params["total_eflux1_spectrum1"] = "free"
# fitter.params["index1_spectrum1"] = "free"
# fitter.params["e_c1_spectrum1"] = "free"

# hel1_spec_fit = fitter.fit(tol=tol, options={"maxiter": 1000})

fitter.energy_fitting_range = [10, 15]

# sort model parameters
fitter.params["T1_spectrum1"] = "free"
fitter.params["EM1_spectrum1"] = "free"
# fitter.params["total_eflux1_spectrum1"] = "free"
# fitter.params["index1_spectrum1"] = "free"
# fitter.params["e_c1_spectrum1"] = "free"

hel1_spec_fit = fitter.fit(tol=tol, options={"maxiter": 5000})

plt.figure()
axes, res_axes = fitter.plot()
axes[0].set_xlim([9,20])
# axes[0].text(0.05, 0.95, f"Red. Chi2 = {red_chi2:.2f}", transform=axes[0].transAxes, fontsize=12, verticalalignment='top')
plt.savefig('case7_nov1_preflarePeak.png',dpi=300)
plt.show()

spec_mcmc = fitter.run_mcmc(number_of_walkers=4, steps_per_walker=1200,)
fitter.burn_mcmc = 250

axes, res_axes = fitter.plot()
for a in axes:
    a.set_xlim([9,20])

plt.savefig('with_mcmc.png',dpi=300)
plt.show()