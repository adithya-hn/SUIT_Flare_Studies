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


lc_filename = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/helios/lightcurve_cdte1.fits' #[i for i in os.listdir() if i.endswith('.lc.gz')][0]
gti_filename = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/helios/gticdte1.fits'

# Loading light curve (LC) file
lc_hdul = fits.open(lc_filename)

# Loading good time intervals (GTI) file
if(gti_filename is not None):
    gti_hdul = fits.open(gti_filename)
    gti_start_arr_UTC = Time(gti_hdul[1].data['tstart'], format='mjd').to_datetime()
    gti_stop_arr_UTC = Time(gti_hdul[1].data['tstop'], format='mjd').to_datetime()

lcurve_idx = 1 # can vary between 1 to 5
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
plt.savefig('light_curve_oct09_1.png',dpi=300)
plt.show(block=False)
plt.pause(1)
plt.close()

#-------------------------

idx = (plot_times_UTC>datetime.datetime(2024,11,1,0,0)) & (plot_times_UTC<datetime.datetime(2024,11,1,2,30))
plt.errorbar(x=plot_times_UTC[idx], y=plot_lcurve[idx], yerr=np.sqrt(plot_lcurve[idx]), fmt='.', capsize=5, c='b', label='Light curve');
if(gti_filename is not None):
    for n,i,j in zip(np.arange(len(gti_start_arr_UTC)), gti_start_arr_UTC, gti_stop_arr_UTC):
        plt.axvspan(i,j, color='green', alpha=0.5, label='_'*n + 'Good Time Intervals')
integration_interval = (datetime.datetime(2024,11,1,1,31), datetime.datetime(2024,11,1,1,33)) # not clear after this step

plt.axvspan(*integration_interval, color='r', alpha=0.5, label='Spectrum Integration Period')
plt.legend(loc='upper right')
plt.grid(); plt.yscale('log')
plt.xlim(plot_times_UTC[idx][0], plot_times_UTC[idx][-1])
plt.xlabel('Time [UT]')
plt.ylabel('Counts')
plt.xticks(rotation=45)
plt.savefig('selected range.png',dpi=300)
plt.show()

#-------------------------

pi_filename = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/helios/hel1os_cdte_spectra_cdte1.fits'
gti_filename = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/helios/gticdte1.fits'
pi_hdul = fits.open(pi_filename)
pi_hdul.info()

plt.imshow(pi_hdul[1].data['COUNTS'].T, cmap='hot', aspect='auto', vmin=0, vmax=5, extent=[pi_hdul[1].data['TSTART'][0], pi_hdul[1].data['TSTART'][-1],511,0]); plt.colorbar();
plt.gca().invert_yaxis()
plt.xlabel('Relative Time [s]')
plt.ylabel('Channel Number')
plt.title(f'{pi_filename} Spectrogram')
plt.show(block=False)
plt.pause(1)
plt.close()

pi_times_MJD = (pi_hdul[1].data['TSTART']/86400) + pi_hdul[1].header['TSTART']
cadence = np.mean(pi_hdul[1].data['TSTOP'] - pi_hdul[1].data['TSTART'])
pi_idx = (pi_times_MJD>Time(integration_interval[0].isoformat(), format='isot').mjd) & (pi_times_MJD<Time(integration_interval[1].isoformat(), format='isot').mjd)
int_time = len(np.where(pi_idx)[0])*cadence
spec = np.sum(pi_hdul[1].data['COUNTS'][pi_idx,:], axis=0)
spec_err = np.sqrt(spec)
channels = np.arange(len(spec))

cdte_arf_filename = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/helios/hel1os_cdte_arf_v03.fits'
cdte_arf_hdul = fits.open(cdte_arf_filename)
inp_enes = np.mean(np.array([cdte_arf_hdul[1].data['ENERG_LO'], cdte_arf_hdul[1].data['ENERG_HI']]).T, axis=1)
cdte_srf_filename = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/helios/hel1os_cdte_srf_v03.fits'
cdte_srf_hdul = fits.open(cdte_srf_filename)
enes = np.mean(np.array([cdte_srf_hdul[2].data['E_MIN'],cdte_srf_hdul[2].data['E_MAX']]).T, axis=1)

plt.errorbar(x=enes, y=spec/int_time, yerr=spec_err/int_time, fmt='.-', capsize=5, c='k')
plt.yscale('log')
plt.xlabel('Energy [keV]')
plt.ylabel('Count rate [cps]')
plt.title('Spectrum')
plt.grid()
plt.show()

# plot the spectral redistribution matrix for CdTe
plt.figure()
plt.imshow(np.log10(cdte_srf_hdul[1].data['MATRIX']), extent=[cdte_srf_hdul[2].data['E_MIN'][0], cdte_srf_hdul[2].data['E_MAX'][-1], cdte_srf_hdul[1].data['ENERG_HI'][-1], cdte_srf_hdul[1].data['ENERG_LO'][0]],cmap='coolwarm');
plt.gca().invert_yaxis()
plt.colorbar()
plt.xlabel('Ouput Energy [keV]');
plt.ylabel('Input Energy [keV]')
plt.title('CdTe SRF')
plt.close()

# plot the effective area curve for CdTe
plt.figure()
plt.semilogy(inp_enes,cdte_arf_hdul[1].data['SPECRESP'], linewidth=3)
plt.xlabel('Energy (keV)')
plt.ylabel('Effective area (cm$^2$)')
plt.title('CdTe Effective area')
plt.xlim([0,60]);
plt.grid();
plt.yscale('log')
plt.close()

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
               "srm":srm}; #

# load the dict to instantiate the Fitter class:
custom_spec = Fitter(custom_dict)
zz = custom_spec.data.loaded_spec_data['spectrum1']._loaded_spec_data

# select the solar model that is appropriate
custom_spec.model = "f_vth"
custom_spec.loglikelihood= 'chi2';

# select the energy range to fit
custom_spec.energy_fitting_range = [11,35]
custom_spec.params["T1_spectrum1"] = {"Value":15, "Bounds":(10, 50)}
custom_spec.params["EM1_spectrum1"] = {"Value":0.1 ,"Bounds":(1, 1000)}

minimiser_results = custom_spec.fit()

print(custom_spec.params)

#Do some plot formatting
spec_plot_size = (10,6)
spec_font_size = 12
default_text = 10

plt.figure(figsize=spec_plot_size)
axes, res_axes = custom_spec.plot()
axes[0].set_xlim([11,40])
axes[0].set_ylim([1e-2,1e3])
plt.show()

# update the model to include a non-thermal component

T_fixed = 7.22069 #custom_spec.params["T1_spectrum1"]['Value']
EM_fixed = 374.992 #custom_spec.params["EM1_spectrum1"]['Value']

custom_spec.update_model = "(f_vth+thick_fn)"
custom_spec.energy_fitting_range = [11,35]
custom_spec.params["T1_spectrum1"] = {"Value":T_fixed, "Bounds":(7.21175, 7.22909)} # 26.06
custom_spec.params["EM1_spectrum1"] = {"Value":EM_fixed,"Bounds":(355.273, 396.135)} # 1.48
#custom_spec.params.param_status['T1_spectrum1'] = 'frozen'
#custom_spec.params.param_status['EM1_spectrum1'] = 'frozen'

custom_spec.params["total_eflux1_spectrum1"] = {"Value":91, "Bounds":(1e-3, 500)} # 0.05
custom_spec.params["index1_spectrum1"] = {"Value":8.5,"Bounds":(1, 70)} # 10
custom_spec.params["e_c1_spectrum1"] = {"Value":9.2, "Bounds":(1e-3, 50)} # 0.46

minimiser_results2 = custom_spec.fit()
print(custom_spec.params)

plt.figure(figsize=spec_plot_size)
axes, res_axes = custom_spec.plot()
axes[0].set_xlim([11,40])
axes[0].set_ylim([1e-3,1e3])
plt.savefig('spectrum_fit_oct09_1.png',dpi=300)
plt.show()