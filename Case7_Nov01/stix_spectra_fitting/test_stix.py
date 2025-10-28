
''''
Created on 17 oct 2025
@author: adithya-hn

A simple example to load and plot STIX spectrum and lightcurve using SunPy and SunKit-Spex

'''


import matplotlib.pyplot as plt
from parfive import Downloader
import astropy.units as u
from astropy.time import Time
from sunkit_spex.extern.stix import STIXLoader
from sunkit_spex.legacy.fitting.fitter import Fitter, load
from datetime import datetime, timedelta


#---------------Input parameters----------------

spec_file="../data/stix/stx_spectrum_2410315184.fits"
srm_file="../data/stix/stx_srm_2410315184.fits"
Start_t="2024-11-01T01:54:30"
End_t="2024-11-01T01:57:00"
case='7_Nov01'
event_id=4

#-----------------------------------------------

time_profile_size = (9, 6)
spec_plot_size = (14, 8)
joint_spec_plot_size = (25, 10)
tol = 1e-20
spec_font_size = 18
xlims, ylims = [3, 100], [1e-1, 1e6]
plt.rcParams["font.size"] = spec_font_size

stix_spec = STIXLoader(spectrum_file=spec_file, srm_file=srm_file)


stix_spec.update_event_times(start=Time(Start_t), end=Time(End_t))

fitter = Fitter(stix_spec)

fitter.model = "(thick_warm)" #+ thick_fn
fitter.loglikelihood = "cstat"
fitter.show_params

fitter.energy_fitting_range = [4, 50]

# sort model parameters
# fitter.params["T1_spectrum1"] = {"Value": 11, "Bounds": (5, 15), "Status": "free"}
# fitter.params["EM1_spectrum1"] = {"Value": 470, "Bounds": (10, 800), "Status": "free"}
# # fitter.params["T2_spectrum1"] = {"Value": 15, "Bounds": (12, 20), "Status": "free"}
# # fitter.params["EM2_spectrum1"] = {"Value": 7, "Bounds": (3, 20), "Status": "free"}
fitter.params["tot_eflux1_spectrum1"] = 'free'
fitter.params["ind1_spectrum1"] = 'free'
fitter.params["ind1_spectrum1"] = 'free'
fitter.params["plasma_d1_spectrum1"] ='free'
fitter.params["loop_temp1_spectrum1"]='free'
fitter.params["length1_spectrum1"]='free'

stix_spec_fit = fitter.fit(tol=tol,options={"maxiter": 5000})

print(fitter.params)

plt.figure(figsize=spec_plot_size)

# the line that actually plots
axes, res_axes = fitter.plot()

# make plot nicer
for a in axes:
    a.set_xlim(xlims)
    a.set_ylim(ylims)
    a.set_xscale("log")
plt.savefig(f'test_non_thermal_{case}_peak_{event_id}.png',dpi=300)
plt.show()

spec_mcmc = fitter.run_mcmc(number_of_walkers=6, steps_per_walker=1200,)
fitter.burn_mcmc = 250

plt.figure()
axes, res_axes = fitter.plot()
for a in axes:
    a.set_xlim([3,30])

plt.savefig(f'test_non_thermal_{case}_peak_{event_id}_with_mcmc.png',dpi=300)
plt.show()