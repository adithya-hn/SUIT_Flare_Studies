
''''
Created on 17 oct 2025
@author: adithya-hn

A simple example to load and plot STIX spectrum and lightcurve using SunPy and SunKit-Spex

'''


from astropy.time import Time, TimeDelta
import pandas as pd
import matplotlib.pyplot as plt
from parfive import Downloader
import astropy.units as u
from astropy.time import Time
from sunkit_spex.extern.stix import STIXLoader
from sunkit_spex.legacy.fitting.fitter import Fitter, load
from datetime import datetime, timedelta
# import astropy.time as atime

#---------------Input parameters----------------

spec_file="../../data/stix/stx_spectrum_2410315184.fits"
srm_file="../../data/stix/stx_srm_2410315184.fits"
Start_t="2024-11-01T00:29:30"
End_t  ="2024-11-01T00:30:30"
start_background_time = "2024-11-01T01:08:00"
end_background_time   = "2024-11-01T01:09:00"
case='7_Nov01'
event_id=2

#-----------------------------------------------

time_profile_size = (9, 6)
spec_plot_size = (14, 8)
joint_spec_plot_size = (25, 10)
tol = 1e-20
spec_font_size = 18
xlims, ylims = [3, 100], [1e-1, 1e6]
plt.rcParams["font.size"] = spec_font_size
stix_spec = STIXLoader(spectrum_file=spec_file, srm_file=srm_file)

# plt.figure(layout="tight",figsize=spec_plot_size)

# the line that actually plots
# stix_spec.lightcurve(energy_ranges=[[4, 15], [15, 30]])
# plt.xticks(rotation=25)
# plt.savefig("stix_lightcurve.png", dpi=300)
# plt.show()

stix_spec.update_event_times(start=Time(Start_t), end=Time(End_t))
stix_spec.update_background_times(Time(start_background_time), Time(end_background_time))


print('Fitting thermal part..')
fitter = Fitter(stix_spec)
fitter.model = "(f_vth)"
fitter.loglikelihood = "cstat"
fitter.show_params

fitter.energy_fitting_range = [6,9]
fitter.params["T1_spectrum1"] = {"Status":"free", "Value":6, "Bounds":(1, 30)}
fitter.params["EM1_spectrum1"] = {"Status":"free", "Value":0.8, "Bounds":(0.1, 100)}
stix_spec_fit = fitter.fit(tol=tol)
fitter1 = Fitter(stix_spec)

print('Fitting non-thermal part..')

fitter1.model = "(f_vth+thick_fn)"
fitter1.loglikelihood = "cstat"
fitter1.show_params

fitter1.energy_fitting_range = [9,25]
fitter1.params["T1_spectrum1"] = {"Status":"fix", "Value":fitter.params["T1_spectrum1"].Value, "Bounds":(1, 30)}
fitter1.params["EM1_spectrum1"] = {"Status":"fix", "Value":fitter.params["EM1_spectrum1"].Value, "Bounds":(0.1, 100)}
fitter1.params["total_eflux1_spectrum1"] = {"Status": "free", "Value": 2, "Bounds": (1e-1, 100)}
fitter1.params["index1_spectrum1"] = {"Status": "free", "Value": 2, "Bounds": (1e-1, 15)}
fitter1.params["e_c1_spectrum1"] = {"Status": "free", "Value": 8, "Bounds": (1e-1, 1e2)}
stix_spec_fit = fitter1.fit(tol=tol) 
fitter = Fitter(stix_spec)


print('Done with initial guessings. satrting final fitting..')

fitter.model = "(f_vth+thick_fn)"
fitter.loglikelihood = "cstat"
fitter.show_params

fitter.energy_fitting_range = [6,25]
fitter.params["T1_spectrum1"] = {"Status":"free", "Value":fitter1.params["T1_spectrum1"].Value, "Bounds":(1, 30)}
fitter.params["EM1_spectrum1"] = {"Status":"free", "Value":fitter1.params["EM1_spectrum1"].Value, "Bounds":(0.1, 100)}
fitter.params["total_eflux1_spectrum1"] = {"Status": "free", "Value": fitter1.params["total_eflux1_spectrum1"].Value, "Bounds": (1e-1, 100)}
fitter.params["index1_spectrum1"] = {"Status": "free", "Value": fitter1.params["index1_spectrum1"].Value, "Bounds": (1e-1, 15)}
fitter.params["e_c1_spectrum1"] = {"Status": "free", "Value": fitter1.params["e_c1_spectrum1"].Value, "Bounds": (1e-1, 1e2)}
stix_spec_fit = fitter.fit(tol=tol)


plt.figure(figsize=(12,8))
axes, res_axes = fitter.plot()
axes[0].set_xlim([5,30])
plt.savefig(f'{case}_stix_preflarePeak_{event_id}.png',dpi=300)
plt.show()

spec_mcmc = fitter.run_mcmc(number_of_walkers=10, steps_per_walker=1200,)
fitter.burn_mcmc = 250

plt.figure(figsize=(12,8))
axes, res_axes = fitter.plot()
for a in axes:
    axes[0].set_xlim([5,30])

plt.savefig(f'{case}_peak_{event_id}_with_mcmc.png',dpi=300)
plt.show()
