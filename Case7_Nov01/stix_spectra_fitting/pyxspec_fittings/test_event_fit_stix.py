
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
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
# set_pub_style()

#---------------Input parameters----------------

spec_file="stx_spectrum_2410315184.fits"
srm_file="stx_srm_2410315184.fits"
Start_t="2024-11-01T01:22:30"
End_t="2024-11-01T01:23:00"
case='7_Nov01'
event_id=11

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

# fitter.model = "(f_vth+f_vth+thick_fn)" #thick_fn f_vth
# fitter.loglikelihood = "cstat"
# fitter.show_params

# fitter.energy_fitting_range = [4,12]
# # sort model parameters
# fitter.params["T1_spectrum1"] = {"Value": 5, "Bounds": (1, 20)}
# fitter.params["EM1_spectrum1"] = {"Value": 2e1, "Bounds": (1e0, 800)}
# fitter.params["T2_spectrum1"] = {"Value": 12, "Bounds": (10, 60)}
# fitter.params["EM2_spectrum1"] = {"Value": 500, "Bounds": (400, 1200)}
# fitter.params["total_eflux1_spectrum1"] = {"Status": "fix", "Value": 9, "Bounds": (1e-2, 1e1)}
# fitter.params["index1_spectrum1"] = {"Status": "fix", "Value": 8, "Bounds": (1e-1, 10)}
# fitter.params["e_c1_spectrum1"] = {"Status": "fix", "Value": 10, "Bounds": (1e-1, 1e2)}
# stix_spec_fit = fitter.fit(tol=tol)

plt.figure()
axes, res_axes = fitter.plot()
axes[0].set_xlim([3,40])
plt.show()


fitter.energy_fitting_range = [12, 18]

# sort model parameters
fitter.params["T1_spectrum1"] = "fix"
fitter.params["EM1_spectrum1"] = "fix"
fitter.params["T2_spectrum1"] = "fix"
fitter.params["EM2_spectrum1"] = "fix"
fitter.params["total_eflux1_spectrum1"] = "free"
fitter.params["index1_spectrum1"] = "free"
fitter.params["e_c1_spectrum1"] = "free"
stix_spec_fit = fitter.fit(tol=tol)

plt.figure()
axes, res_axes = fitter.plot()
axes[0].set_xlim([3,30])
plt.show()

# # define energy fitting range
# fitter.energy_fitting_range = [4, 18]

# # sort model parameters
# fitter.params["T1_spectrum1"] = "free"
# fitter.params["EM1_spectrum1"] = "free"
# fitter.params["T2_spectrum1"] = "free"
# fitter.params["EM2_spectrum1"] = "free"
# fitter.params["total_eflux1_spectrum1"] = "free"
# fitter.params["index1_spectrum1"] = "free"
# fitter.params["e_c1_spectrum1"] = "free"

# stix_spec_fit = fitter.fit(tol=tol, options={"maxiter": 5000})

# plt.figure()
# axes, res_axes = fitter.plot()
# axes[0].set_xlim([4,30])
# plt.savefig(f'{case}_stix_preflarePeak_{event_id}.png',dpi=300)
# plt.show()

# spec_mcmc = fitter.run_mcmc(number_of_walkers=14, steps_per_walker=1200,)
# fitter.burn_mcmc = 250

plt.figure()
axes, res_axes = fitter.plot()
for a in axes:
    a.set_xlim([3,30])

plt.savefig(f'{case}_peak_{event_id}_with_mcmc.png',dpi=300)
plt.show()
