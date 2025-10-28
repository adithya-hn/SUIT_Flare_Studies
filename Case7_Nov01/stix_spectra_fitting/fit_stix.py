
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


#-----------------------------------------------

time_profile_size = (9, 6)
spec_plot_size = (14, 8)
joint_spec_plot_size = (25, 10)
tol = 1e-20
spec_font_size = 18
xlims, ylims = [3, 100], [1e-1, 1e6]

plt.rcParams["font.size"] = spec_font_size

stix_spec = STIXLoader(spectrum_file=spec_file, srm_file=srm_file)
plt.figure(layout="tight",figsize=spec_plot_size)

# the line that actually plots
stix_spec.lightcurve(energy_ranges=[[4, 15], [15, 30]])
plt.xticks(rotation=25)
plt.savefig("stix_lightcurve.png", dpi=300)
plt.show()


Start_t="2024-11-01T01:44:30"
End_t="2024-11-01T01:47:00"
event_id=1
stix_spec.update_event_times(start=Time(Start_t), end=Time(End_t))
# Alternatively, you can select  the start and end event times in separate lines. e.g.
# stix_spec.start_event_time = "2024-10-01T22:10:10"
# stix_spec.end_event_time = "2024-10-01T22:10:18"


plot_range=[datetime.fromisoformat(Start_t)- timedelta(minutes=10),datetime.fromisoformat(End_t)+ timedelta(minutes=10)]
plt.figure(layout="tight")
stix_spec.lightcurve(energy_ranges=[[4, 15], [15, 30]])
plt.xlim(plot_range[0], plot_range[1])
plt.xticks(rotation=45)
plt.savefig(f"stix_event_{event_id}_{Start_t}_{End_t}.png", dpi=300)
plt.show()

fitter = Fitter(stix_spec)

fitter.model = "(f_vth+thick_fn)"
fitter.loglikelihood = "cstat"
fitter.show_params

fitter.params["T1_spectrum1"] = {"Status":"free", "Value":10, "Bounds":(1, 100)}
#fitter.params["EM1_spectrum1"] = {"Status":"free", "Value":28100, "Bounds":(10, 100000)}
fitter.energy_fitting_range = [10, 30]
fitter.params["T1_spectrum1"] = "free"
fitter.params["EM1_spectrum1"] = "free"
hel1_spec_fit = fitter.fit(tol=tol, options={"maxiter": 5000})

plt.figure()
axes, res_axes = fitter.plot()
axes[0].set_xlim([9,20])
plt.savefig('case7_nov1_stix_preflarePeak1.png',dpi=300)
plt.show()

spec_mcmc = fitter.run_mcmc(number_of_walkers=4, steps_per_walker=1200,)
fitter.burn_mcmc = 250

axes, res_axes = fitter.plot()
for a in axes:
    a.set_xlim([9,20])

plt.savefig('with_mcmc.png',dpi=300)
plt.show()
