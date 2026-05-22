from astropy.time import Time, TimeDelta
import pandas as pd
import matplotlib.pyplot as plt
from parfive import Downloader
import astropy.units as u
from astropy.time import Time
from sunkit_spex.extern.stix import STIXLoader
from sunkit_spex.legacy.fitting.fitter import Fitter, load
from datetime import datetime, timedelta


flare_start = Time("2024-11-01T14:15:00")
flare_end   = Time("2024-11-01T14:40:00")
#---------------Input parameters----------------

spec_file="stx_spectrum_2411012243.fits"
srm_file="stx_srm_2411012243.fits"

case='8_Nov01'


#-----------------------------------------------


bin_size = TimeDelta(20, format='sec')
tol = 1e-20
time_bins = []
t = flare_start
while t < flare_end:
    time_bins.append((t, t+bin_size))
    t += bin_size

outfile = f"{case}_stix_timeResolved_fit.csv"

columns = [
"time_start","time_end",
"T","T_err",
"EM","EM_err",
"Flux","Flux_err",
"Index","Index_err",
"Ecut","Ecut_err"
]

df = pd.DataFrame(columns=columns)

for i,(t0,t1) in enumerate(time_bins):

    print(f"Fitting bin {i} : {t0.isot} - {t1.isot}")

    stix_spec = STIXLoader(spectrum_file=spec_file, srm_file=srm_file)

    stix_spec.update_event_times(start=t0,end=t1)

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
    title_str = f"{t0.isot}  –  {t1.isot}"
    axes[0].set_title(title_str)
    axes[0].set_xlim([5,30])
    plt.savefig(f'{case}_stix_onset_step_{i}.png',dpi=300)
    plt.close()

    spec_mcmc = fitter.run_mcmc(number_of_walkers=10, steps_per_walker=1200,)
    fitter.burn_mcmc = 250

    plt.figure(figsize=(12,8))
    axes, res_axes = fitter.plot()
    for a in axes:
        axes[0].set_xlim([5,30])
        axes[0].set_title(title_str)

    plt.savefig(f'{case}_step_{i}_with_mcmc.png',dpi=300)
    plt.close()
    print('Statistic: ',stix_spec_fit)
    tab = fitter.mcmc_table

    def get_val_err(tab, pname):

        mid  = tab[tab['Param']==pname]['Mid'][0]
        low  = tab[tab['Param']==pname]['LowB'][0]
        high = tab[tab['Param']==pname]['HighB'][0]

        err = (high - low)/2
        return mid, err
    T,     T_err     = get_val_err(tab,'T1_spectrum1')
    EM,    EM_err    = get_val_err(tab,'EM1_spectrum1')
    Flux,  Flux_err  = get_val_err(tab,'total_eflux1_spectrum1')
    Index, Index_err = get_val_err(tab,'index1_spectrum1')
    Ecut,  Ecut_err  = get_val_err(tab,'e_c1_spectrum1')

    print(f'{i} fit vals: ',T, EM)
    # plot_cstat = stix_spec_fit[0].statistic
    # print("Fit Cstat =", plot_cstat)


    
    # log_likelihood = stix_spec_fit.fun 

    # print(log_likelihood)

    # print(T,T_err)#, Index,Ecut,loglike)
    
    row = [
    t0.isot,t1.isot,
    T,T_err,
    EM,EM_err,
    Flux,Flux_err, 
    Index,Index_err,
    Ecut,Ecut_err
    ]
 
    df.loc[len(df)] = row

df.to_csv(outfile,index=False)

