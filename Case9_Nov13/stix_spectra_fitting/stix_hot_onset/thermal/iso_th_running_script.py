



from astropy.time import Time, TimeDelta
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from parfive import Downloader
import astropy.units as u
from astropy.time import Time
from sunkit_spex.extern.stix import STIXLoader
from sunkit_spex.legacy.fitting.fitter import Fitter, load
from datetime import datetime, timedelta
import astropy.time as atime



flare_start = Time("2024-11-13T00:00:30")
flare_end   = Time("2024-11-13T00:22:00")
#---------------Input parameters----------------


spec_file="../stx_spectrum_2411125825.fits"
srm_file="../stx_srm_2411125825.fits"

case='9_Nov13'

start_background_time = "2024-11-12T23:08:00"
end_background_time   = "2024-11-12T23:11:00"

#-----------------------------------------------



bin_size_20 = TimeDelta(20, format='sec')
bin_size_30 = TimeDelta(30, format='sec')
tol = 1e-20

outfile = f"{case}_stix_timeResolved_nth_fit.csv"

columns = ["time_start","time_end",
"T","T_er1","T_er2",
"EM","EM_er1","EM_er2","L"]

df = pd.DataFrame(columns=columns)
log = {"20s": [], "30s": [], "skipped": []}

t = flare_start
i = 0

while t < flare_end:

    t0 = t
    t1 = t0 + bin_size_20

    print(f"\nFitting bin {i} : {t0.isot} - {t1.isot}")

    stix_spec = STIXLoader(spectrum_file=spec_file, srm_file=srm_file)

    bin_used = None
    try:
        stix_spec.update_event_times(start=t0, end=t1)
        bin_used = "20s"
        t = t0 + bin_size_20        # ← next bin starts after 20s
        # print(f"  ✓ Using 20s bin")

    except IndexError:
        t1 = t0 + bin_size_30
        # print(f"  ⚠ 20s failed → trying 30s: {t0.isot} - {t1.isot}")
        try:
            stix_spec.update_event_times(start=t0, end=t1)
            bin_used = "30s"
            t = t0 + bin_size_30    # ← next bin starts after 30s, no overlap
            print(f"  ✓ Using 30s bin")

        except IndexError:
            print(f"  ✗ Skipping bin {i} — no STIX data found")
            log["skipped"].append(t0.isot)
            t = t0 + bin_size_30    # ← still advance to avoid infinite loop
            i += 1
            continue

    log[bin_used].append(t0.isot)
    

    stix_spec.update_background_times(atime.Time(start_background_time), atime.Time(end_background_time))


    # print('Fitting thermal part..')
    if i==0:
        fitter = Fitter(stix_spec)
        fitter.model = "(f_vth)"
        fitter.loglikelihood = "cstat"
        fitter.show_params

        fitter.energy_fitting_range = [6,25]
        fitter.params["T1_spectrum1"] = {"Status":"free", "Value":6, "Bounds":(1, 30)}
        fitter.params["EM1_spectrum1"] = {"Status":"free", "Value":0.8, "Bounds":(0.1, 100)}
        stix_spec_fit = fitter.fit(tol=tol)
        fitter1 = Fitter(stix_spec)


    # print('Done with initial guessings. satrting final fitting..')
    else:
        fitter = Fitter(stix_spec)
        fitter.model = "(f_vth)"
        fitter.loglikelihood = "cstat"
        fitter.show_params

        fitter.energy_fitting_range = [6,25]
        fitter.params["T1_spectrum1"] = {"Status":"free", "Value":T, "Bounds":(1, 30)}
        fitter.params["EM1_spectrum1"] = {"Status":"free", "Value":EM, "Bounds":(0.1, 100)}

        stix_spec_fit = fitter.fit(tol=tol)

    
    plt.figure(figsize=(12,8))
    axes, res_axes = fitter.plot()
    title_str = f"{t0.isot}  –  {t1.isot}"
    axes[0].set_title(title_str)
    axes[0].set_xlim([5,30])
    plt.savefig(f'{case}_stix_onset_step_{i}.png',dpi=300)
    plt.close()

    spec_mcmc = fitter.run_mcmc(number_of_walkers=4, steps_per_walker=1200,)
    fitter.burn_mcmc = 250

    plt.figure(figsize=(12,8))
    axes, res_axes = fitter.plot()
    for a in axes:
        axes[0].set_xlim([5,30])
        axes[0].set_title(title_str)

    plt.savefig(f'{case}_step_{i}_with_mcmc.png',dpi=300)
    plt.close()
    
    tab = fitter.mcmc_table
    th_param=fitter
    T,(T_er1,T_er2)=th_param.show_params[0]['Value'],th_param.show_params[0]['Error']
    EM,(EM_er1,EM_er2)=th_param.show_params[1]['Value'],th_param.show_params[1]['Error']
    L=th_param.show_params[2]['Value']
    row_th = [t0.isot,t1.isot,T,T_er1,T_er2,EM,EM_er1,EM_er2,L]
    df.loc[len(df)] = row_th
    i += 1
df.to_csv(f'{outfile}',index=False)

