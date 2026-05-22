
''''
Created on 4 apr 2026
@author: adithya-hn

Descrip: Fitting helios peaks with STIX with thermal and non-thermal model and comapre

'''


from astropy.time import Time, TimeDelta
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from parfive import Downloader
import astropy.units as u
from astropy.time import Time
from sunkit_spex.extern.stix import STIXLoader
from sunkit_spex.legacy.fitting.fitter import Fitter, load
from datetime import datetime, timedelta
from sys import path as sys_path
# sys_path.append('/home/adithya/Adithya_repos')
# from plots_styl import set_pub_style
# set_pub_style()

#---------------Input parameters----------------

spec_file="stx_spectrum_2411127496.fits"
srm_file="stx_srm_2411127496.fits"
# Start_t="2024-11-12T22:22:00" #no date before this time
# End_t="2024-11-12T23:12:00"

# spec_file="stx_spectrum_2411125825.fits"
# srm_file="stx_srm_2411125825.fits"
# Start_t="2024-11-12T23:08:00" 
# End_t="2024-11-13T00:20:30"
case='9_Nov13'

start_background_time = "2024-11-12T23:08:00"
end_background_time   = "2024-11-12T23:11:00"

helios_pks=np.genfromtxt('suit_helios_time_count_c09.csv',delimiter=',', names=True, dtype=None)
pks_dt=helios_pks['date_time']
fit_mode='both' #'th' 'both' 'th_nth'


columns_nth = ["time_start","time_end","T","T_er1","T_er2","EM","EM_er1","EM_er2","Flux","Flux_er1","Flux_er2","Index","Index_er1","Index_er2",
"Ecut","Ecut_er1","Ecut_er2","L"]
columns_th = ["time_start","time_end","T","T_er1","T_er2","EM","EM_er1","EM_er2","L"]
columns_th_nth = ["time_start","time_end","T1","T1_er1","T1_er2","EM1","EM1_er1","EM1_er2","L1","T","T_er1","T_er2","EM","EM_er1","EM_er2","Flux","Flux_er1","Flux_er2","Index","Index_er1","Index_er2",
"Ecut","Ecut_er1","Ecut_er2","L"]
if fit_mode=='th':
    df = pd.DataFrame(columns=columns_th)
    out_file='thermal_result.csv'
if fit_mode=='th_nth':
    df = pd.DataFrame(columns=columns_nth)
    out_file='non_thermal_result.csv'
if fit_mode=='both':
    df = pd.DataFrame(columns=columns_th_nth)
    out_file='2_model_result.csv'


#-----------------------------------------------
# try:
for i in range(len(pks_dt)):
    bin_size = TimeDelta(30, format='sec')
    print(pks_dt[i])
    t0=Time(pks_dt[i])-bin_size
    t1=t0+bin_size*2
    print(f"Fitting bin {i} : {t0.isot} - {t1.isot}")

    event_id=i
    time_profile_size = (9, 6)
    spec_plot_size = (14, 8)
    joint_spec_plot_size = (25, 10)
    tol = 1e-20
    spec_font_size = 18
    xlims, ylims = [3, 100], [1e-1, 1e6]
    plt.rcParams["font.size"] = spec_font_size
    stix_spec = STIXLoader(spectrum_file=spec_file, srm_file=srm_file)
    stix_spec.update_event_times(start=t0, end=t1)
    stix_spec.update_background_times(Time(start_background_time), Time(end_background_time))

    #-------------------------------------
    def fit_th(stix_spec):
        print('Fitting thermal part..')
        fitter = Fitter(stix_spec)
        fitter.model = "(f_vth)"
        fitter.loglikelihood = "cstat"
        fitter.show_params
        fitter.energy_fitting_range = [6,25]
        fitter.params["T1_spectrum1"] = {"Status":"free", "Value":6, "Bounds":(1, 30)}
        fitter.params["EM1_spectrum1"] = {"Status":"free", "Value":0.8, "Bounds":(0.1, 100)}
        stix_spec_fit = fitter.fit(tol=tol)
        # print(fitter.show_params)
        
        # print(stix_spec_fit)

        plt.figure(figsize=(12,8))
        axes, res_axes = fitter.plot()
        axes[0].set_xlim([5,30])
        axes[0].set_title(f'{t0.to_datetime().date()} {t0.to_datetime().time()} - {t1.to_datetime().time()}')
        for txt in axes[0].texts:
            label = txt.get_text().lower()
            if any(k in label for k in ["slope", "offset", "m=", "c="]):
                txt.set_visible(False)
        res_axes[0].set_ylim(-2,2)
        plt.savefig(f'{case}_stix_preflarePeak_{event_id}_th.png',dpi=300)
        plt.close()

        spec_mcmc = fitter.run_mcmc(number_of_walkers=4, steps_per_walker=1200,)
        fitter.burn_mcmc = 250
        # print(fitter.show_params)

        plt.figure(figsize=(12,8))
        axes, res_axes = fitter.plot()
        axes[0].set_xlim([5,30])
        axes[0].set_title(f'{t0.to_datetime().date()} {t0.to_datetime().time()} - {t1.to_datetime().time()}')
        for txt in axes[0].texts:
            label = txt.get_text().lower()
            if any(k in label for k in ["slope", "offset", "m=", "c="]):
                txt.set_visible(False)

        plt.savefig(f'{case}_peak_{event_id}_th_with_mcmc.png',dpi=300)
        plt.close()
        return fitter

    def fit_th_nth(stix_spec):
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
        fitter1.params["index1_spectrum1"] = {"Status": "free", "Value": 2, "Bounds": (1, 15)}
        fitter1.params["e_c1_spectrum1"] = {"Status": "free", "Value": 11, "Bounds": (10, 30)}
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
        fitter.params["index1_spectrum1"] = {"Status": "free", "Value": fitter1.params["index1_spectrum1"].Value, "Bounds": (1, 15)}
        fitter.params["e_c1_spectrum1"] = {"Status": "free", "Value": fitter1.params["e_c1_spectrum1"].Value, "Bounds": (10, 30)}
        stix_spec_fit = fitter.fit(tol=tol)


        plt.figure(figsize=(12,8))
        axes, res_axes = fitter.plot()
        axes[0].set_xlim([5,30])
        axes[0].set_title(f'{t0.to_datetime().date()} {t0.to_datetime().time()} - {t1.to_datetime().time()}')
        for txt in axes[0].texts:
            label = txt.get_text().lower()
            if any(k in label for k in ["slope", "offset", "m=", "c="]):
                txt.set_visible(False)
        res_axes[0].set_ylim(-2,2)
        plt.savefig(f'{case}_stix_preflarePeak_{event_id}_nth.png',dpi=300)
        plt.close()

        spec_mcmc = fitter.run_mcmc(number_of_walkers=10, steps_per_walker=1200,)
        fitter.burn_mcmc = 250

        plt.figure(figsize=(12,8))
        axes, res_axes = fitter.plot()
        for a in axes:
            axes[0].set_xlim([5,30])
        axes[0].set_title(f'{t0.to_datetime().date()} {t0.to_datetime().time()} - {t1.to_datetime().time()}')
        for txt in axes[0].texts:
            label = txt.get_text().lower()
            if any(k in label for k in ["slope", "offset", "m=", "c="]):
                txt.set_visible(False)
        plt.savefig(f'{case}_peak_{event_id}_th_nth_with_mcmc.png',dpi=300)
        plt.close()
        return fitter
    #-------------------------------------
    if fit_mode=='th':
        th_param=fit_th(stix_spec)
        # print(th_param.show_params)
        T,(T_er1,T_er2)=th_param.show_params[['Param']=='T1_spectrum1']['Value'],th_param.show_params[['Param']=='T1_spectrum1']['Error']
        EM,(EM_er1,EM_er2)=th_param.show_params[['Param']=='EM1_spectrum1']['Value'],th_param.show_params[['Param']=='EM1_spectrum1']['Error']
        L=th_param.show_params[['Param']=='Fit Stat.']['Value']
        row_th = [t0.isot,t1.isot,T,T_er1,T_er2,EM,EM_er1,EM_er2,L]
        df.loc[len(df)] = row_th
        
    if fit_mode=='th_nth':
        th_nth_param=fit_th_nth(stix_spec)
        print(th_nth_param.show_params)
        T,(T_er1,T_er2)=th_nth_param.show_params[0]['Value'],th_nth_param.show_params[0]['Error']
        EM,(EM_er1,EM_er2)=th_nth_param.show_params[1]['Value'],th_nth_param.show_params[1]['Error']
        Flux,(F_er1,F_er2)=th_nth_param.show_params[2]['Value'],th_nth_param.show_params[2]['Error']
        Index,(I_er1,I_er2)=th_nth_param.show_params[3]['Value'],th_nth_param.show_params[3]['Error']
        Ecut,(Ec_er1,Ec_er2)=th_nth_param.show_params[4]['Value'],th_nth_param.show_params[4]['Error']
        L=th_nth_param.show_params[5]['Value']
        row_nth = [t0.isot,t1.isot,T,T_er1,T_er2,EM,EM_er1,EM_er2,Flux,F_er1,F_er2, Index,I_er1,I_er2,Ecut,Ec_er1,Ec_er2,L]
        print(row_nth)
        df.loc[len(df)] = row_nth

    if fit_mode=='both':
        th_param=fit_th(stix_spec)
        T1,(T1_er1,T1_er2)=th_param.show_params[0]['Value'],th_param.show_params[0]['Error']
        EM1,(EM1_er1,EM1_er2)=th_param.show_params[1]['Value'],th_param.show_params[1]['Error']
        L1=th_param.show_params[2]['Value']

        th_nth_param=fit_th_nth(stix_spec)
        # print(th_nth_param.show_params)
        T,(T_er1,T_er2)=th_nth_param.show_params[0]['Value'],th_nth_param.show_params[0]['Error']
        EM,(EM_er1,EM_er2)=th_nth_param.show_params[1]['Value'],th_nth_param.show_params[1]['Error']
        Flux,(F_er1,F_er2)=th_nth_param.show_params[2]['Value'],th_nth_param.show_params[2]['Error']
        Index,(I_er1,I_er2)=th_nth_param.show_params[3]['Value'],th_nth_param.show_params[3]['Error']
        Ecut,(Ec_er1,Ec_er2)=th_nth_param.show_params[4]['Value'],th_nth_param.show_params[4]['Error']
        L=th_nth_param.show_params[5]['Value']
        row_comb=[t0.isot,t1.isot,T1,T1_er1,T1_er2,EM1,EM1_er1,EM1_er2,L1,T,T_er1,T_er2,EM,EM_er1,EM_er2,Flux,F_er1,F_er2, Index,I_er1,I_er2,Ecut,Ec_er1,Ec_er2,L]
        df.loc[len(df)] = row_comb
        # print(row_comb)
        df.to_csv(out_file,index=False)
df.to_csv(out_file,index=False)
# except:
#     print('could not exicute')
#     pass


