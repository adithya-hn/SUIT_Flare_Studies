''''
Created on 4 apr 2026
@author: adithya-hn

Descrip: Fitting helios peaks with STIX with thermal and non-thermal model and compare

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
import matplotlib.dates as mdates
import matplotlib as mpl
import scienceplots
import re

plt.style.use('science')

mpl.rcParams.update({
    'axes.linewidth'     : 1.5,
    'axes.titleweight'   : 'bold',
    'xtick.major.size'   : 7,
    'xtick.minor.size'   : 4,
    'ytick.major.size'   : 7,
    'ytick.minor.size'   : 4,
    'xtick.major.width'  : 1.8,
    'xtick.minor.width'  : 1.2,
    'ytick.major.width'  : 1.8,
    'ytick.minor.width'  : 1.2,
    'lines.linewidth'    : 1.5,
    
})

plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

# ─────────────────────────────────────────────────────────────
# NUMBER FORMATTING UTILITIES
# ─────────────────────────────────────────────────────────────

def fmt_plain_asym(value, upper, lower):
    """
    Format value with asymmetric errors as plain decimal (no 10^x).

    Automatically picks decimal places from the smaller error magnitude.

    Examples:
        (15.2,  1.77, 1.02) → '$15.20^{+1.77}_{-1.02}$'
        (3.24,  0.43, 0.43) → '$3.24^{+0.43}_{-0.43}$'
        (0.113, 0.04, 0.05) → '$0.11^{+0.04}_{-0.05}$'
        (10.3,  0.40, 0.35) → '$10.30^{+0.40}_{-0.35}$'
    """
    import math
    try:
        smallest = min(abs(upper), abs(lower))
        if smallest > 0:
            # decimal places = position of first significant digit + 1
            decimals = max(2, -int(math.floor(math.log10(smallest))) + 1)
        else:
            decimals = 2
        decimals = min(decimals, 4)   # cap at 4 dp
    except Exception:
        decimals = 2

    return (
        rf'${value:.{decimals}f}'
        rf'^{{+{upper:.{decimals}f}}}'
        rf'_{{-{lower:.{decimals}f}}}$'
    )


def clear_and_annotate_params(fig, ax, fitter, param_config,
                              x=0.97, y_start=0.95, y_step=0.07,
                              fontsize=18):
    """
    1. Hide ALL SunKitSpex-generated text objects (figure-level AND axes-level).
    2. Redraw parameter annotations cleanly inside the spectrum axes
       using plain decimal formatting (no 10^x notation).

    Parameters
    ----------
    fig          : matplotlib Figure
    ax           : spectrum Axes
    fitter       : SunKitSpex Fitter (after fit/MCMC)
    param_config : list of dicts with keys:
                     param_key, label, color, (optional) decimals
    """
    # ── Step 1: hide everything SunKitSpex drew ───────────────
    # Figure-level texts (this is where SunKitSpex puts params)
    for txt in fig.texts:
        txt.set_visible(False)
    
    # Axes-level texts in ALL axes
    keep_fragments = ['bg', 'ln(l)', 'ln(', 'cstat', 'likelihood', 'minimiser']
    for a in fig.get_axes():
        for txt in a.texts:
            t = txt.get_text().strip()
            if t in ['bg', 'b'] or any(k in t for k in keep_fragments):
                        # keep the BG corner label
                if any(k in t for k in ['ln(l)', 'ln(', 'cstat', 'likelihood', 'minimiser']):
                    txt.set_position((0.03, 0.05))          # bottom-left
                    txt.set_transform(a.transAxes)
                    txt.set_ha('left')
                    txt.set_va('bottom')
                continue  
            txt.set_visible(False)

    # ── Step 2: get fitted values from fitter ────────────────
    # Try multiple access methods — show_params structure varies
    # between SunKitSpex versions and pre/post MCMC
    def _get_param_val_err(fitter, key):
        """
        Robustly extract (value, upper_err, lower_err) from fitter.
        Tries three methods in order:
          1. fitter.params[key].Value / .Error  (most direct)
          2. fitter.show_params DataFrame lookup (legacy)
          3. fitter._params dict fallback
        """
        # Method 1: direct fitter.params attribute (works pre/post MCMC)
        try:
            p     = fitter.params[key]
            val   = float(p.Value)
            err   = p.Error
            if err is None:
                upper = lower = 0.0
            elif hasattr(err, '__len__') and len(err) == 2:
                upper, lower = float(err[0]), float(err[1])
            else:
                upper = lower = float(err)
            return val, upper, lower
        except Exception:
            pass

        # Method 2: show_params DataFrame
        try:
            df = fitter.show_params
            # show_params may use 'Param' or index-based access
            # Try column 'Param' first
            if 'Param' in df.columns:
                # key might be 'T1_spectrum1' or just 'T1'
                short_key = key.replace('_spectrum1', '')
                for search_key in [key, short_key]:
                    mask = df['Param'].str.contains(search_key, na=False)
                    if mask.any():
                        row   = df[mask].iloc[0]
                        val   = float(row['Value'])
                        err   = row.get('Error', (0, 0))
                        if hasattr(err, '__len__') and len(err) == 2:
                            upper, lower = float(err[0]), float(err[1])
                        else:
                            upper = lower = float(err) if err else 0.0
                        return val, upper, lower
        except Exception:
            pass

        # Method 3: print show_params for debugging and return None
        try:
            print(f"[DEBUG] show_params columns: {fitter.show_params.columns.tolist()}")
            print(f"[DEBUG] show_params:{fitter.show_params}")
        except Exception as e:
            print(f"[DEBUG] Cannot read show_params: {e}")
        return None, None, None

    # ── Step 3: redraw cleanly inside spectrum axes ───────────
    for i, cfg in enumerate(param_config):
        key   = cfg['param_key']
        label = cfg['label']
        color = cfg.get('color', 'black')

        val, upper, lower = _get_param_val_err(fitter, key)

        if val is not None:
            val_str = fmt_plain_asym(val, upper, lower)
        else:
            val_str = '—'

        ax.text(
            x,
            y_start - i * y_step,
            f'{label}: {val_str}',
            transform = ax.transAxes,
            ha        = 'right',
            va        = 'top',
            fontsize  = fontsize,
            color     = color,
        )


# ─────────────────────────────────────────────────────────────
# PARAM CONFIG — label, color and show_params key per parameter
# ─────────────────────────────────────────────────────────────

PARAM_CONFIG_NTH = [
    {'param_key': 'T1_spectrum1',
     'label'   : 'Temperature (MK)',
     'color'   : '#1f77b4'},
    {'param_key': 'EM1_spectrum1',
     'label'   : r'Emission Measure ($10^{46}$ cm$^{-3}$)',
     'color'   : '#1f77b4'},
    {'param_key': 'total_eflux1_spectrum1',
     'label'   : r'Electron Flux ($10^{35}$ e$^{-1}$ s$^{-1}$)',
     'color'   : '#2ca02c'},
    {'param_key': 'index1_spectrum1',
     'label'   : r'Spectral Index ($\delta$)',
     'color'   : '#2ca02c'},
    {'param_key': 'e_c1_spectrum1',
     'label'   : r'Low Energy Cutoff $E_c$ (keV)',
     'color'   : '#2ca02c'},
]

PARAM_CONFIG_TH = [
    {'param_key': 'T1_spectrum1',
     'label'   : 'Temperature (MK)',
     'color'   : '#1f77b4'},
    {'param_key': 'EM1_spectrum1',
     'label'   : r'Emission Measure ($10^{46}$ cm$^{-3}$)',
     'color'   : '#1f77b4'},
]

# ─────────────────────────────────────────────────────────────
# INPUT PARAMETERS
# ─────────────────────────────────────────────────────────────

spec_file="../../stx_spectrum_2410315184.fits"
srm_file="../../stx_srm_2410315184.fits"
start_background_time = "2024-11-01T01:59:00"
end_background_time   = "2024-11-01T02:03:00"
case                  = '7_Nov01'
# helios_pks            = np.genfromtxt('helios_peaks.csv', delimiter=',', names=True, dtype=None)
# pks_dt                = helios_pks['date_time']
fit_mode              = 'th'   # 'th' | 'th_nth' | 'both'

columns_nth    = ["time_start","time_end","T","T_er1","T_er2","EM","EM_er1","EM_er2",
                  "Flux","Flux_er1","Flux_er2","Index","Index_er1","Index_er2",
                  "Ecut","Ecut_er1","Ecut_er2","L"]
columns_th     = ["time_start","time_end","T","T_er1","T_er2","EM","EM_er1","EM_er2","L"]
columns_th_nth = ["time_start","time_end","T1","T1_er1","T1_er2","EM1","EM1_er1","EM1_er2","L1",
                  "T","T_er1","T_er2","EM","EM_er1","EM_er2",
                  "Flux","Flux_er1","Flux_er2","Index","Index_er1","Index_er2",
                  "Ecut","Ecut_er1","Ecut_er2","L"]

if fit_mode == 'th':
    df = pd.DataFrame(columns=columns_th);     out_file = 'thermal_result.csv'
if fit_mode == 'th_nth':
    df = pd.DataFrame(columns=columns_nth);    out_file = 'non_thermal_result.csv'
if fit_mode == 'both':
    df = pd.DataFrame(columns=columns_th_nth); out_file = '2_model_result.csv'

# ─────────────────────────────────────────────────────────────
# MAIN LOOP
# ─────────────────────────────────────────────────────────────

# for i in range(len(pks_dt)):
i = 3
bin_size = TimeDelta(30, format='sec')
# print(pks_dt[i])
t0 = Time("2024-11-01T02:09:40")#Time(pks_dt[i]) - bin_size
t1 = Time("2024-11-01T02:10:00")#t0 + bin_size * 2
print(f"Fitting bin {i} : {t0.isot} - {t1.isot}")

event_id       = i
tol            = 1e-20
spec_font_size = 18
plt.rcParams["font.size"] = spec_font_size

stix_spec = STIXLoader(spectrum_file=spec_file, srm_file=srm_file)
stix_spec.update_event_times(start=t0, end=t1)
stix_spec.update_background_times(Time(start_background_time), Time(end_background_time))

# plot_range = [
#     datetime.fromisoformat(pks_dt[i]) - timedelta(minutes=7),
#     datetime.fromisoformat(pks_dt[i]) + timedelta(minutes=4),
# ]
stix_spec.lightcurve(energy_ranges=[[4,8],[9,12],[13,25],[22,30]])
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=3))
# plt.xlim(plot_range[0], plot_range[1])
plt.ylim(0.8, 1e6)
plt.xlabel('Time (UT)')
plt.savefig(f"stix_{case}_{i}_lc.png", dpi=300)
plt.close()

# ─────────────────────────────────────────────────────
def fit_th(stix_spec):
    print('Fitting thermal part..')
    fitter = Fitter(stix_spec)
    fitter.model          = "(f_vth)"
    fitter.loglikelihood  = "cstat"
    fitter.energy_fitting_range = [6, 25]
    fitter.params["T1_spectrum1"]  = {"Status":"free", "Value":6,   "Bounds":(1, 30)}
    fitter.params["EM1_spectrum1"] = {"Status":"free", "Value":0.8, "Bounds":(0.1, 100)}
    fitter.fit(tol=tol)
    fitter.run_mcmc(number_of_walkers=4, steps_per_walker=1200)
    fitter.burn_mcmc = 250

    plt.figure(figsize=(8, 9))
    axes, res_axes = fitter.plot()
    axes[0].set_xlim([6, 25])
    axes[0].set_title(
        f'{t0.to_datetime().date()} '
        f'{t0.to_datetime().time()} - {t1.to_datetime().time()}'
    )
    # ── Apply renaming + 10^x formatting ──
    clear_and_annotate_params(plt.gcf(), axes[0], fitter, PARAM_CONFIG_TH)

    plt.savefig(f'{case}_peak_{event_id}_th_with_mcmc.pdf', dpi=300)
    plt.close()
    return fitter

# ─────────────────────────────────────────────────────
def fit_th_nth(stix_spec):
    # Step 1: thermal-only fit for initial guess
    print('Initial thermal fit (6-9 keV)..')
    fitter = Fitter(stix_spec)
    fitter.model         = "(f_vth)"
    fitter.loglikelihood = "cstat"
    fitter.energy_fitting_range = [6, 9]
    fitter.params["T1_spectrum1"]  = {"Status":"free", "Value":6,   "Bounds":(1, 30)}
    fitter.params["EM1_spectrum1"] = {"Status":"free", "Value":0.8, "Bounds":(0.1, 100)}
    fitter.fit(tol=tol)

    # Step 2: fix thermal, fit nonthermal (9-25 keV)
    print('Initial nonthermal fit (9-25 keV)..')
    fitter1 = Fitter(stix_spec)
    fitter1.model         = "(f_vth+thick_fn)"
    fitter1.loglikelihood = "cstat"
    fitter1.energy_fitting_range = [9, 25]
    fitter1.params["T1_spectrum1"]              = {"Status":"fix",  "Value":fitter.params["T1_spectrum1"].Value,  "Bounds":(1, 30)}
    fitter1.params["EM1_spectrum1"]             = {"Status":"fix",  "Value":fitter.params["EM1_spectrum1"].Value, "Bounds":(0.1, 100)}
    fitter1.params["total_eflux1_spectrum1"]    = {"Status":"free", "Value":2,   "Bounds":(1e-1, 100)}
    fitter1.params["index1_spectrum1"]          = {"Status":"free", "Value":2,   "Bounds":(1, 15)}
    fitter1.params["e_c1_spectrum1"]            = {"Status":"free", "Value":11,  "Bounds":(10, 30)}
    fitter1.fit(tol=tol)

    # Step 3: full joint fit (6-25 keV)
    print('Final joint fit (6-25 keV)..')
    fitter = Fitter(stix_spec)
    fitter.model         = "(f_vth+thick_fn)"
    fitter.loglikelihood = "cstat"
    fitter.energy_fitting_range = [6, 25]
    fitter.params["T1_spectrum1"]           = {"Status":"free", "Value":fitter1.params["T1_spectrum1"].Value,           "Bounds":(1, 30)}
    fitter.params["EM1_spectrum1"]          = {"Status":"free", "Value":fitter1.params["EM1_spectrum1"].Value,          "Bounds":(0.1, 100)}
    fitter.params["total_eflux1_spectrum1"] = {"Status":"free", "Value":fitter1.params["total_eflux1_spectrum1"].Value, "Bounds":(1e-1, 100)}
    fitter.params["index1_spectrum1"]       = {"Status":"free", "Value":fitter1.params["index1_spectrum1"].Value,       "Bounds":(1, 15)}
    fitter.params["e_c1_spectrum1"]         = {"Status":"free", "Value":fitter1.params["e_c1_spectrum1"].Value,         "Bounds":(10, 30)}
    fitter.fit(tol=tol)

    # ── Pre-MCMC plot ──────────────────────────────────
    plt.figure(figsize=(12, 9))
    axes, res_axes = fitter.plot()
    axes[0].set_xlim([6, 25])
    axes[0].set_ylim(1e-2, 1e3)
    res_axes[0].set_ylim(-2, 2)
    axes[0].set_title(
        f'{t0.to_datetime().date()} '
        f'{t0.to_datetime().time()} - {t1.to_datetime().time()}'
    )
    clear_and_annotate_params(plt.gcf(), axes[0], fitter, PARAM_CONFIG_NTH)
    plt.savefig(f'{case}_stix_preflarePeak_{event_id}_nth.png', dpi=300)
    plt.close()

    # ── MCMC ──────────────────────────────────────────
    fitter.run_mcmc(number_of_walkers=10, steps_per_walker=1200)
    fitter.burn_mcmc = 250

    # ── Post-MCMC plot ─────────────────────────────────
    plt.figure(figsize=(12, 10))
    axes, res_axes = fitter.plot()
    axes[0].set_xlim([6, 25])
    axes[0].set_ylim(5e-4, 1e3)
    axes[0].set_title(
        f'{t0.to_datetime().date()} '
        f'{t0.to_datetime().time()} - {t1.to_datetime().time()}'
    )
    clear_and_annotate_params(plt.gcf(), axes[0], fitter, PARAM_CONFIG_NTH)
    plt.savefig(f'{case}_peak_{event_id}_th_nth_with_mcmc.png', dpi=300)
    plt.close()
    return fitter

# ─────────────────────────────────────────────────────
if fit_mode == 'th':
    th_param = fit_th(stix_spec)
    T,  (T_er1,  T_er2)  = th_param.show_params[0]['Value'], th_param.show_params[0]['Error']
    EM, (EM_er1, EM_er2) = th_param.show_params[1]['Value'], th_param.show_params[1]['Error']
    L                    = th_param.show_params[2]['Value']
    df.loc[len(df)]      = [t0.isot, t1.isot, T, T_er1, T_er2, EM, EM_er1, EM_er2, L]

if fit_mode == 'th_nth':
    th_nth_param         = fit_th_nth(stix_spec)
    T,    (T_er1,  T_er2)  = th_nth_param.show_params[0]['Value'], th_nth_param.show_params[0]['Error']
    EM,   (EM_er1, EM_er2) = th_nth_param.show_params[1]['Value'], th_nth_param.show_params[1]['Error']
    Flux, (F_er1,  F_er2)  = th_nth_param.show_params[2]['Value'], th_nth_param.show_params[2]['Error']
    Index,(I_er1,  I_er2)  = th_nth_param.show_params[3]['Value'], th_nth_param.show_params[3]['Error']
    Ecut, (Ec_er1, Ec_er2) = th_nth_param.show_params[4]['Value'], th_nth_param.show_params[4]['Error']
    L                      = th_nth_param.show_params[5]['Value']
    df.loc[len(df)] = [t0.isot, t1.isot, T, T_er1, T_er2, EM, EM_er1, EM_er2,
                        Flux, F_er1, F_er2, Index, I_er1, I_er2, Ecut, Ec_er1, Ec_er2, L]

if fit_mode == 'both':
    th_param             = fit_th(stix_spec)
    T1,  (T1_er1,  T1_er2)  = th_param.show_params[0]['Value'], th_param.show_params[0]['Error']
    EM1, (EM1_er1, EM1_er2) = th_param.show_params[1]['Value'], th_param.show_params[1]['Error']
    L1                       = th_param.show_params[2]['Value']

    th_nth_param         = fit_th_nth(stix_spec)
    T,    (T_er1,  T_er2)  = th_nth_param.show_params[0]['Value'], th_nth_param.show_params[0]['Error']
    EM,   (EM_er1, EM_er2) = th_nth_param.show_params[1]['Value'], th_nth_param.show_params[1]['Error']
    Flux, (F_er1,  F_er2)  = th_nth_param.show_params[2]['Value'], th_nth_param.show_params[2]['Error']
    Index,(I_er1,  I_er2)  = th_nth_param.show_params[3]['Value'], th_nth_param.show_params[3]['Error']
    Ecut, (Ec_er1, Ec_er2) = th_nth_param.show_params[4]['Value'], th_nth_param.show_params[4]['Error']
    L                      = th_nth_param.show_params[5]['Value']
    df.loc[len(df)] = [t0.isot, t1.isot,
                        T1, T1_er1, T1_er2, EM1, EM1_er1, EM1_er2, L1,
                        T,  T_er1,  T_er2,  EM,  EM_er1,  EM_er2,
                        Flux, F_er1, F_er2, Index, I_er1, I_er2, Ecut, Ec_er1, Ec_er2, L]
    df.to_csv(out_file, index=False)

df.to_csv(out_file, index=False)