import numpy as np
from datetime import datetime, timedelta  
import matplotlib.pyplot as plt
from astropy.io import fits

# def running_mean(x, window):
#     """Simple running mean with NaN handling."""
#     out = np.full_like(x, np.nan, dtype=float)
#     for i in range(window-1, len(x)):
#         out[i] = np.nanmean(x[i-window+1:i+1])
#     return out

# import numpy as np

def weighted_mean(x, sx):
    """Error-weighted mean and its 1-sigma uncertainty."""
    w = 1.0 / sx**2
    mean = np.nansum(w * x) / np.nansum(w)
    sigma = np.sqrt(1.0 / np.nansum(w))
    return mean, sigma


def rolling_weighted_mean(x, sx, window):
    """
    Rolling error-weighted mean.
    Returns mean and uncertainty arrays.
    """
    n = len(x)
    xm = np.full(n, np.nan)
    sxm = np.full(n, np.nan)

    for i in range(window - 1, n):
        xi = x[i - window + 1:i + 1]
        sxi = sx[i - window + 1:i + 1]

        good = np.isfinite(xi) & np.isfinite(sxi)
        if np.sum(good) < 2:
            continue

        xm[i], sxm[i] = weighted_mean(xi[good], sxi[good])

    return xm, sxm


def hope_trigger_with_errors(time,
                             T, sT,
                             EM, sEM,
                             cadence_sec,
                             avg_sec=60,
                             diff_sec=180,
                             T_thr=5.0,
                             EM_norm_thr=0.005,
                             nsig=1,
                             persistence=1):
    
    T = np.asarray(T, float)
    EM = np.asarray(EM, float)
    sT = np.asarray(sT, float)
    sEM = np.asarray(sEM, float)

    # -----------------------------
    # 1. Error-weighted 60 s mean
    # -----------------------------
    avg_pts = max(2, int(avg_sec / cadence_sec))

    Tm, sTm = rolling_weighted_mean(T, sT, avg_pts)
    EMm, sEMm = rolling_weighted_mean(EM, sEM, avg_pts)

    # -----------------------------
    # 2. 180 s running difference
    # -----------------------------
    diff_pts = int(diff_sec / cadence_sec)

    Tdel = Tm - np.roll(Tm, diff_pts)
    EMdel = EMm - np.roll(EMm, diff_pts)

    sTdel = np.sqrt(sTm**2 + np.roll(sTm, diff_pts)**2)
    sEMdel = np.sqrt(sEMm**2 + np.roll(sEMm, diff_pts)**2)

    Tdel[:diff_pts] = np.nan
    EMdel[:diff_pts] = np.nan
    sTdel[:diff_pts] = np.nan
    sEMdel[:diff_pts] = np.nan

    # -----------------------------
    # 3. Positivity + significance
    # -----------------------------
    sig_T = Tdel - nsig * sTdel
    sig_EM = EMdel - nsig * sEMdel

    bad = sig_EM <= 0
    Tdel[bad] = 0.0
    EMdel[bad] = 0.0

    # -----------------------------
    # 4. Normalized EM
    # -----------------------------
    EMmax = np.nanmax(EMdel)
    EMnorm = EMdel / EMmax if EMmax > 0 else np.zeros_like(EMdel)

    # -----------------------------
    # 5. HOPE trigger condition
    # -----------------------------
    condition = (
        (sig_T >= T_thr) &
        (EMnorm >= EM_norm_thr)
    )

    trigger_index = None
    count = 0
    for i in range(len(condition)):
        if condition[i]:
            count += 1
            if count >= persistence:
                trigger_index = i
                break
        else:
            count = 0

    return trigger_index, {
        "T_avg": Tm,
        "sT_avg": sTm,
        "EM_avg": EMm,
        "sEM_avg": sEMm,
        "T_delta": Tdel,
        "sT_delta": sTdel,
        "EM_delta": EMdel,
        "sEM_delta": sEMdel,
        "EM_norm": EMnorm,
        "condition": condition
    }


def plot_hope_diagnostics(time, out, trigger_index=None,
                          T_thr=5.0, EM_norm_thr=0.005):
    """
    Diagnostic plots for HOPE trigger with uncertainties.
    """

    fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    # --------------------------------------------------
    # Panel 1: Averaged Temperature
    # --------------------------------------------------
    axs[0].plot(time, out["T_avg"], color="k", lw=1.5, label="T (60 s avg)")
    axs[0].fill_between(
        time,
        out["T_avg"] - out["sT_avg"],
        out["T_avg"] + out["sT_avg"],
        color="gray", alpha=0.3, label=r"$\pm 1\sigma$"
    )
    axs[0].set_ylabel("Temperature (MK)")
    axs[0].legend(loc="upper left")

    # --------------------------------------------------
    # Panel 2: Averaged Emission Measure
    # --------------------------------------------------
    axs[1].plot(time, out["EM_avg"], color="k", lw=1.5, label="EM (60 s avg)")
    axs[1].fill_between(
        time,
        out["EM_avg"] - out["sEM_avg"],
        out["EM_avg"] + out["sEM_avg"],
        color="gray", alpha=0.3
    )
    axs[1].set_ylabel(r"Emission Measure")
    axs[1].legend(loc="upper left")

    # --------------------------------------------------
    # Panel 3: Running-difference Temperature
    # --------------------------------------------------
    axs[2].plot(time, out["T_delta"], color="tab:red", lw=1.5,
                label=r"$T_\Delta$")
    axs[2].fill_between(
        time,
        out["T_delta"] - out["sT_delta"],
        out["T_delta"] + out["sT_delta"],
        color="tab:red", alpha=0.25
    )
    axs[2].axhline(T_thr, color="k", ls="--", lw=1,
                   label=f"{T_thr:.1f} MK threshold")
    axs[2].set_ylabel(r"$T_\Delta$ (MK)")
    axs[2].legend(loc="upper left")

    # --------------------------------------------------
    # Panel 4: Normalized EM running difference
    # --------------------------------------------------
    axs[3].plot(time, out["EM_norm"], color="tab:blue", lw=1.5,
                label=r"$EM_\Delta$ (normalized)")
    axs[3].axhline(EM_norm_thr, color="k", ls="--", lw=1,
                   label=f"{EM_norm_thr}")
    axs[3].set_ylabel(r"$EM_\Delta$ (norm.)")
    axs[3].set_xlabel("Time")
    axs[3].legend(loc="upper left")

    # --------------------------------------------------
    # HOPE trigger marker
    # --------------------------------------------------
    if trigger_index is not None:
        for ax in axs:
            ax.axvline(time[trigger_index], color="purple",
                       lw=2, ls=":", label="HOPE trigger")
        axs[0].legend(loc="upper right")

    plt.tight_layout()
    plt.show()



hdul=fits.open('/home/adithya/Downloads/AL1_SOLEXS_20241009_SDD2_L1_puc_tb_fit_results_T_EM.fits')
Solexs = hdul[1].data 
T=[float(tp) for tp in Solexs['TEMPERATURE']]
sT=[float(tpe) for tpe in Solexs['TEMPERATURE_ERR']]
EM=[float(em) for em in Solexs['EM']]
sEM=[float(eme) for eme in Solexs['EM_ERR']]
slt=Solexs[0]

Solexs=(np.loadtxt(f'csv_files/AL1_SOLEXS_20241113_SDD2_L1_puc_tb_fit_results_TEMP_EM.txt',skiprows=1,dtype='str')).transpose()
time_array4=[]

T=[float(tp) for tp in Solexs[1]]
sT=[float(tpe) for tpe in Solexs[2]]
EM=[float(em) for em in Solexs[3]]
sEM=[float(eme) for eme in Solexs[4]]
slt=Solexs[0]
time=np.array([datetime.fromtimestamp(float(tp)) for tp in Solexs[0]])
print(time)
# sltime=np.array([float(tp) for tp in Solexs[0]])
# base_time = datetime(2024, 11, 12, 23, 45, 0)  # Jun 1, 2025 07:00:00 UTC
# time_seconds = sltime-sltime[0]  # Convert string times to float seconds
# time = [base_time + timedelta(seconds=int(t)) for t in time_seconds]

idx, out = hope_trigger_with_errors(
    time,
    T, sT,
    EM, sEM,
    cadence_sec=120.0,
    nsig=2
)

if idx is not None:
    print("HOPE trigger at:", time[idx])
plot_hope_diagnostics(time, out, trigger_index=idx)
