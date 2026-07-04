"""
@Author      : Adithya H N
@Created On  : 2026-06-25
@Last Updated: 2026-06-25
@Project     : Pre-flare study using Aditya-L1
@Version     : 1.0

@Description
-----------
Brief description
"""
import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
from sunpy import timeseries as ts
from sunpy.data.sample import GOES_XRS_TIMESERIES
import numpy as np
from sunkit_instruments import goes_xrs
import matplotlib.dates as mdates
import datetime
from sunpy.timeseries.sources.goes import XRSTimeSeries
import astropy.units as u
import pandas as pd

# #Input
# #-----------------------
# start_time = '2024-11-01 02:00'
# end_time =   '2024-11-01 02:20'
# impulsive_phase_start = '2024-11-01 02:12'


# df=pd.read_csv('7_Nov01_stix_hot_onset_phase.csv')

# plt.figure(figsize=(8, 6))
# time_nums =  mdates.date2num(df.index.to_pydatetime())
# sc = plt.scatter( df['EM'],df['temperature'],c=time_nums, cmap='rainbow', edgecolor='k', s=10)

# # Add colorbar with time formatting
# cbar = plt.colorbar(sc, label='Time (UTC)')
# cbar.ax.yaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
# #plt.plot(gf["emission_measure"],gf["temperature"])
# plt.xscale("log")
# plt.ylabel('Temperature (in MK)')
# plt.xlabel('Emission measure ($cm^{-3}$)')
# plt.grid(True)

# plt.savefig('GOES_em_v_temp.png',dpi=300)
# plt.show()

"""
EM vs Temperature plot for HOPE (Hot Onset Pre-flare Event) phase analysis.
Reads a CSV with columns: time_start, time_end, T, T_er1, T_er2, EM, EM_er1, EM_er2
Error convention: er1 = lower uncertainty, er2 = upper uncertainty
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import matplotlib.ticker as ticker

# ── USER CONFIG ─────────────────────────────────────────────────────────────
CSV_FILE    = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/solexs_hot_onset/7_Nov01_stix_hot_onset_phase.csv"   # path to your CSV
OUTPUT_FILE = "em_vs_T_hope.png"    # set None to only display

LABEL_T  = r"Temperature $T$ (MK)"
LABEL_EM = r"Emission Measure EM ($10^{46}$ cm$^{-3}$)"
TITLE    = "HOPE Phase — EM–$T$ Evolution (2024 Nov 01)"

CMAP     = "plasma"          # colormap for time progression
ALPHA    = 0.85
MS       = 6                 # marker size
CAPSIZE  = 3
EW       = 1.2               # error bar linewidth
LW       = 0.6               # connecting line linewidth
# ────────────────────────────────────────────────────────────────────────────


def load_data(path):
    df = pd.read_csv(path)
    df["time_mid"] = pd.to_datetime(df["time_start"]) + (
        pd.to_datetime(df["time_end"]) - pd.to_datetime(df["time_start"])
    ) / 2
    return df


def plot_em_T(df):
    fig, ax = plt.subplots(figsize=(7, 5.5))

    # Numeric time for colormap
    t_num   = mdates.date2num(df["time_mid"])
    t_norm  = Normalize(vmin=t_num.min(), vmax=t_num.max())
    cmap    = plt.get_cmap(CMAP)
    colors  = cmap(t_norm(t_num))

    # Asymmetric errorbars: er1=lower, er2=upper
    T_errlo  = df["T_er1"].values
    T_errhi  = df["T_er2"].values
    EM_errlo = df["EM_er1"].values
    EM_errhi = df["EM_er2"].values

    T  = df["T"].values
    EM = df["EM"].values

    # Draw connecting line (faint, shows temporal evolution)
    ax.plot(EM, T, color="gray", lw=LW, zorder=1, alpha=0.4, ls="--")

    # Plot each point with its color
    for i in range(len(df)):
        ax.errorbar(
            EM[i], T[i],
            xerr=[[EM_errlo[i]], [EM_errhi[i]]],
            yerr=[[T_errlo[i]], [T_errhi[i]]],
            fmt="o",
            color=colors[i],
            ecolor=colors[i],
            elinewidth=EW,
            capsize=CAPSIZE,
            capthick=EW,
            ms=MS,
            alpha=ALPHA,
            zorder=2,
        )

    # Mark first and last points
    ax.annotate("Start", xy=(EM[0], T[0]),
                xytext=(4, 6), textcoords="offset points",
                fontsize=7, color=cmap(0.05))
    ax.annotate("End",   xy=(EM[-1], T[-1]),
                xytext=(4, -10), textcoords="offset points",
                fontsize=7, color=cmap(0.95))

    # Colorbar — time axis
    sm = ScalarMappable(cmap=cmap, norm=t_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, aspect=30)
    cbar.set_label("Time (UT)", fontsize=10)

    # Format colorbar ticks as HH:MM
    tick_locs = np.linspace(t_num.min(), t_num.max(), 5)
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels([
        mdates.num2date(t).strftime("%H:%M") for t in tick_locs
    ])
    cbar.ax.tick_params(labelsize=8)

    # Axes formatting
    ax.set_xlabel(LABEL_EM, fontsize=11)
    ax.set_ylabel(LABEL_T, fontsize=11)
    ax.set_title(TITLE, fontsize=11, pad=8)
    ax.tick_params(labelsize=9)

    # Log scale on EM (now X) if dynamic range is large
    em_range = EM.max() / max(EM.min(), 1e-10)
    if em_range > 10:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())

    ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.5)
    fig.tight_layout()

    if OUTPUT_FILE:
        fig.savefig(OUTPUT_FILE, dpi=200, bbox_inches="tight")
        print(f"Saved → {OUTPUT_FILE}")

    plt.show()


if __name__ == "__main__":
    df = load_data(CSV_FILE)
    print(f"Loaded {len(df)} time bins.")
    print(df[["time_mid", "T", "EM"]].to_string(index=False))
    plot_em_T(df)