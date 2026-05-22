"""
@Author      : Adithya H N
@Created On  : 2026-04-23
@Last Updated: 2026-04-23
@Project     : SUIT pre-flare study
@Version     : 1.0

@Description
-----------
Brief description: Plotting hot onset phase with STIX and SUIT

plot energy wise light curve from stix in panel1, indicate the thermal and nothermal fitting range
plot temp and em in same plot
third row with excess intenisity and light curve of the time range

"""
import matplotlib.pyplot as plt
from parfive import Downloader
import astropy.units as u
from astropy.time import Time
from datetime import datetime , timedelta
import matplotlib.dates as mdates
import pandas as pd
import numpy as np
import scienceplots
plt.style.use('science')


#---------------Input parameters----------------

Start_t = "2024-11-01T02:00:00"
End_t   = "2024-11-01T02:15:00"

thermal_range    = [datetime.fromisoformat("2024-11-01T02:00:00"),
                    datetime.fromisoformat("2024-11-01T02:10:20")]

nonthermal_range = [datetime.fromisoformat("2024-11-01T02:10:20"),
                    datetime.fromisoformat("2024-11-01T02:16:00")]

hot_onset_time   =  datetime.fromisoformat("2024-11-01T02:07:40")

df=pd.read_csv('stix_lightcurves.csv')
band_labels=["4-8 keV", "9-12 keV", "13-22 keV", "22-30 keV"]
colors = ["steelblue", "darkorange", "mediumseagreen", "crimson"]

t_em=pd.read_csv('7_Nov01_stix_hot_onset_phase.csv')
suit_dif=pd.read_csv('Diff_img_data_NB04.csv')
suit_int=pd.read_csv('c7_NB04_lc_data.csv')

#-----------------------------------------------

temp=t_em['T']
t_dt=t_em['time_start']
em=t_em['EM']
t_er_l=t_em['T_er1']
t_er_h=t_em['T_er2']
em_er_l=t_em['EM_er1']
em_er_h=t_em['EM_er2']


excess_int=suit_dif['diff_count']
suit_dt=pd.to_datetime(suit_dif['Date'])
intensity=suit_int['Total']


fig, (ax, ax1,ax2) = plt.subplots(3, 1, sharex=True, figsize=(5,6), gridspec_kw={'hspace': 0}) 


for label, color in zip(band_labels, colors):
    ax.step(pd.to_datetime(df["time"]), df[label], where="mid", color=color, lw=1.2, label=label)

ax.axvline(datetime.fromisoformat("2024-11-01T02:05:00"),ls='--',lw=1,color='b',label='GOES flare start time',alpha=0.7)
ax.axvline(datetime.fromisoformat("2024-11-01T02:13:25"),ls='-',lw=1,color='b',label='Impuslive phase start time',alpha=0.7)
plot_range=[datetime.fromisoformat(Start_t),datetime.fromisoformat(End_t)+timedelta(minutes=1)]
ax.set_xlim(plot_range)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_major_locator(mdates.AutoDateLocator())
ax.set_ylabel(r"Counts $s^{-1}$")
ax.set_title(f"Flare onset phase (SOL2024-11-01T02:16)")
ax.legend(fontsize=9, framealpha=0.6)
ax.grid(alpha=0.3, ls="--")
ax.set_ylim(5,1e6)
ax.set_yscale('log')

# fig.autofmt_xdate(rotation=30, ha="right")

ax11=ax1.twinx()
ax21=ax2.twinx()
xe=10/86400
ax1.errorbar(pd.to_datetime(t_dt),temp,xerr=timedelta(seconds=10),yerr=[t_er_h,t_er_l],ls='dotted',color='tab:red',label='Temperature')
ax11.errorbar(pd.to_datetime(t_dt),em,xerr=timedelta(seconds=10),yerr=[em_er_h,em_er_l],ls='dotted',color='tab:blue',label='Emission measure')
ax.axvspan(1.25, 1.55, facecolor='g', alpha=0.5)
ax11.set_yscale('log')
ax1.set_ylabel('Temperature (MK)')
ax11.set_ylabel(r'EM ($\times 10^{46}~cm^{-3}$)')

# Shaded fit ranges
ax1.axvspan(*thermal_range,    alpha=0.1, color="tomato",      zorder=0)#, label="Thermal fit range")
ax1.axvspan(*nonthermal_range, alpha=0.1, color="mediumpurple", zorder=0)#, label="Non-thermal fit range")
# Hot onset vertical line
ax1.axvline(hot_onset_time, color="black", lw=1, ls="--", zorder=3, label="Hot onset start",alpha=0.8)

ax1.annotate("Thermal fit", xy=(thermal_range[0], ax1.get_ylim()[1]),
            xytext=(5, -12), textcoords="offset points",
            fontsize=8, color="tomato", fontstyle="italic")

ax1.annotate("Non-thermal fit", xy=(nonthermal_range[0], ax1.get_ylim()[1]),
            xytext=(5, -12), textcoords="offset points",
            fontsize=8, color="mediumpurple", fontstyle="italic")

ax1.annotate("Hot onset", xy=(hot_onset_time, ax1.get_ylim()[1]),
            xytext=(5, -12), textcoords="offset points",
            fontsize=8, color="black", fontstyle="italic")
ln1, lbl1 = ax1.get_legend_handles_labels()
ln2, lbl2 = ax11.get_legend_handles_labels()
ax1.legend(ln1 + ln2, lbl1 + lbl2, loc='lower left',frameon=True, framealpha=0.8, facecolor='white', edgecolor='none')

#----------------------------------------------------

exposure= np.array(suit_int['Exposure'],dtype=float)
nb3_counts_er=np.sqrt(intensity*(exposure/1000))/(exposure/1000)
lc1_mean_er= np.sqrt(excess_int*(exposure/1000))/(exposure/1000)

ax21.errorbar(suit_dt,excess_int,yerr=lc1_mean_er,color='tab:blue',marker='o',markersize=3,label='Excess intensity')
ax21.set_ylabel(r'Excess intensity DN/s')
ax21.set_xlabel("Time (UT)")

ax2.errorbar(suit_dt,intensity/1e8,yerr=nb3_counts_er/1e8,color='k',marker='o',markersize=3,label='Mg II h counts')
ax2.set_ylabel(r'Mg II h counts $\times 10^{8}$ DN/s ')

ln3, lbl3 = ax2.get_legend_handles_labels()
ln4, lbl4 = ax21.get_legend_handles_labels()

ax2.legend(ln3 + ln4, lbl3 + lbl4, loc='upper left')
ax21.set_yscale('log')
# ax2.set_yscale('log')
# ax1.legend(loc='lower right')
plt.tight_layout()
plt.savefig("c7_stix_onset_lc.pdf", dpi=300)
plt.show()

