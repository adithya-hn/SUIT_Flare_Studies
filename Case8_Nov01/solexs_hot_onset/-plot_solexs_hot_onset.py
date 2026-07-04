"""
@Author      : Adithya H N
@Created On  : 2026-04-23
@Last Updated: 2026-04-23
@Project     : SUIT pre-flare study
@Version     : 1.0

@Description
-----------
Brief description: Plotting hot onset phase with SoLEXS and SUIT
Modification:
30-06-2026: Modifing STIX plot to SoLEXS plot, ther is no non thermal fit, non thermal start is indicated by impulsive phase start time only.

"""
import matplotlib.pyplot as plt
from parfive import Downloader
import astropy.units as u
from astropy.time import Time
from datetime import datetime , timedelta
import matplotlib.dates as mdates
import pandas as pd
import numpy as np
import astropy.io.fits as fits
import scienceplots
plt.style.use('science')

#---------------Input parameters----------------
import glob
from pathlib import Path

# files = glob.glob("/path/to/folder/*abc.fits")
fit_file= next(Path(".").glob("*.fits"))#glob.glob("*SDD2_L1_puc_tb_fit_results_with_bkg_T_EM_case2.fits")
lc1_file =next(Path(".").glob("*SDD2_L1_2.03_11.99keV_30sec.lc"))#solexs 2-12 kev
# lc2_file= "AL1_SOLEXS_20241101_SDD2_L1_12.08_21.99keV_30sec.lc" #solexs 12-22 kev
lc3_file= 'LightCurve.csv' #helios above 22 kev
suit_dif=pd.read_csv('Diff_img_data_NB04.csv')
suit_int=pd.read_csv('c8_NB04_lc_data.csv')



Start_t = "2024-11-01T14:13:00"
End_t   = "2024-11-01T14:23:00"

thermal_range    = [datetime.fromisoformat("2024-11-01T14:13:00"),
                    datetime.fromisoformat("2024-11-01T14:19:20")]

nonthermal_range = [datetime.fromisoformat("2024-11-01T14:19:20"),
                    datetime.fromisoformat("2024-11-01T14:24:00")]

hot_onset_time   =  datetime.fromisoformat("2024-11-01T14:18:00")

lc1_data=fits.open(lc1_file)
# lc2_data=fits.open(lc2_file)
lc3_data=pd.read_csv(lc3_file)
t_em=fits.open(fit_file)

lc1_dt=np.array(lc1_data[1].data['TIME'],dtype='datetime64[s]')
lc1_rate=np.array(lc1_data[1].data['COUNTS'],dtype='float')

# lc2_dt=np.array(lc2_data[1].data['TIME'],dtype='datetime64[s]')
# lc2_rate=np.array(lc2_data[1].data['COUNTS'],dtype='float')

lc3_dt=np.array(lc3_data['TIME'],dtype='datetime64[s]')
lc3_rate=np.array(lc3_data['COUNTS'],dtype='float')

# print(t_em.info())
# df=pd.read_csv('stix_lightcurves.csv')
# band_labels=["4-8 keV", "9-12 keV", "13-22 keV", "22-30 keV"]
colors = ["steelblue", "darkorange", "mediumseagreen", "crimson"]

# t_em=pd.read_csv('7_Nov01_stix_hot_onset_phase.csv')


#-----------------------------------------------

temp=t_em[1].data['TEMPERATURE']
time_array=t_em[1].data['TIME']
em=t_em[1].data['EM']
t_er=t_em[1].data['TEMPERATURE_ERR']
em_er=t_em[1].data['EM_ERR']
t_dt=[]
# print(t_er)

for t in time_array:
    t_dt.append(datetime.fromtimestamp(t))#-timedelta(hours=5,minutes=30))
# print(t_dt)
excess_int=suit_dif['diff_count']
suit_dt=pd.to_datetime(suit_dif['Date'])
intensity=suit_int['Total']


fig, (ax,ax1,ax2) = plt.subplots(3, 1, sharex=True, figsize=(5,6), gridspec_kw={'hspace': 0}) 

ax.step(pd.to_datetime(lc1_dt), lc1_rate, where="mid", color='tab:blue', lw=1.2, label='SoLEXS 2-12 keV')
# ax.step(pd.to_datetime(lc2_dt), lc2_rate, where="mid", color='tab:orange', lw=1.2, label='SoLEXS 12-22 keV')
ax.step(pd.to_datetime(lc3_dt), lc3_rate, where="mid", color='tab:green', lw=1.2, label=r'HEL1OS $>$22 keV')

ax.axvline(datetime.fromisoformat("2024-11-01T14:18:00"),ls='--',lw=1,color='b',label='GOES flare start time',alpha=0.7)
ax.axvline(datetime.fromisoformat("2024-11-01T14:23:00"),ls='-',lw=1,color='b',label='Impuslive phase start time',alpha=0.7)
plot_range=[datetime.fromisoformat(Start_t),datetime.fromisoformat(End_t)+timedelta(minutes=4)]
ax.set_xlim(plot_range)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_major_locator(mdates.AutoDateLocator())
ax.set_ylabel(r"Count rate (s$^{-1}$)")
ax.set_title(f"Flare onset phase (SOL2024-11-01T14:23:00L)", fontsize=10)
ax.legend(fontsize=9, framealpha=0.6)
ax.grid(alpha=0.3, ls="--")
ax.set_ylim(5,1e6)
ax.set_yscale('log')
ax.axvspan(1.25, 1.55, facecolor='g', alpha=0.5)
fig.autofmt_xdate(rotation=30, ha="right")
ax1.set_xlim(plot_range)
ax11=ax1.twinx()
ax21=ax2.twinx()
xe=10/86400
ax1.errorbar(pd.to_datetime(t_dt[11:]),temp[11:],xerr=timedelta(seconds=10),yerr=t_er[11:],ls='dotted',color='tab:red',label='Temperature')
ax11.errorbar(pd.to_datetime(t_dt[10:]),em[10:],xerr=timedelta(seconds=10),yerr=em_er[10:],ls='dotted',color='tab:blue',label='Emission measure')

ax11.set_yscale('log')
ax1.set_ylabel('Temperature (MK)')
ax11.set_ylabel(r'EM ($\times 10^{46}$cm$^{-3}$)')

ax1.axvline(hot_onset_time, color="black", lw=1, ls="--", zorder=3, label="Hot onset start",alpha=0.8)
# # Shaded fit ranges
# ax1.axvspan(*thermal_range,    alpha=0.1, color="tomato",      zorder=0)#, label="Thermal fit range")
# ax1.axvspan(*nonthermal_range, alpha=0.1, color="mediumpurple", zorder=0)#, label="Non-thermal fit range")
# # Hot onset vertical line
# 

# ax1.annotate("Thermal fit", xy=(thermal_range[0], ax1.get_ylim()[1]),
#             xytext=(5, -12), textcoords="offset points",
#             fontsize=8, color="tomato", fontstyle="italic")

# ax1.annotate("Non-thermal fit", xy=(nonthermal_range[0], ax1.get_ylim()[1]),
#             xytext=(5, -12), textcoords="offset points",
#             fontsize=8, color="mediumpurple", fontstyle="italic")

# ax1.annotate("Hot onset", xy=(hot_onset_time, ax1.get_ylim()[1]),
#             xytext=(5, -12), textcoords="offset points",
#             fontsize=8, color="black", fontstyle="italic")

ln1, lbl1 = ax1.get_legend_handles_labels()
ln2, lbl2 = ax11.get_legend_handles_labels()
leg=ax11.legend(ln1 + ln2, lbl1 + lbl2, loc='lower left',frameon=True, framealpha=0.8, facecolor='white', edgecolor='none')
# leg.set_zorder(10)
#----------------------------------------------------

exposure= np.array(suit_int['Exposure'],dtype=float)
nb3_counts_er=np.sqrt(intensity*(exposure/1000))/(exposure/1000)
lc1_mean_er= np.sqrt(excess_int*(exposure/1000))/(exposure/1000)

ax21.errorbar(suit_dt,excess_int,yerr=lc1_mean_er,color='tab:blue',marker='o',markersize=3,label='Excess intensity')
ax21.set_ylabel(r'Excess intensity (DN/s)')
ax2.set_xlabel("Time (UT)")

ax2.errorbar(suit_dt,intensity/1e8,yerr=nb3_counts_er/1e8,color='k',marker='o',markersize=3,label='Mg II h counts')
ax2.set_ylabel(r'Mg II h counts ($\times 10^{8}$ DN/s)')

ln3, lbl3 = ax2.get_legend_handles_labels()
ln4, lbl4 = ax21.get_legend_handles_labels()

ax2.legend(ln3 + ln4, lbl3 + lbl4, loc='upper left')
ax21.set_yscale('log')
# ax2.set_yscale('log')
# ax1.legend(loc='lower right')
plt.tight_layout()
plt.savefig("c8_stix_onset_lc.pdf", dpi=300)
plt.show()

