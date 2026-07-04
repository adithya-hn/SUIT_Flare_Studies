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
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import matplotlib.ticker as ticker
import scienceplots
plt.style.use('science')

#---------------Input parameters----------------

impulsive_phase_start_time = "2024-10-09T01:31:40"
goes_flare_start_time = "2024-10-09T01:25:00"


thermal_range    = [datetime.fromisoformat("2024-10-09T01:25:30"),
                    datetime.fromisoformat("2024-10-09T01:30:40")]

nonthermal_range = [datetime.fromisoformat("2024-10-09T01:30:40"),
                    datetime.fromisoformat("2024-10-09T01:33:00")]

# hot_onset_time   =  [datetime.fromisoformat("2024-10-09T01:29:00"),
#                     datetime.fromisoformat("2024-10-09T01:30:40")]

fit_file='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/solexs_hot_onset/AL1_SOLEXS_20241009_SDD2_L1_puc_tb_fit_results_with_bkg_T_EM.fits'
lc1_file ="/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/solexs_hot_onset/AL1_SOLEXS_20241009_SDD2_L1_2.03_11.99keV_30sec.lc" #solexs 2-12 kev
# lc2_file= "AL1_SOLEXS_20241101_SDD2_L1_12.08_21.99keV_30sec.lc" #solexs 12-22 kev
lc3_file= 'LightCurve.csv' #helios above 22 kev

stix_t_em=pd.read_csv('c6_stix_hotonset_phase.csv')
suit_dif=pd.read_csv('Diff_img_data_NB04.csv')
suit_int=pd.read_csv('c6_NB04_lc_data.csv')
df=pd.read_csv('stix_lightcurves.csv')

lc1_data=fits.open(lc1_file)
# lc2_data=fits.open(lc2_file)
lc3_data=pd.read_csv(lc3_file)
t_em=fits.open(fit_file)

band_labels=["4-8 keV", "9-12 keV"]
colors = ["darkorange",  "crimson","steelblue", ]


lc1_dt=np.array(lc1_data[1].data['TIME'],dtype='datetime64[s]')
lc1_rate=np.array(lc1_data[1].data['COUNTS'],dtype='float')

lc3_dt=np.array(lc3_data['TIME'],dtype='datetime64[s]')
lc3_rate=np.array(lc3_data['COUNTS'],dtype='float')

#-----------------------------------------------

stix_temp=stix_t_em['T']
stix_em=stix_t_em['EM']
stix_time_array=stix_t_em['time_start']
stix_temp_er1=stix_t_em['T_er1']
stix_temp_er2=stix_t_em['T_er2'] 
stix_em_er1=stix_t_em['EM_er1']
stix_em_er2=stix_t_em['EM_er2']

temp=t_em[1].data['TEMPERATURE']
time_array=t_em[1].data['TIME']
em=t_em[1].data['EM']
t_er=t_em[1].data['TEMPERATURE_ERR']
em_er=t_em[1].data['EM_ERR']
t_dt=[]
for t in time_array:
    t_dt.append(datetime.fromtimestamp(t))#-timedelta(hours=5,minutes=30))
# print(t_dt)
excess_int=suit_dif['diff_count']
suit_dt=pd.to_datetime(suit_dif['Date'])
intensity=suit_int['Total']


fig, (ax,ax1,ax2) = plt.subplots(3, 1, sharex=True, figsize=(5,6), gridspec_kw={'hspace': 0}) 
for label, color in zip(band_labels, colors):
    ax.step(pd.to_datetime(df["time"]), df[label], where="mid", color=color, lw=1.2, label='STIX '+ label)
ax.step(pd.to_datetime(lc1_dt), lc1_rate, where="mid", color='tab:blue', lw=1.2, label='SoLEXS 2-12 keV')
# ax.step(pd.to_datetime(lc2_dt), lc2_rate, where="mid", color='tab:orange', lw=1.2, label='SoLEXS 12-22 keV')
ax.step(pd.to_datetime(lc3_dt), lc3_rate, where="mid", color='tab:green', lw=1.2, label=r'HEL1OS $>$22 keV')
ax.axvline(datetime.fromisoformat(goes_flare_start_time),ls='--',lw=1,color='b',label='GOES flare start time',alpha=0.7)
ax.axvline(datetime.fromisoformat(impulsive_phase_start_time) ,ls='-',lw=1,color='b',label='Impulsive phase start time',alpha=0.7)
plot_range=[datetime.fromisoformat(goes_flare_start_time)-timedelta(minutes=3),datetime.fromisoformat(impulsive_phase_start_time)+timedelta(minutes=1)]
ax.set_xlim(plot_range)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.xaxis.set_major_locator(mdates.AutoDateLocator())
ax.set_ylabel(r"Count rate (s$^{-1}$)")
ax.set_title(f"Flare onset phase (SOL2024-10-09T01:56)")
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
ax1.errorbar(pd.to_datetime(t_dt),temp,xerr=timedelta(seconds=10),yerr=t_er,ls='dotted',color='crimson',label='SoLEXS T')
ax11.errorbar(pd.to_datetime(t_dt),em,xerr=timedelta(seconds=10),yerr=em_er,ls='dotted',color='royalblue',label='SoLEXS EM')
ax1.errorbar(pd.to_datetime(stix_time_array),stix_temp,xerr=timedelta(seconds=10),yerr=[stix_temp_er2,stix_temp_er1],ls='dotted',color='forestgreen',label='STIX T')
ax11.errorbar(pd.to_datetime(stix_time_array),stix_em,xerr=timedelta(seconds=10),yerr=[stix_em_er2,stix_em_er1],ls='dotted',color='darkorange',label='STIX EM')
ax11.set_yscale('log')
ax1.set_ylabel('Temperature (MK)')
ax11.set_ylabel(r'EM ($\times 10^{46}$cm$^{-3}$)')

# Shaded fit ranges
ax1.axvspan(*thermal_range,    alpha=0.08, color="tomato",      zorder=0)#, label="Thermal fit range")
ax1.axvspan(*nonthermal_range, alpha=0.08, color="mediumpurple", zorder=0)#, label="Non-thermal fit range")
# ax1.axvspan(*hot_onset_time, alpha=0.2, color="#D9D9D9",      zorder=0)
# Hot onset vertical line

ax1.annotate("TH", xy=(thermal_range[0], ax1.get_ylim()[1]),
            xytext=(5, -12), textcoords="offset points",
            fontsize=8, color="tomato", fontstyle="italic")

ax1.annotate("NTH", xy=(nonthermal_range[0], ax1.get_ylim()[1]),
            xytext=(5, -12), textcoords="offset points",
            fontsize=8, color="mediumpurple", fontstyle="italic")

# ax1.annotate("Hot onset", xy=(hot_onset_time[0], ax1.get_ylim()[1]*0.8),rotation=90,
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
plt.savefig("c6_solex_onset_lc.pdf", dpi=300)
plt.close()

#--------------Emission measure Temperature plot-----------

OUTPUT_FILE = "em_vs_T_hope.pdf"    # set None to only display

LABEL_T  = r"Temperature (MK)"
LABEL_EM = r"Emission Measure ($10^{46}$ cm$^{-3}$)"
TITLE    = "EM-T Evolution (SOL2024-10-09T01:56)"

CMAP     = "plasma"          # colormap for time progression
ALPHA    = 0.85
MS       = 4                 # marker size
CAPSIZE  = 3
EW       = 1.2               # error bar linewidth
LW       = 0.6               # connecting line linewidth

EM=em
T=temp
# EM_errlo=em_er
# EM_errhi=em_er
# T_errlo=t_er
# T_errhi=t_er
# print(t_dt)
t_num   = mdates.date2num(t_dt)
t_norm  = Normalize(vmin=t_num.min(), vmax=t_num.max())
cmap1    = plt.get_cmap(CMAP)
cmap2    = plt.get_cmap('cividis')
colors1  = cmap1(t_norm(t_num))
colors2  = cmap2(t_norm(t_num))
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(EM, T, color="gray", lw=LW, zorder=1, alpha=0.4, ls="--")
# Plot each point with its color
for i in range(len(T)):
    ax.errorbar(EM[i],T[i] ,
        xerr=em_er[i],
        yerr=t_er[i],
        fmt="o",
        color=colors1[i],
        ecolor=colors1[i],
        elinewidth=EW,
        capsize=CAPSIZE,
        capthick=EW,
        ms=MS,
        alpha=ALPHA, zorder=2,)
# for i in range(4,len(stix_temp)):
#     ax.errorbar(stix_em[i],stix_temp[i],
#         xerr=[[stix_em_er1[i]], [stix_em_er2[i]]],
#         yerr=[[stix_temp_er1[i]], [stix_temp_er2[i]]],
#         fmt="o",
#         color=colors2[i],
#         ecolor=colors2[i],
#         elinewidth=EW,
#         capsize=CAPSIZE,
#         capthick=EW,
#         ms=MS,
#         alpha=ALPHA,
#         zorder=2,)

ax.annotate("Start", xy=(EM[0], T[0]),
                xytext=(4, 6), textcoords="offset points",
                fontsize=7, color=cmap1(0.05))
ax.annotate("End",   xy=(EM[-1], T[-1]),
                xytext=(4, -10), textcoords="offset points",
                fontsize=7, color=cmap1(0.95))
# Colorbar — time axis
sm = ScalarMappable(cmap=cmap1, norm=t_norm)
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
fig.savefig(OUTPUT_FILE, dpi=300, bbox_inches="tight")
print(f"Saved → {OUTPUT_FILE}")

plt.show()