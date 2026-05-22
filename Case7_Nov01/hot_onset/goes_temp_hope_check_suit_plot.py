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

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()
import logging
logging.getLogger("matplotlib.texmanager").setLevel(logging.ERROR)
logging.getLogger("sunpy").setLevel(logging.ERROR)
from matplotlib.ticker import MultipleLocator


plt.rcParams['legend.fontsize'] = 19
ml = MultipleLocator(5)
# 1. Define the time range
start_time = '2024-11-01 02:00'
end_time =   '2024-11-01 02:20'
impulsive_phase_start = '2024-11-01 02:12'

result = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"),a.goes.SatelliteNumber(16),a.Resolution("avg1m"))

# 3. Download the data
downloaded_files = Fido.fetch(result)

# 4. Load as TimeSeries
goes_ts = ts.TimeSeries(downloaded_files,source='XRS')
goes_df = goes_ts.to_dataframe()
df_rd = goes_df.diff(periods=3)

df_rd["xrsa"] = df_rd["xrsa"].values
df_rd["xrsb"] = df_rd["xrsb"].values

diff_ts = XRSTimeSeries(df_rd, meta=goes_ts.meta)
diff_ts.units = goes_ts.units


goes_flare = goes_ts.truncate(start_time, end_time)
goes_temp_em =goes_xrs.calculate_temperature_em(goes_flare) #goes_xrs.calculate_temperature_em(goes_flare)
flare_df=goes_flare.to_dataframe()
goes_diff_flare = diff_ts.truncate(start_time, end_time)
goes_diff_temp_em =goes_xrs.calculate_temperature_em(goes_diff_flare)
goes_flare_df = goes_diff_flare.to_dataframe()
neg_mask = (goes_flare_df["xrsa"].values <= 0) | (goes_flare_df["xrsb"].values <= 0)
df1 = goes_temp_em.to_dataframe()
df2 = goes_diff_temp_em.to_dataframe()
msk2= (df2["temperature"].values > 100)
df2.loc[neg_mask, "temperature"] = 0.0
df2.loc[msk2, "temperature"] = 0.0
df2.loc[neg_mask, "emission_measure"] = 0.0

condition = (df2["temperature"] >= 5.0) & (df2["emission_measure"] >= 5e46)
if condition.any():
    trigger_time = df2.index[condition.argmax()]
else:
    trigger_time = None
df2["HOPE_trigger"] = False
if trigger_time is not None:
    df2.loc[df2.index == trigger_time, "HOPE_trigger"] = True
    print("HOPE trigger at:", trigger_time)

# ts_modified = XRSTimeSeries(df2,meta=goes_ts.meta,units=goes_ts.units) 
# df2 = ts_modified.to_dataframe()

suit_diff=pd.read_csv('Diff_img_data_NB04.csv')
suit_lc=pd.read_csv('c7_NB04_lc_data.csv')
dt   =pd.to_datetime(suit_diff['Date'], errors='coerce')
lc_dt=pd.to_datetime(suit_lc['Time'], errors='coerce')
# print(lc_dt)
lc1_tot = np.array(suit_lc['Total'],dtype=float)
lc1_mean_er= np.sqrt(lc1_tot*(suit_lc['Exposure']/1000))/(suit_lc['Exposure']/1000)

# fig, (ax1, ax2, ax3,ax6) = plt.subplots(4,figsize=(16, 16), sharex=True)
fig, (ax1, ax2, ax3,ax6) = plt.subplots(4, 1, sharex=True, figsize=(12,16), gridspec_kw={'hspace': 0})  # no vertical spacing 
#fig.suptitle("X5.8 Class Flare on 2024-05-11 01:10:00 UTC")
ax4 = ax2.twinx()
ax5 = ax3.twinx()
# ax7 = ax6.twinx()
ax11= ax1.twinx()

ax11.plot(goes_flare_df.index, goes_flare_df["xrsa"],markersize=7,marker='o', color="r", label=r"$\Delta$ XRS-A")
ax11.plot(goes_flare_df.index, goes_flare_df["xrsb"],markersize=7,marker='s', color="orange", label=r"$\Delta$ XRS-B")
ax1.plot(flare_df.index, flare_df["xrsa"], color="k",markersize=7,marker='o',linestyle='--', label="XRS-A")
ax1.plot(flare_df.index, flare_df["xrsb"],markersize=7,marker='s',linestyle='--', color="gray", label="XRS-B")
ax1.set_ylabel("XRS Flux (W/m$^2$)")
ax11.set_ylabel(r"$\Delta$ XRS Flux (MK) (W/m$^2$)",color='red')
ax6.set_ylabel(r'Excess intensity (DN/s) ',color='k')

ax1.set_yscale("log")
ax11.set_yscale("log")
ax1.legend(loc="upper right",fontsize=2)
ax1.minorticks_on()
# ax1.xaxis.set_minor_locator(MultipleLocator(10))
# ax1.grid(True)
ax2.plot(df1.index, df1["temperature"], color="k",markersize=7,marker='o',linestyle='--', label="Temperature")
ax4.plot(df2.index, df2["temperature"], color="r",markersize=7,marker='o', label=r"$\Delta$ Temperature")
ax4.axhline(5, color="g", ls="--", lw=2)
ax2.set_ylabel("Temperature (MK)")
ax2.legend(loc="upper right")
# ax2.tick_params(axis='x', which='both',  pad=-28)   # shift labels inward
ax2.tick_params(labelbottom=False)
if trigger_time is not None:
    ax2.axvline(trigger_time, color="b", ls="--", lw=2)
    ax2.axvline(impulsive_phase_start, color="b", ls="-", lw=1)
    ax3.axvline(trigger_time, color="b", ls="--", lw=2,label="HOPE trigger")
    ax3.axvline(impulsive_phase_start, color="b", ls="-", lw=1,label="Impulsive phase start time")
    ax6.axvline(trigger_time, color="b", ls="--", lw=2)
    ax6.axvline(impulsive_phase_start, color="b", ls="-", lw=1)

ax3.plot(df1.index, df1["emission_measure"]/1e48,markersize=7,marker='o',linestyle='--', color="k",label="Emission Measure")
ax5.plot(df2.index, df2["emission_measure"]/1e48,markersize=7,marker='o', color="r",label=r"$\Delta$ Emission Measure")
ax6.errorbar(dt,suit_diff['diff_count'].values/1e7, yerr=suit_diff['Diff_error'].values/1e7,marker='o',markersize=7,color='k',label=r'SUIT Mg II h Excess intensity')
# ax7.errorbar(lc_dt, lc1_tot/10e8,yerr=lc1_mean_er/10e5,fmt='k', marker="o",capsize=2,markersize=2,linewidth=0.5, label=r"SUIT Mg II h (errors multiplied by $ 10^3$)")
ax3.set_ylabel("EM ($10^{48}$ cm$^{-3}$)")
ax5.axhline(5e-2, color="g", ls="--", lw=2)
ax3.legend(loc="upper right")
# ax6.legend(loc="upper right")
#ax3.set_yscale("log")
ax5.set_yscale("log")
ax6.legend()

#ax2.set_ylim(3,20)
#ax3.set_ylim(1e-3,1e1)
# ax2.grid(True)
# ax3.grid(True)
ax4.set_ylabel(r"$\Delta$ Temperature (MK)",color='red')
ax5.set_ylabel(r"$\Delta$ EM($10^{48}$ cm$^{-3}$)",color='red')
# ax3.axvline(datetime.datetime(2024, 5, 11, 1, 13), color="b", ls="-", lw=1)
# ax2.axvline(datetime.datetime(2024, 5, 11, 1, 13), color="b", ls="-", lw=1)

lines1, labels1 = ax2.get_legend_handles_labels()
lines2, labels2 = ax4.get_legend_handles_labels()
ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper left")


lines3, labels3 = ax3.get_legend_handles_labels()
lines4, labels4 = ax5.get_legend_handles_labels()
ax3.legend(lines3 + lines4, labels3 + labels4, loc="upper left")

lines11, labels11 = ax11.get_legend_handles_labels()
lines12, labels12 = ax1.get_legend_handles_labels()
ax1.legend(lines12 + lines11, labels12 + labels11, loc="upper left")


ax4.set_ylim(-0.7,16)
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
ax6.set_xlabel(f"Time (UT)")
panel_labels = ['a)', 'b)', 'c)', 'd)']

ax1.text(0.03, 0.21,'a)',transform=ax1.transAxes, fontsize=24, fontweight='bold',va='top', ha='left')
ax2.text(0.03, 0.21,'b)',transform=ax2.transAxes, fontsize=24, fontweight='bold',va='top', ha='left')
ax3.text(0.03, 0.21,'c)',transform=ax3.transAxes, fontsize=24, fontweight='bold',va='top', ha='left')
ax6.text(0.03, 0.21,'d)',transform=ax6.transAxes, fontsize=24, fontweight='bold',va='top', ha='left')

plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig('c7_GOES_temp_em1.pdf',dpi=300)
plt.close()

df1.index.name = "date_time"
df1.to_csv('goes_temp_em.csv')
df2.index.name = "date_time"
df2.to_csv('goes_diff_temp_em.csv')
###########################################


# plt.figure(figsize=(8, 6))
# time_nums =  mdates.date2num(df1.index.to_pydatetime())
# sc = plt.scatter( df1['emission_measure'],df1['temperature'],c=time_nums, cmap='rainbow', edgecolor='k', s=10)

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


# fig,ax1=plt.subplots(1,1, figsize=(8,4))
# plt.rcParams["font.size"]=12
# plt.rcParams["axes.labelsize"]=12
# plt.rcParams["figure.titlesize"]=12

# ax2=ax1.twinx()

# ax1.errorbar(dt,suit_diff['diff_count'].values/1e7, yerr=suit_diff['Diff_error'].values/1e7,marker='o',markersize=4,linestyle='None',color='k',label=r'SUIT Mg II h Excess intensity')
# ax2.plot(df2.index, df2["temperature"], color="r",markersize=4,marker='s',linestyle='None', label=r"$\Delta$ Temperature")
# plt.xticks(rotation=45)
# ax1.set_ylabel(r'Mg II h counts ($\times 10^7$ DN/s) ',fontsize=12)
# ax2.set_ylabel(r"$\Delta$ Temperature (MK)",fontsize=12)
# ax1.set_xlabel(f"Time (UT)",fontsize=12)
# ax1.tick_params(axis='x', labelsize=12)
# ax1.tick_params(axis='y', labelsize=12)
# ax2.tick_params(axis='y', labelsize=12)
# plt.title('SUIT Excess intensity and GOES $\Delta$ Temperature',fontsize=14)
# ax2.axvline(trigger_time, color="b", ls="--", lw=2,label="HOPE trigger")
# ax2.axvline(impulsive_phase_start, color="b", ls="-", lw=1,label="Impulsive phase start time")
# time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()
# ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left',fontsize=10)

# # plt.legend(loc='upper left',fontsize=10)
# plt.gca().xaxis.set_major_formatter(time_formatter)
# plt.savefig('suit_nb04_goes_diff_temp.png',dpi=300)
# plt.close()

