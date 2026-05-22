import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
from sunpy import timeseries as ts
from sunpy.data.sample import GOES_XRS_TIMESERIES
import numpy as np
from sunkit_instruments import goes_xrs
import matplotlib.dates as mdates
import datetime
import pandas as pd
from datetime import datetime, timedelta   
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()
# import sc

# 1. Define the time range
start_time = '2024-11-01 01:50'
end_time = '2024-11-01 02:20'

ost="2024-11-01 01:59"
est="2024-11-01 02:15"

# 2. Search GOES data
#result = Fido.search(a.Time(start_time, end_time), a.Instrument('goes'))
#results = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"), a.Resolution("flx1s"))
#result_goes15 = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"), a.goes.SatelliteNumber(15), a.Resolution("flx1s"))

result = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"),a.goes.SatelliteNumber(18),a.Resolution("flx1s"))

# 3. Download the data
downloaded_files = Fido.fetch(result)

# 4. Load as TimeSeries
#goes_ts = TimeSeries(downloaded_files, concatenate=True)

goes_ts = ts.TimeSeries(downloaded_files)

goes_flare = goes_ts.truncate(start_time, end_time)
goes_temp_em =goes_xrs.calculate_temperature_em(goes_flare) #goes_xrs.calculate_temperature_em(goes_flare)


fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
goes_flare.plot(axes=ax1)
goes_temp_em.plot(columns=["temperature"], axes=ax2)
goes_temp_em.plot(columns=["emission_measure"], axes=ax3)
ax3.set_yscale("log")
plt.savefig('GOES_temp_em.png',dpi=300)
plt.close()

fig, ( ax2, ax3) = plt.subplots(2, sharex=True)
goes_temp_em.plot(columns=["temperature"], axes=ax2)
goes_temp_em.plot(columns=["emission_measure"], axes=ax3)
ax3.set_yscale("log")
ax2.grid(True)
ax3.grid(True)
plt.savefig('GOES_temp_em1.png',dpi=300)
plt.close()

#----------------------------------
stix_t=pd.read_csv('7_Nov01_stix_hot_onset_phase.csv')
Solexs=(np.loadtxt('fit_results_AL1_SOLEXS_20241101_SDD2_L1_2411010000_2411010230_TEMP_EM.txt',skiprows=1,dtype='str')).transpose()
time_array4=[]

sl_temp=[float(tp) for tp in Solexs[1]]
sl_temp_er=[float(tpe) for tpe in Solexs[2]]
sl_Em=[float(em) for em in Solexs[3]]
sl_Em_er=[float(eme) for eme in Solexs[4]]
slt=Solexs[0]
time_array4=[datetime.strptime(str(ts)[:19], "%Y-%m-%dT%H:%M:%S") for ts in slt]


mask=(goes_ts.time>=start_time) & (goes_ts.time<=end_time)
gf=goes_temp_em.to_dataframe()

onset_ts=gf[ost:est]
print(onset_ts.columns)
time_nums =  mdates.date2num(onset_ts.index.to_pydatetime())

temp=stix_t['T']
t_dt=stix_t['time_start']
em=stix_t['EM']
t_er_l=stix_t['T_er1']
t_er_h=stix_t['T_er2']
em_er_l=stix_t['EM_er1']
em_er_h=stix_t['EM_er2']


fig,  ax4 = plt.subplots(1, sharex=True, figsize=(10,7))
# ax5.plot(onset_ts['emission_measure'])
ax4.plot(onset_ts['temperature'],label='GOES Temperature')
ax4.set_title('Temperature from GOES - SoLEXS - STIX')
ax4.errorbar(time_array4,sl_temp,yerr=sl_temp_er, marker="o",capsize=2,markersize=2,linewidth=0.5, label='SoLEXS Temperature'); #axs[2].legend(loc='upper center')
# ax4.plot(pd.to_datetime(t_dt),temp)#,xerr=timedelta(seconds=10),yerr=[t_er_h,t_er_l],ls='dotted',color='tab:red',label='Temperature')
ax4.errorbar(pd.to_datetime(t_dt),temp,xerr=timedelta(seconds=10),yerr=[t_er_h,t_er_l],ls='dotted',color='tab:red',label='STIX Temperature')
# ax5.errorbar(pd.to_datetime(t_dt),em,xerr=timedelta(seconds=10),yerr=[em_er_h,em_er_l],ls='dotted',color='tab:blue',label='Emission measure')
# ax5.set_title('Emission measure')
ax4.set_xlabel('Time')
ax4.set_ylabel('Temperature (in MK)')
# ax5.set_ylabel('Emission measure (in $cm^{-3}$)')
plt.legend()

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig('GOES_temp_em1_onset.png',dpi=300)
plt.show()
#------------------------------------------------



# plt.figure(figsize=(8, 6))
# sc = plt.scatter( onset_ts['emission_measure'],onset_ts['temperature'],c=time_nums, cmap='rainbow', edgecolor='k', s=10)

# # Add colorbar with time formatting
# cbar = plt.colorbar(sc, label='Time (UTC)')
# cbar.ax.yaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
# #plt.plot(gf["emission_measure"],gf["temperature"])
# plt.xscale("log")
# plt.ylabel('Temperature (in MK)')
# plt.xlabel('Emission measure ($cm^{-3}$)')
# plt.grid(True)

# plt.savefig('GOES_em_v_temp_onset.png',dpi=300)
# plt.show()
