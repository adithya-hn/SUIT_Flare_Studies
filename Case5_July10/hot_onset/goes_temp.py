import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
from sunpy import timeseries as ts
from sunpy.data.sample import GOES_XRS_TIMESERIES
import numpy as np
from sunkit_instruments import goes_xrs
import matplotlib.dates as mdates
import datetime

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()


# 1. Define the time range
start_time = '2024-06-02 6:30'
end_time = '2024-06-02 09:30'
ost="2024-06-02 08:30"
est="2024-06-02 08:45"

# 2. Search GOES data
#result = Fido.search(a.Time(start_time, end_time), a.Instrument('goes'))
#results = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"), a.Resolution("flx1s"))
#result_goes15 = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"), a.goes.SatelliteNumber(15), a.Resolution("flx1s"))

result = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"),a.goes.SatelliteNumber(16),a.Resolution("flx1s"))

# 3. Download the data
downloaded_files = Fido.fetch(result)

# 4. Load as TimeSeries
#goes_ts = TimeSeries(downloaded_files, concatenate=True)

goes_ts = ts.TimeSeries(downloaded_files)

goes_flare = goes_ts.truncate(start_time, end_time)
goes_temp_em =goes_xrs.calculate_temperature_em(goes_flare) #goes_xrs.calculate_temperature_em(goes_flare)
goes_tp_em=goes_temp_em.to_dataframe()
goes_tp_em_clean = goes_tp_em.dropna(subset=['temperature', 'emission_measure'])# Drop rows where either temperature or EM is missing
goes_tp_em_clean.to_csv('goes_emperature_emission.csv', sep=',', index_label='time')
#goes_temp_em.save_csv('GOES_temp_em.csv')

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



mask=(goes_ts.time>=start_time) & (goes_ts.time<=end_time)
gf=goes_temp_em.to_dataframe()
time_nums =  mdates.date2num(gf.index.to_pydatetime())

onset_ts=gf[ost:est]
print(onset_ts.columns)

fig, ( ax4, ax5) = plt.subplots(2, sharex=True, figsize=(10,7))
ax5.plot(onset_ts['emission_measure'])
ax4.plot(onset_ts['temperature'])
ax4.set_title('Temperature')
ax5.set_title('Emission measure')
ax5.set_xlabel('Time')
ax4.set_ylabel('Temperature (in MK)')
ax5.set_ylabel('Emission measure (in $cm^{-3}$)')

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig('GOES_temp_em1_onset.png',dpi=300)
plt.show()


plt.figure(figsize=(8, 6))
sc = plt.scatter( gf['emission_measure'],gf['temperature'],c=time_nums, cmap='rainbow', edgecolor='k', s=10)

# Add colorbar with time formatting
cbar = plt.colorbar(sc, label='Time (UTC)')
cbar.ax.yaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
#plt.plot(gf["emission_measure"],gf["temperature"])
plt.xscale("log")
plt.ylabel('Temperature (in MK)')
plt.xlabel('Emission measure ($cm^{-3}$)')
plt.grid(True)

plt.savefig('GOES_em_v_temp.png',dpi=300)
plt.show()