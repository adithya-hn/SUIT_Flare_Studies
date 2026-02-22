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

df=pd.read_csv('fit_results_AL1_SOLEXS_20240602_SDD2_L1_2406020630_2406020915_TEMP_EM.txt',delim_whitespace=True,index_col=0,parse_dates=True)
# start_time = datetime.datetime.strptime('2024-06-02 08:40', '%Y-%m-%d %H:%M')
# end_time =   datetime.datetime.strptime('2024-06-02 08:55', '%Y-%m-%d %H:%M')
start_time = datetime.datetime(2024,6,2,8,30)
end_time =   datetime.datetime(2024,6,2,8,55)
impulsive_phase_start = datetime.datetime(2024,6,2,8,48)
print(type(df.index))
print(df.index.dtype)
df.index = pd.to_datetime(df.index, format="ISO8601")#pd.to_datetime(df.index)
# print(pd.to_datetime(df.index))

df1 = df[(df.index >= start_time) & (df.index <= end_time)]

plt.figure(figsize=(8, 6))
time_nums =  mdates.date2num(df1.index)
#print(time_nums)
sc = plt.scatter( df1['EM'],df1['Temperature'],c=time_nums, cmap='rainbow', edgecolor='k', s=10)

# Add colorbar with time formatting
cbar = plt.colorbar(sc, label='Time (UTC)')
cbar.ax.yaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
#plt.plot(gf["emission_measure"],gf["temperature"])
plt.xscale("log")
plt.ylabel('Temperature (in MK)')
plt.xlabel('Emission measure ($cm^{-3}$)')
plt.grid(True)

plt.savefig('solexs_em_v_temp.png',dpi=300)
plt.show()