
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import astropy.units as u
from sunpy.map import Map
import glob
import datetime
from datetime import timedelta
import timeit
import pathlib
import numpy as np
from tqdm import tqdm
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

data=(np.loadtxt('QS_data.csv',delimiter=',',dtype='str',skiprows=1)).transpose()

dt=np.array(data[0],dtype='datetime64')
t_start =  np.datetime64("2024-11-13T15:08:00")
t_end   =  np.datetime64("2024-11-13T16:57:00")

mask = (dt >= t_start) & (dt <= t_end)
qs_dt=dt[mask]

qs_data=np.array(data[1:,mask],dtype=float)

plt.figure(figsize=(14,8))

plt.errorbar(qs_dt,qs_data[3],yerr=qs_data[4],capsize=2,capthick=1,fmt='o-',markersize=1,lw=1,label='Mean value')
plt.plot(qs_dt,qs_data[0],label='Image mode')
plt.plot(qs_dt,qs_data[1],label='qs mode')
plt.plot(qs_dt,qs_data[2],label='median')
plt.axhline(np.mean(qs_data[2]),label=f'Mean of median ({np.mean(qs_data[2]):.2f})',linestyle=':',color='k')
#plt.plot(dt,qs_data[3],label='Mean value')
plt.ylabel('Intensity')
plt.xlabel('Time (UT)')
plt.title('Quiet Sun parameters in pre flare phase')
plt.legend()
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig('QS_estimators.png',dpi=300)
plt.close()

print('mean of median: ',np.mean(qs_data[2]))