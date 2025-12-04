
import os
import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
import glob
import datetime
from datetime import timedelta
import timeit
import pathlib
import numpy as np
from tqdm import tqdm


data=(np.loadtxt('QS_data.csv',delimiter=',',dtype='str',skiprows=1)).transpose()

dt=np.array(data[0],dtype='datetime64')
t_start =  np.datetime64("2024-11-13T15:08:00")
t_end   =  np.datetime64("2024-11-13T16:57:00")

mask = (dt >= t_start) & (dt <= t_end)
qs_data=np.array(data[1:,:],dtype=float)

plt.figure(figsize=(14,8))

plt.errorbar(dt,qs_data[3],yerr=qs_data[4],capsize=2,capthick=1,fmt='o-',lw=1,label='Mean value')
plt.plot(dt,qs_data[0],label='Image mode')
plt.plot(dt,qs_data[1],label='qs mode')
plt.plot(dt,qs_data[2],label='median')
#plt.plot(dt,qs_data[3],label='Mean value')
plt.legend()
plt.show()