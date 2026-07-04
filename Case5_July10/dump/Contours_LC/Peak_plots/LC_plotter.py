

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import *
from astropy.io import fits
import scipy.misc
import math as mt
from datetime import datetime
import os
import timeit
from scipy import stats as S
import scipy as sp
import pathlib
import pandas as pd
from subprocess import call
from matplotlib import colors
import matplotlib.dates as mdates

pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'NB08_results.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB08_Light_curve_data.dat'
data2=(np.loadtxt(f'NB03_results.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'

date_array=data[0]
date_array2=data2[0] 


time_array=[]
print(len(date_array))
for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time)


time_array2=[]
print(len(date_array2))
for i in range(len(date_array2)):
    parsed_time = datetime.fromisoformat(date_array2[i])
    time_array2.append(parsed_time)


rc('axes', linewidth=1.2)
plt.rcParams["xtick.major.size"] = 10
plt.rcParams['font.family'] = 'serif'  
plt.rcParams['font.sans-serif'] = 'Times New Roman'
fig,axs=plt.subplots(1,1, figsize=(10,5))

axs.xaxis.set_tick_params(size=0.5)
axs.yaxis.set_tick_params(size=0.5)
axs.tick_params(axis='both', direction='in', length=6, width=1.5)
axs.tick_params(which='minor', direction='in', length=3, width=1.5)
axs.yaxis.set_ticks_position('both')
axs.xaxis.set_ticks_position('both')
axs.minorticks_on()
ax2 = axs.twinx()

float_array = [float(string) for string in data[1]]
float_array2 = [float(string) for string in data2[1]]

contour_area1=2816 #ca_ii_h
contour_area2=7463 #mg_ii_k
float_array=np.array(float_array)
float_array2=np.array(float_array2)
axs.plot(time_array,float_array/contour_area1,'ko-',markersize=2,linewidth=0.5,label='Total count of ca ii h peak time contour')
ax2.plot(time_array2,float_array2/contour_area2,'bo--',markersize=2,linewidth=0.5,label='Total count of mg ii k peak time contour')
'''
ax2.plot(g_time_array,g_float_array,'bo--',markersize=0.1,linewidth=0.5)
ax2.set_ylabel("X-ray flux [1-8 A] (Wm$^{-2}$$s^{-1}$)")
axs.set_ylabel('Total count (NUV)')
ax2.set_yscale('log')'''

m_cls=datetime.fromisoformat('2024-07-10T15:25:00.000')
m_cls_p=datetime.fromisoformat('2024-07-10T15:37:00.000')
#m_cls=datetime.fromisoformat('2024-06-01T08:29:00.000')


axis_title='Total count'
img_nm='NB03_NB08_light_curve.png'

plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
#plt.axvline(m_cls,color='r',label='M class Flare start time',linestyle='dotted')
#plt.axvline(x_cls,color='b',linestyle='dotted',label='GOES Flare start time')
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time',alpha=0.5)
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time',  alpha=0.5)
plt.title('NB03 and NB08 Light Curve (area normalized)')
plt.legend(loc='best')
#plt.axhline(2.58e8,color='g',linestyle='dotted')

plt.legend(loc='best')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show() #close()