

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

Filters=['NB08']
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'nb08_contours.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'
date_array=data[0] 

time_array=[]
print(len(date_array))
for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time)
#g_time_array=[]
'''
for i in range(len(goes_date)):
    parsed_time = datetime.fromisoformat(goes_date[i])
    g_time_array.append(parsed_time)'''

rc('axes', linewidth=1.2)
plt.rcParams["xtick.major.size"] = 10
plt.rcParams['font.family'] = 'serif'  
plt.rcParams['font.sans-serif'] = 'Times New Roman'
fig,axs=plt.subplots(1,1, figsize=(10,5))
ax2 = axs.twinx()
axs.xaxis.set_tick_params(size=0.5)
axs.yaxis.set_tick_params(size=0.5)
axs.tick_params(axis='both', direction='in', length=6, width=1.5)
axs.tick_params(which='minor', direction='in', length=3, width=1.5)
axs.yaxis.set_ticks_position('both')
axs.xaxis.set_ticks_position('both')
axs.minorticks_on()


float_array = [float(string) for string in data[4]]
float_array_er = [float(string) for string in data[2]]
#y_er=np.std(float_array_er)

contour_area=7463 #ca_ii_h
float_array=np.array(float_array)
axs.plot(time_array,float_array,'ko-',markersize=2,linewidth=0.5,label='CMD_EXPT')

ax2.plot(time_array,float_array_er,'bo--',markersize=1,linewidth=0.5)
ax2.set_ylabel("QS mean")


#m_cls=datetime.fromisoformat('2024-07-10T15:25:00.000')
#m_cls_p=datetime.fromisoformat('2024-07-10T15:37:00.000')
#m_cls=datetime.fromisoformat('2024-06-01T08:29:00.000')


Flt=param
axis_title='Total count'
img_nm=Flt+'Check_light_curve.png'

axs.set_ylabel('Exposure',fontsize=13)
axs.set_xlabel('Time (2024-oct-08)',fontsize=13)
#plt.axvline(m_cls,color='r',label='M class Flare start time',linestyle='dotted')
#plt.axvline(x_cls,color='b',linestyle='dotted',label='GOES Flare start time')
#plt.axvline(m_cls,color='b',linestyle='--',label='GOES Flare start time')
#plt.axvline(m_cls_p,color='b',linestyle='-',label='GOES Flare peak time')
#plt.axhline(2.58e8,color='g',linestyle='dotted')
plt.title('Mean QS Light Curve Ca II h')
plt.legend(loc='best')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show() #close()

