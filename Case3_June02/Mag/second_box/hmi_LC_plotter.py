

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

Filters=['magnetogram']
param=Filters[0]

data=(np.loadtxt(f'secn_cabox_900thmagnetogram_M2.1_Light_curve_data.csv',delimiter=',',dtype='str')).transpose()
date_array=data[0]
print(date_array.shape)

time_array=[]
print(len(data))
for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time)

float_array = [float(string) for string in data[3]]

rc('axes', linewidth=1.2)
plt.rcParams["xtick.major.size"] = 10
fig,axs=plt.subplots(1,1, figsize=(10,5))
axs.xaxis.set_tick_params(size=0.5)
axs.yaxis.set_tick_params(size=0.5)
axs.tick_params(axis='both', direction='in', length=6, width=1)
axs.tick_params(which='minor', direction='in', length=3, width=1)
axs.yaxis.set_ticks_position('both')
axs.xaxis.set_ticks_position('both')
axs.minorticks_on()
axs.plot(time_array,float_array,'go--',markersize=2,linewidth=0.5,label='LOS Magnetic field')

x_cls=datetime.fromisoformat('2024-06-02T08:40:00.000')
x_cls_p=datetime.fromisoformat('2024-06-02T08:50:00.000')

#axs2[0,0].plot(AR_I,AR_M,'ko',markersize=1.5)
Flt=param
axis_title='Total count'
img_nm='caboxII_500th_Mag_light_curve.png'
plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
#plt.axvline(m_cls,color='r',label='M class Flare start time',linestyle='dotted')
plt.axvline(x_cls,color='orange',linestyle='dotted',label='GOES Flare start time')
#plt.axvline(m_cls_p,color='r',linestyle='-',label='M class Flare peak time')
plt.axvline(x_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
#plt.axhline(2.58e8,color='g',linestyle='dotted')
plt.title('500G and above Light Curve')
#plt.ylim(57e4,68e4)#(66e4,700000)
plt.legend(loc='best')
# Format x-axis to show only time
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
#plt.yscale('log')
#mpld3.save_html(fig, '12th_June_ROI_CRval.html')
plt.savefig(img_nm,dpi=300)
plt.show() #close()