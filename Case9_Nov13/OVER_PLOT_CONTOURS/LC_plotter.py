

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


param='caIIh'
fig_title='Ca II h'
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'{param}_contours.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'
date_array=data[0] 

time_array=[]
print(len(date_array))
for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time)

fig,axs=plt.subplots(1,1, figsize=(10,5))
#ax2 = axs.twinx()
axs.xaxis.set_tick_params(size=0.5)
axs.yaxis.set_tick_params(size=0.5)
axs.tick_params(axis='both', direction='in', length=6, width=1.5)
axs.tick_params(which='minor', direction='in', length=3, width=1.5)
axs.yaxis.set_ticks_position('both')
axs.xaxis.set_ticks_position('both')
axs.minorticks_on()

flare_count = [float(string) for string in data[1]]
cont_area=[float(string) for string in data[2]]

flare_count=np.array(flare_count)
cont_area=np.array(cont_area)

plt.plot(time_array,flare_count/cont_area,'ko-',markersize=2,linewidth=0.5,label='Mean count of the threshold 1')
m_cls=datetime.fromisoformat('2024-11-13T00:10:00.000')
m_cls_p=datetime.fromisoformat('2024-11-13T00:16:00.000')
#m_cls=datetime.fromisoformat('2024-06-01T08:29:00.000')

Flt=param
axis_title='Mean count'
img_nm=Flt+'_th1_light_curve.png'

plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
plt.axvline(m_cls,color='b',linestyle='--',label='Onset start')
plt.axvline(m_cls_p,color='b',linestyle='-',label='Impulsive phase start')
plt.title(f'{fig_title} Light Curve')
plt.legend(loc='best')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show() #close()
#---------------------------------------------------------
fig,axs=plt.subplots(1,1, figsize=(10,5))

plage_count = [float(string) for string in data[3]]
plage_area=[float(string) for string in data[4]]


plage_count=np.array(plage_count)
plage_area=np.array(plage_area)

plt.plot(time_array,plage_count/plage_area,'ko-',markersize=2,linewidth=0.5,label='Mean count of threshold 2')

Flt=param
axis_title='Mean count'
img_nm=Flt+'_th2_light_curve.png'

plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
plt.axvline(m_cls,color='b',linestyle='--',label='Onset start')
plt.axvline(m_cls_p,color='b',linestyle='-',label='Impulsive phase start')
plt.title(f'{fig_title} plage light curve')
plt.legend(loc='best')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show() #close()

#-----------------------------------------
