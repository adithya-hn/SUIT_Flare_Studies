

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import *
from astropy.io import fits
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

Filters=['NB05']
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/light_curve/csv_files/c4_{param}_lc_data.csv',skiprows=1,delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
data1=(np.loadtxt(f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/light_curve/case1_crpix_values_2k.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 

time_array=np.array(data[0], dtype='datetime64')
lc=np.array(data[1],dtype=float)

time_array_2=np.array(data1[0],dtype='datetime64')
lc2=np.array(data1[1],dtype=float)

rc('axes', linewidth=1.2)
plt.rcParams["xtick.major.size"] = 10
fig,axs=plt.subplots(1,1, figsize=(10,5))
ax2 = axs.twinx()
axs.xaxis.set_tick_params(size=0.5)
axs.yaxis.set_tick_params(size=0.5)
axs.tick_params(axis='both', direction='in', length=6, width=1)
axs.tick_params(which='minor', direction='in', length=3, width=1)
axs.yaxis.set_ticks_position('both')
axs.xaxis.set_ticks_position('both')
axs.minorticks_on()


float_array = [float(string) for string in data[1]]
float_array_er = [float(string) for string in data[2]]
y_er=np.std(float_array_er)
axs.plot(time_array,lc,'k',linewidth=0.5,label=f'{param}')
ax2.plot(time_array_2,lc2,label='2K CRPIX Values')
Flt=param
axis_title='Total count'
img_nm=Flt+'_light_curve.png'
axs.set_xlabel('Time',fontsize=13)
ax2.set_ylabel('2K_CRPIX1 value')
axs.set_ylabel(f'{Flt} total count ')
ax2.legend(loc='upper right')
axs.legend(loc='upper left')
plt.title(Flt+' Light Curve'+f' ({str(time_array[0])[:10]})')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show() #close()