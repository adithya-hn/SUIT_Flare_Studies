

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from astropy.io import fits
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


Filters=['NB03']
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'csv_files/NB03_c2_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'

time_array=np.array(data[0], dtype='datetime64')

rc('axes', linewidth=1.2)
plt.rcParams["xtick.major.size"] = 10
fig,axs=plt.subplots(1,1, figsize=(10,5))
#ax2 = axs.twinx()
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
flare_lc=np.array(data[1],dtype=float)
axs.plot(time_array,flare_lc,markersize=2,linewidth=0.5)

Flt=param
axis_title='Total count'
img_nm=Flt+'_light_curve.png'

plt.xlabel('Time',fontsize=13)
plt.title(Flt+' Light Curve')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show()