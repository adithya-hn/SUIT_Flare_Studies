

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
import matplotlib.dates as mdates

Filters=['NB08']
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'{param}_M2.1_Light_curve_data.csv',delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
#print(data[0])
#goes_data=(np.loadtxt(f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June01_Flare/Goes_data.csv',delimiter=',',dtype='str')).transpose()
#print(goes_data)
date_array=data[0] #np.loadtxt(f'{param}_date_data.dat',dtype='str')
#goes_date=goes_data[0]
#print(goes_data[1])
#g_float_array = [float(string) for string in goes_data[1]]
#g_data=list(map(int,g_float_array))
#date_array=np.array(date_array)
#print(date_array.shape)

time_date=[]
time_array=[]
#print(len(date_array))
for i in range(len(date_array)):

    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time.timestamp())
    time_date.append(parsed_time)
#g_time_array=[]
'''
for i in range(len(goes_date)):
    parsed_time = datetime.fromisoformat(goes_date[i])
    g_time_array.append(parsed_time)'''



float_array = [float(string) for string in data[1]]
float_array_er = [float(string) for string in data[2]]
y_er=np.std(float_array_er)*3*np.sqrt(63455)
er_bars=np.ones(len(float_array))*y_er
#axs.errorbar(time_array,list(map(int,float_array)),yerr=y_er,fmt='ko',capsize=2,markersize=2,linewidth=0.5,label=Filters[0])
'''
ax2.plot(g_time_array,g_float_array,'bo--',markersize=0.1,linewidth=0.5)
ax2.set_ylabel("X-ray flux [1-8 A] (Wm$^{-2}$$s^{-1}$)")
axs.set_ylabel('Total count (NUV)')
ax2.set_yscale('log')'''


Flt=param
axis_title='Total count'
img_nm=Flt+'_light_curve.png'
#print(np.diff(time_array))

np.savetxt(f'{param}_time_series.csv',np.c_[time_date,time_array,float_array,er_bars],delimiter=',',fmt='%s')