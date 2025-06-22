

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
import mpld3

Filters=['NB08']
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'{param}_M1.2_Light_curve_data.csv',delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
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
axs.errorbar(time_array,list(map(int,float_array)),yerr=float_array_er,fmt='ko-',capsize=2,markersize=2,linewidth=0.5)
'''
ax2.plot(g_time_array,g_float_array,'bo--',markersize=0.1,linewidth=0.5)
ax2.set_ylabel("X-ray flux [1-8 A] (Wm$^{-2}$$s^{-1}$)")
axs.set_ylabel('Total count (NUV)')
ax2.set_yscale('log')'''

#m_cls_p=datetime.fromisoformat('2024-06-01T08:46:00.000')
#x_cls_p=datetime.fromisoformat('2024-06-01T08:48:00.000')
#m_cls=datetime.fromisoformat('2024-06-01T08:29:00.000')
#x_cls=datetime.fromisoformat('2024-06-01T08:25:00.000')
#axs2[0,0].plot(AR_I,AR_M,'ko',markersize=1.5)
m_cls=datetime.fromisoformat('2024-06-02T04:41:00.000')
m_cls_p=datetime.fromisoformat('2024-06-02T04:50:00.000')

Flt=param
axis_title='Total count'
img_nm=Flt+'_light_curve.png'

#plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
#plt.axvline(m_cls,color='r',label='M class Flare start time',linestyle='dotted')
#plt.axvline(x_cls,color='b',linestyle='dotted',label='GOES Flare start time')
plt.axvline(m_cls_p,color='b',linestyle='-',label='GOES Flare peak time')
plt.axvline(m_cls,color='b',linestyle='--',label='GOES Flare start time')
#plt.axhline(2.58e8,color='g',linestyle='dotted')
plt.title(Flt+' Light Curve')
#plt.ylim(57e4,68e4)#(66e4,700000)
#plt.legend(loc='best')

#mpld3.save_html(fig, '12th_June_ROI_CRval.html')
plt.savefig(img_nm,dpi=300)
plt.show() #close()