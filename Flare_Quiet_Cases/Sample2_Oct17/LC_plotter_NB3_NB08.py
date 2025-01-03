

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


Filters=['NB08']
param1='magnetogram'
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'NB08_Quiet_AR_Light_curve_data.csv',delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
date_array=data[0]

NB3_data=(np.loadtxt(f'NB03_Quiet_AR_Light_curve_data.csv',delimiter=',',dtype='str')).transpose()
NB3_date_array=NB3_data[0]

time_array=[]
#print(len(data[1]))
for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time)

nb3_time_array=[]
for i in range(len(NB3_date_array)):
    parsed_time = datetime.fromisoformat(NB3_date_array[i])
    nb3_time_array.append(parsed_time)

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

nb3float_array = [float(string) for string in NB3_data[1]]
nb3float_array_er = [float(string) for string in NB3_data[2]]

axs.plot(time_array,list(map(int,float_array)),'go-',markersize=2,linewidth=0.5,label='NB08')
ax2.plot(nb3_time_array,list(map(int,nb3float_array)),'ro-',markersize=2,linewidth=0.5,label='NB03')
#ax2.plot(nb3_time_array,hmi_data,'bo--',markersize=2,linewidth=0.5)
ax2.set_ylabel("NB03 Total count ")
axs.set_ylabel('NB08 Total count ')
ax2.set_yscale('log')
axs.legend(loc='best')
ax2.legend()
#ax2.set_ylim(7.2e6,7.3e6)
#axs.set_ylim(2.1e8,2.25e8) #NB08
#axs.set_ylim(5.5e8,6.5e8) #NB03

m_cls=datetime.fromisoformat('2024-11-13T16:57:00.000')
#m_cls=datetime.fromisoformat('2024-06-01T08:29:00.000')
m_cls_p=datetime.fromisoformat('2024-11-13T17:08:00.000')

#axs2[0,0].plot(AR_I,AR_M,'ko',markersize=1.5)
Flt=param
axis_title='Total count'
img_nm=Flt+'_light_curve.png'

#plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
#plt.axvline(m_cls,color='r',label='M class Flare start time',linestyle='dotted')
#plt.axvline(m_cls,color='b',linestyle='dotted',label='GOES Flare start time')
#plt.axvline(m_cls_p,color='r',linestyle='-',label='M class Flare peak time')
#plt.axvline(m_cls_p,color='b',linestyle='-',label='GOES Flare peak time')
#plt.axhline(2.58e8,color='g',linestyle='dotted')
plt.title('Light Curve_Overplot')
#plt.ylim(57e4,68e4)#(66e4,700000)
plt.legend(loc='best')

#mpld3.save_html(fig, '12th_June_ROI_CRval.html')
plt.savefig(img_nm,dpi=300)
plt.show() #close()