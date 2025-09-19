

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

Filters=['NB03'] #'magnetogram','Shear','TOTPOT','USJH','Temp']
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=np.loadtxt(f'{param}_M2.1_Light_curve_data.dat') #'NB03_Light_curve_data.dat'
date_array=np.loadtxt(f'{param}_M2.1_date_data.dat',dtype='str')


#date_array=np.array(date_array)
print(date_array.shape)

time_array=[]
print(len(data))
for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time)

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

axs.plot(time_array,data, 'k',linewidth=1)
#m_cls_p=datetime.fromisoformat('2024-06-01T08:46:00.000')
x_cls=datetime.fromisoformat('2024-06-02T08:40:00.000')
#m_cls=datetime.fromisoformat('2024-06-01T08:29:00.000')
x_cls_p=datetime.fromisoformat('2024-06-02T08:50:00.000')

#axs2[0,0].plot(AR_I,AR_M,'ko',markersize=1.5)
Flt=param
axis_title='Total count'
img_nm=Flt+'M2.1_light_curve.png'

plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
#plt.axvline(m_cls,color='r',label='M class Flare start time',linestyle='dotted')
plt.axvline(x_cls,color='b',linestyle='dotted',label='GOES Flare start time')
#plt.axvline(m_cls_p,color='r',linestyle='-',label='M class Flare peak time')
plt.axvline(x_cls_p,color='b',linestyle='-',label='GOES Flare peak time')
plt.title(Flt+' Light Curve')
#plt.ylim(57e4,68e4)#(66e4,700000)
plt.legend(loc='best')
plt.yscale('log')
plt.ylim(220000,350000)
mpld3.save_html(fig, f'2nd_June_M2.1_Flare{param}.html')
plt.savefig(img_nm,dpi=300)
plt.show() #close()