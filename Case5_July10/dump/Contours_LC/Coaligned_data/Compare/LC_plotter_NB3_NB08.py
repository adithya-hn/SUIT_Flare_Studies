

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
import matplotlib.dates as mdates
from plots_styl import set_pub_style
set_pub_style()

Filters=['NB08']
param1='magnetogram'
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'NB08_results.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'
date_array=data[0]

NB3_data=(np.loadtxt(f'NB03_results.csv',delimiter=',',skiprows=1,dtype='str')).transpose()
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

fig,axs=plt.subplots(1,1, figsize=(11,5))
#fig.subplots_adjust(right=0.83)
ax2 = axs.twinx()
float_array = [float(string) for string in data[1]]

nb3float_array = [float(string) for string in NB3_data[1]]

axs.plot(time_array,list(map(int,float_array)),'ko-',markersize=2,linewidth=0.5,label='Ca II h')
ax2.plot(nb3_time_array,list(map(int,nb3float_array)),'bo-',markersize=2,linewidth=0.5,label='Mg II k')
#ax2.plot(nb3_time_array,hmi_data,'bo--',markersize=2,linewidth=0.5)
ax2.set_ylabel("Mg II k Total count ")
axs.set_ylabel('Ca II h Total count ')
axs.set_xlabel('Time')

x_cls=datetime.fromisoformat('2024-07-10T15:25:00.000')
x_cls_p=datetime.fromisoformat('2024-07-10T15:37:00.000')

#axs2[0,0].plot(AR_I,AR_M,'ko',markersize=1.5)
Flt=param
axis_title='Total count'
img_nm='light_curve.png'

#plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
#plt.axvline(m_cls,color='r',label='M class Flare start time',linestyle='dotted')
plt.axvline(x_cls,color='orange',linestyle='dotted',label='GOES Flare start time')
#plt.axvline(m_cls_p,color='r',linestyle='-',label='M class Flare peak time')
plt.axvline(x_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
#plt.axhline(2.58e8,color='g',linestyle='dotted')
plt.title('Light Curves')
plt.figlegend(bbox_to_anchor=(0.0001, 0.35, 0.34, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)

plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show() #close()