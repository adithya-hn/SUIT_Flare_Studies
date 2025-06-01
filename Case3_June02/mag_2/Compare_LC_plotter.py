

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

data=(np.loadtxt(f'NB08_M2.1_Light_curve_data.csv',delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'

NB3_data=(np.loadtxt(f'cabox_900thmagnetogram_M2.1_Light_curve_data.csv',delimiter=',',dtype='str')).transpose()

Helios=(np.load("cdte_data_flare_3.npy", allow_pickle=True)).transpose()
Helios_er=Helios[2]
date_array=data[0]

NB3_date_array=NB3_data[0]

cdte1=Helios[1]
cdte2=Helios[2]
datetime_objects = pd.to_datetime(Helios[0])
helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]


time_array=[]

for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time)

nb3_time_array=[]
for i in range(len(NB3_date_array)):
    parsed_time = datetime.fromisoformat(NB3_date_array[i])
    nb3_time_array.append(parsed_time)

rc('axes', linewidth=1.2)
plt.rcParams["xtick.major.size"] = 10
fig,axs=plt.subplots(1,1, figsize=(11,5))
fig.subplots_adjust(right=0.85)
ax2 = axs.twinx()
#ax3 = axs.twinx()
#ax3.spines.right.set_position(("axes", 1.1))

axs.xaxis.set_tick_params(size=0.5)
axs.yaxis.set_tick_params(size=0.5)
axs.tick_params( axis='both', direction='in', length=6, width=1)
axs.tick_params(which='minor', direction='in', length=3, width=1)
axs.yaxis.set_ticks_position('both')
axs.xaxis.set_ticks_position('both')
axs.minorticks_on()

float_array = np.array([float(string) for string in data[1]])
float_array_er = [float(string) for string in data[2]]
y_err=np.sqrt(float_array/78336)

nb3float_array = [float(string) for string in NB3_data[2]]
#nb3float_array_er = [float(string) for string in NB3_data[2]]

ax2.errorbar(time_array,list(map(int,float_array)),yerr=y_err,fmt='ko-',capsize=2,markersize=2,linewidth=0.5,label='NB08')
axs.plot(nb3_time_array,list(map(int,nb3float_array)),'go-',markersize=2,linewidth=0.5,label='HMI LOS')
#ax3.errorbar(helio_time_array,cdte1,yerr=Helios_er,fmt='ro-',capsize=2,markersize=2,linewidth=0.5,label='Helios-cdte1')
#ax2.plot(nb3_time_array,hmi_data,'bo--',markersize=2,linewidth=0.5)
ax2.set_ylabel('Ca II h Total count ')
axs.set_ylabel("Magnetic flux (Mx)")
#ax3.set_ylabel("HEL1OS count ")
#ax3.set_yscale('log')
#axs.legend(loc='upper right')
#ax2.legend()

m_cls=datetime.fromisoformat('2024-06-02T08:50:00.000')
m_cls_p=datetime.fromisoformat('2024-06-02T08:56:00.000')

Flt=param
axis_title='Total count'
img_nm='Comb_cabox_900G_Light_curve.png'
plt.title('100G and above Light Curve')
#
plt.xlabel('Time',fontsize=13)
plt.axvline(m_cls,color='orange',linestyle='dotted',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))

plt.savefig(img_nm,dpi=300)
plt.show() #close()