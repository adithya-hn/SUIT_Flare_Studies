

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
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()


pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'nb08_contours.csv',skiprows=1,delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
date_array=data[0]

NB3_data=(np.loadtxt(f'nb03_contours.csv',skiprows=1,delimiter=',',dtype='str')).transpose()
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

fig,axs=plt.subplots(1,1, figsize=(11,5))
fig.subplots_adjust(right=0.83)
ax2 = axs.twinx()

qs_area=int(data[4][0])

float_array = np.array([float(string) for string in data[1]])
float_array_er = np.array([float(string) for string in data[2]])
#float_array_er=np.sqrt(float_array_er/qs_area)

nb3float_array = np.array([float(string) for string in NB3_data[1]])
nb3float_array_er = np.array([float(string) for string in NB3_data[2]])
#nb3float_array_er=np.sqrt(nb3float_array_er/qs_area)


axs.errorbar(time_array,list(map(int,float_array)),yerr=np.sqrt(float_array),fmt='ko-',capsize=2,markersize=2,linewidth=0.5,label='Ca II h')
ax2.errorbar(nb3_time_array,list(map(int,nb3float_array)),yerr=np.sqrt(nb3float_array),fmt='bo-',capsize=2,markersize=2,linewidth=0.5,label='Mg II k')
#ax2.plot(nb3_time_array,hmi_data,'bo--',markersize=2,linewidth=0.5)
ax2.set_ylabel("Mg II k Total count ")
axs.set_ylabel('Ca II h Total count ')


x_cls=datetime.fromisoformat('2024-07-10T15:25:00.000')
x_cls_p=datetime.fromisoformat('2024-07-10T15:37:00.000')

#axs2[0,0].plot(AR_I,AR_M,'ko',markersize=1.5)

axis_title='Total count'
img_nm='core_nb03_nb08_lc.png'

axs.set_xlabel('Time',fontsize=13)
plt.axvline(x_cls,color='orange',linestyle='dotted',label='GOES Flare start time')
plt.axvline(x_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.title('Ca II h and Mg II k light curves')
plt.figlegend(bbox_to_anchor=(0.00001, 0.38, 0.34, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()

fig2,axs2=plt.subplots(1,1, figsize=(11,5))
fig2.subplots_adjust(right=0.83)
ax2 = axs2.twinx()
axs2.errorbar(time_array,float_array_er,yerr=np.sqrt(float_array_er/qs_area),fmt='ko-',capsize=2,markersize=2,linewidth=0.5,label='Ca II h QS')
ax2.errorbar(nb3_time_array,nb3float_array_er,yerr=np.sqrt(nb3float_array_er/qs_area),fmt='bo-',capsize=2,markersize=2,linewidth=0.5,label='Mg II k QS')
axs2.set_ylabel('Ca II h mean count ')
ax2.set_ylabel("Mg II k mean count ")
axs2.set_xlabel('Time',fontsize=13)
plt.axvline(x_cls,color='orange',linestyle='dotted',label='GOES Flare start time')
plt.axvline(x_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.title('Ca II h and Mg II k QS light curves')
plt.figlegend(bbox_to_anchor=(0.00001, 0.38, 0.34, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig('core_nb03_nb08_qs_lc.png',dpi=300)
plt.close()