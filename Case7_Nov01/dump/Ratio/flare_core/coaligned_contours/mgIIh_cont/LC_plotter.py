

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

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

Flt='NB04'
Title='Mg II h'
data=(np.loadtxt(f'nb04_contours.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'


m_cls=datetime.fromisoformat('2024-11-01T02:05:00.000')
m_cls_p=datetime.fromisoformat('2024-11-01T02:16:00.000')
date_array=data[0] 
time_array=[]



print(len(date_array))
for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time)

fig,axs=plt.subplots(1,1, figsize=(10,5))
ax2 = axs.twinx()


float_array = [float(string) for string in data[1]]
float_array_er = [float(string) for string in data[2]]
qs_area=[float(string) for string in data[4]]
float_array=np.array(float_array)
float_array_er=np.array(float_array_er)
qs_err=np.sqrt(float_array_er)/qs_area

axs.plot(time_array,float_array,'ro-',markersize=2,linewidth=0.5,label='Mean count of peak time contour')
ax2.errorbar(time_array,float_array_er,yerr=qs_err,fmt='bo-',markersize=2,linewidth=0.5,label='Mean QS')
axis_title='Total count'
img_nm=Flt+'_qs_ar_light_curve.png'

axs.set_ylabel('Total_count of AR',fontsize=13)
ax2.set_ylabel('Mean qount of QS')
axs.set_xlabel('Time',fontsize=13)
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
plt.axvline(m_cls,color='gray',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='gray',linestyle='-',label='GOES Flare peak time')
plt.title(f'{Title} Light Curve')
#plt.legend(loc='best')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show()

#------------------------- AR -------

fig2,ax3=plt.subplots(1,1, figsize=(10,5))

ax3.errorbar(time_array,float_array,yerr=np.sqrt(float_array),fmt='ro-',markersize=2,linewidth=0.5,label='Total count of peak time contour')

img_nm=Flt+'_ar_lc.png'

ax3.set_ylabel('Total_count of AR',fontsize=13)
ax3.set_xlabel('Time',fontsize=13)
plt.axvline(m_cls,color='gray',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='gray',linestyle='-',label='GOES Flare peak time')
plt.title(f'{Title} Light Curve')
plt.legend(loc='best')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()


#--------------------------QS-----------------------

fig4,ax4=plt.subplots(1,1, figsize=(10,5))

ax4.errorbar(time_array,float_array_er,yerr=qs_err,fmt='bo-',markersize=2,linewidth=0.5,label='Mean QS')

axis_title='Total count'
img_nm=Flt+'_qs_lc.png'

ax4.set_ylabel('Mean qount of QS')
ax4.set_xlabel('Time',fontsize=13)
plt.axvline(m_cls,color='gray',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='gray',linestyle='-',label='GOES Flare peak time')
plt.title(f'{Title} Light Curve')
plt.legend(loc='best')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()
