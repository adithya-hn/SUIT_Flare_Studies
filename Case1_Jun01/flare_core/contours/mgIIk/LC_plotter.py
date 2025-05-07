

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

Filters=['NB03']
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'nb03_contours.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'
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
plt.rcParams['font.family'] = 'serif'  
plt.rcParams['font.sans-serif'] = 'Times New Roman'

fig,axs=plt.subplots(1,1, figsize=(10,5))
#ax2 = axs.twinx()
axs.xaxis.set_tick_params(size=0.5)
axs.yaxis.set_tick_params(size=0.5)
axs.tick_params(axis='both', direction='in', length=6, width=1.5)
axs.tick_params(which='minor', direction='in', length=3, width=1.5)
axs.yaxis.set_ticks_position('both')
axs.xaxis.set_ticks_position('both')
axs.minorticks_on()

flare_count = [float(string) for string in data[1]]
cont_area=[float(string) for string in data[2]]

flare_count=np.array(flare_count)
cont_area=np.array(cont_area)

plt.plot(time_array,flare_count/cont_area,'ko-',markersize=2,linewidth=0.5,label='Mean count of peak time contour')
m_cls=datetime.fromisoformat('2024-06-01T08:25:00.000')
m_cls_p=datetime.fromisoformat('2024-06-01T08:49:00.000')
#m_cls=datetime.fromisoformat('2024-06-01T08:29:00.000')

Flt=param
axis_title='Mean count'
img_nm=Flt+'_flare_light_curve.png'

plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
plt.axvline(m_cls,color='b',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='b',linestyle='-',label='GOES Flare peak time')
plt.title('Mg II k Light Curve')
plt.legend(loc='best')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show() #close()
#---------------------------------------------------------
fig,axs=plt.subplots(1,1, figsize=(10,5))

plage_count = [float(string) for string in data[3]]
plage_area=[float(string) for string in data[4]]
qs_count1=[float(string) for string in data[5]]
qs_count2=[float(string) for string in data[6]]
qs_count3=[float(string) for string in data[7]]

plage_count=np.array(plage_count)
plage_area=np.array(plage_area)
qs_count1=np.array(qs_count1)
qs_count2=np.array(qs_count2)
#qs_count3=np.array(qs_count3)
qs_count=(qs_count1)#+qs_count2+qs_count3)/3
           
           

flare_count=np.array(flare_count)
cont_area=np.array(cont_area)

plt.plot(time_array,plage_count/plage_area,'ko-',markersize=2,linewidth=0.5,label='Mean count of peak time contour')

Flt=param
axis_title='Mean count'
img_nm=Flt+'_plage_light_curve.png'

plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
plt.axvline(m_cls,color='b',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='b',linestyle='-',label='GOES Flare peak time')
plt.title('Mg II k plage light curve')
plt.legend(loc='best')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show() #close()

#-----------------------------------------
fig,axs=plt.subplots(1,1, figsize=(10,5))

plt.plot(time_array,qs_count,label='Mean QS')
plt.plot(time_array,qs_count1,label=' QS1')
plt.plot(time_array,qs_count2,label=' QS2')
plt.plot(time_array,qs_count3,label=' QS3')
Flt=param
axis_title='Mean count'
img_nm=Flt+'_qs_light_curve.png'

plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
plt.axvline(m_cls,color='b',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='b',linestyle='-',label='GOES Flare peak time')
plt.title('Mg II k QS light curve')
plt.legend(loc='best')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show() #close()