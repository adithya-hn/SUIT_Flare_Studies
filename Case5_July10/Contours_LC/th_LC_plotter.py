
#Plotting code for threshold light curve
#Author: Adithya HN
# To plot area normalized light curve for a given threshold value



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
import matplotlib.dates as mdates
from Plots_styl import set_publication_style


Filters=['NB08']
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'NB08_threshold_count.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'
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
plt.rcParams['font.family'] = 'serif'  
plt.rcParams['font.sans-serif'] = 'Times New Roman'
plt.rcParams["xtick.major.size"] = 10

fig,axs=plt.subplots(1,1, figsize=(10,5))
#ax2 = axs.twinx()

axs.xaxis.set_tick_params(size=0.5)
axs.yaxis.set_tick_params(size=0.5)
axs.tick_params(axis='both', direction='in', length=6, width=1.5)
axs.tick_params(which='minor', direction='in', length=3, width=1.5)
axs.yaxis.set_ticks_position('both')
axs.xaxis.set_ticks_position('both')
axs.minorticks_on()

float_array1 = np.array([float(string) for string in data[1]]) # area
float_array2 = np.array([float(string) for string in data[2]])
float_array3 = np.array([float(string) for string in data[3]])

float_array_area1 =np.array( [float(string) for string in data[4]]) #total count
float_array_area2 = np.array([float(string) for string in data[5]])
float_array_area3 = np.array([float(string) for string in data[6]])


Th=4300 #[3900, 4100, 4300]


#float_array_er = [float(string) for string in data[2]]
#y_er=np.std(float_array_er)
#axs.plot(time_array,float_array_area1/float_array1,'ko-',markersize=2,linewidth=0.5)
#axs.plot(time_array,float_array_area2/float_array2,'bo-',markersize=2,linewidth=0.5)
axs.plot(time_array,float_array_area3/float_array3,'ko-',markersize=2,linewidth=0.5,label='Ca II h Normailzed counts')
'''
ax2.plot(g_time_array,g_float_array,'bo--',markersize=0.1,linewidth=0.5)
ax2.set_ylabel("X-ray flux [1-8 A] (Wm$^{-2}$$s^{-1}$)")
axs.set_ylabel('Total count (NUV)')
ax2.set_yscale('log')'''

m_cls=datetime.fromisoformat('2024-07-10T15:25:00.000')
m_cls_p=datetime.fromisoformat('2024-07-10T15:37:00.000')

Flt=param
axis_title='Total count'
img_nm=Flt+f'_{Th}_Threshold_light_curve.png'

#plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
plt.ylabel('Normalized total count',fontsize=13)
plt.axvline(m_cls,color='b',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='b',linestyle='-',label='GOES Flare peak time')
plt.title(f'Ca II h Threshold Light Curve ({Th})')
plt.legend()
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show() #close()