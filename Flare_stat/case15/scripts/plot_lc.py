

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
import seaborn as sns
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

#palette = sns.color_palette("deep")

pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data1=(np.loadtxt(f'NB03_c15_lc_data.csv',delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
data2=(np.loadtxt(f'NB08_c15_lc_data.csv',delimiter=',',dtype='str')).transpose() #'NB08_Light_curve_data.dat'
data3=(np.loadtxt(f'NB04_c15_lc_data.csv',delimiter=',',dtype='str')).transpose() #'NB04_Light_curve_data.dat'

m_cls=datetime.fromisoformat('2025-02-03T07:32:00.000')
m_cls_p=datetime.fromisoformat('2025-02-03T07:44:00.000')

date_array3=data3[0] #np.loadtxt(f'{param}_date_data.dat',dtype='str')
date_array2=data2[0] #np.loadtxt(f'{param}_date_data.dat',dtype='str')
date_array1=data1[0] #np.loadtxt(f'{param}_date_data.dat',dtype='str')

time_array1=[]

for i in range(len(date_array1)):
    parsed_time = datetime.fromisoformat(date_array1[i])
    time_array1.append(parsed_time)

date=time_array1[0].strftime('%Y-%m-%d')

fig,axs=plt.subplots(1,1, figsize=(10,5))
float_array1 = [float(string)  for string in data1[1]]
float_array_er1 = [float(string)  for string in data1[2]]
float_array_er1_=np.std(float_array_er1)*3*np.sqrt(int(data1[3,0]))
axs.errorbar(time_array1,list(map(int,float_array1)),yerr=float_array_er1_,color='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Mg II k light curve')

img_nm='nb3_core_lc.png'
plt.ylabel('Total counts',fontsize=13)

plt.xlabel(f"Time ",fontsize=13)
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
#print(time_array1[0].strftime('%y-%m-%d'))

plt.title(f'Mg II k light curve ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.9, 0.5))

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()

#---------------------------------------------------------%%%%%%%%%%%%%------------------------------------------#

time_array2=[]

for i in range(len(date_array2)):
    parsed_time = datetime.fromisoformat(date_array2[i])
    time_array2.append(parsed_time)

fig,axs=plt.subplots(1,1, figsize=(10,5))
float_array2 = [float(string)  for string in data2[1]]
float_array_er2 = [float(string)  for string in data2[2]]
float_array_er2_=np.std(float_array_er2)*3*np.sqrt(int(data2[3,0]))
axs.errorbar(time_array2,list(map(int,float_array2)),yerr=float_array_er2_,fmt='ko-',capsize=2,markersize=2,linewidth=0.5, label='Ca II h light curve')

img_nm='nb8_core_lc.png'
plt.ylabel('Total counts',fontsize=13)
plt.xlabel(f"Time ",fontsize=13)
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.title(f'Ca II h light curve ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.9, 0.5))

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()

#---------------------------------------------------------%%%%%%-NB04-%%%%%%%------------------------------------------#

time_array3=[]

for i in range(len(date_array3)):
    parsed_time = datetime.fromisoformat(date_array3[i])
    time_array3.append(parsed_time)

fig,axs=plt.subplots(1,1, figsize=(10,5))
float_array3 = [float(string)  for string in data3[1]]
float_array_er3 = [float(string)  for string in data3[2]]
float_array_er3_=np.std(float_array_er3)*3*np.sqrt(int(data3[3,0]))
axs.errorbar(time_array3,list(map(int,float_array3)),yerr=float_array_er3_,fmt='co-',capsize=2,markersize=2,linewidth=0.5, label='Mg II h light curve')

img_nm='nb4_core_lc.png'
plt.ylabel('Total counts',fontsize=13)
plt.xlabel(f"Time ",fontsize=13)
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.title(f'Mg II h light curve ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.9, 0.5))

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()