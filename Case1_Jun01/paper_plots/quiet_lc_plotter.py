

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


pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data1=(np.loadtxt(f'NB03_X1.4_Light_curve_data.csv',delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
date_array1=data1[0] #np.loadtxt(f'{param}_date_data.dat',dtype='str')

time_array1=[]

for i in range(len(date_array1)):
    parsed_time = datetime.fromisoformat(date_array1[i])
    time_array1.append(parsed_time)

fig,axs=plt.subplots(1,1, figsize=(10,5))
float_array1 = [float(string)  for string in data1[1]]
float_array_er1 = [float(string)  for string in data1[2]]

axs.plot(time_array1,float_array_er1,'bo-',markersize=2,linewidth=0.5, label='Mg II k light curve')

m_cls=datetime.fromisoformat('2024-06-01T08:44:00.000')
m_cls_p=datetime.fromisoformat('2024-06-01T08:58:00.000')

img_nm='nb3_qs_lc_sc.png'
plt.ylabel('Total counts',fontsize=13)
plt.xlabel('Time',fontsize=13)
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.title('Mg II k quiet sun light curve')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.37, 0.5))

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()

#---------------------------------------------------------%%%%%%%%%%%%%------------------------------------------#


data2=(np.loadtxt(f'NB08_X1.4_Light_curve_data.csv',delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
date_array2=data2[0] #np.loadtxt(f'{param}_date_data.dat',dtype='str')

time_array2=[]

for i in range(len(date_array2)):
    parsed_time = datetime.fromisoformat(date_array2[i])
    time_array2.append(parsed_time)

fig,axs=plt.subplots(1,1, figsize=(10,5))

float_array_er2 = [float(string)  for string in data2[2]]

axs.plot(time_array2,float_array_er2,'ko-',markersize=2,linewidth=0.5, label='Ca II h light curve')

img_nm='nb8_qs_lc_sc.png'
plt.ylabel('Total counts',fontsize=13)
plt.xlabel('Time',fontsize=13)
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.title('Ca II h quiet sun light curve')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.37, 0.5))

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()