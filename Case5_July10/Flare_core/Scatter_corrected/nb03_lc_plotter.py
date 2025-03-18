

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


Filters=['NB03']
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'{param}_M1.0_flare_core_Light_curve_data.csv',delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
date_array=data[0] #np.loadtxt(f'{param}_date_data.dat',dtype='str')

time_array=[]
print(len(date_array))
for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time)

fig,axs=plt.subplots(1,1, figsize=(10,5))
float_array = [float(string)  for string in data[1]]
float_array_er = [float(string)  for string in data[2]]
float_array_er_=np.std(float_array_er)*3*np.sqrt(20736)
axs.errorbar(time_array,list(map(int,float_array)),yerr=float_array_er_,fmt='bo',capsize=2,markersize=2,linewidth=0.5, label='Mg II k light curve')

m_cls=datetime.fromisoformat('2024-07-10T15:25:00.000')
m_cls_p=datetime.fromisoformat('2024-07-10T15:37:00.000')

img_nm='_light_curve.png'
plt.ylabel('Total counts',fontsize=13)
plt.xlabel('Time',fontsize=13)
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.title('Mg II k light Curve')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.37, 0.5))

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()