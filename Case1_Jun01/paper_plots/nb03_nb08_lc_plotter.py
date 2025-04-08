

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
data=(np.loadtxt(f'core_nb08_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB08_Light_curve_data.dat'
date_array=data[0]

NB3_data=(np.loadtxt(f'core_nb03_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose()
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
ax2 = axs.twinx()


float_array = [float(string) for string in data[2]]
nb3float_array = [float(string) for string in NB3_data[2]]

axs.plot(time_array,list(map(int,float_array)),'ko-',markersize=2,linewidth=0.5,label='Ca II h')
ax2.plot(nb3_time_array,list(map(int,nb3float_array)),'bo-',markersize=2,linewidth=0.5,label='Mg II k')

ax2.set_ylabel("Mg II k total count ")
axs.set_ylabel('Ca II h total count ')
axs.set_xlabel('Time')
axs.set_ylim(2580,2980)
ax2.set_ylim(4000,4400)
x_cls=datetime.fromisoformat('2024-07-10T15:25:00.000')
x_cls_p=datetime.fromisoformat('2024-07-10T15:37:00.000')


img_nm='Qs_light_curve.png'

plt.axvline(x_cls,color='orange',linestyle='dotted',label='GOES Flare start time')
plt.axvline(x_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.title('Ca II h and Mg II k quiet sun light curves')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()