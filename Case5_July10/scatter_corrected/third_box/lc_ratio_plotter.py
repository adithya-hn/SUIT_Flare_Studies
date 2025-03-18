

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
from scipy.interpolate import interp1d

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()


colors_list = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e', '#000000', '#7f7f7f']

pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'nb08_m1_box3_lc.csv',delimiter=',',dtype='str')).transpose() #'NB08_Light_curve_data.dat'
NB3_data=(np.loadtxt(f'nb03_m1_box3_lc.csv',delimiter=',',dtype='str')).transpose()

date_array=data[0]
NB3_date_array=NB3_data[0]



time_array=[]
t1_stamps=[]
t2_stamps=[]
#print(len(data[1]))
for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    t1_stamps.append(parsed_time.timestamp())
    time_array.append(parsed_time)

nb3_time_array=[]
for i in range(len(NB3_date_array)):
    parsed_time = datetime.fromisoformat(NB3_date_array[i])
    t2_stamps.append(parsed_time.timestamp())
    nb3_time_array.append(parsed_time)
#___________________________________________________________

# Interpolating Light Curve 2 to match timestamps of Light Curve 1
print('-------->',len(data[0]),len(data[1]))
interp_func = interp1d(t2_stamps,NB3_data[1] , kind='linear',fill_value="extrapolate")
interp_func_qs = interp1d(t2_stamps,NB3_data[2] , kind='linear',fill_value="extrapolate")
counts_2_interp = interp_func(t1_stamps)
qs_intrep_counts=interp_func_qs(t1_stamps)

fig,axs=plt.subplots(1,1, figsize=(11,5))
#fig.subplots_adjust(right=0.83)
ax2 = axs.twinx()

float_array = [float(string) for string in data[1]]
float_array_er = [float(string) for string in data[2]]


nb3float_array = [float(string) for string in NB3_data[1]]
nb3float_array_er = [float(string) for string in NB3_data[2]]


counts_2_interp=np.array(counts_2_interp)
float_array=np.array(float_array)
print(float_array_er[0],qs_intrep_counts[0],float_array_er[0]/qs_intrep_counts[0])
axs.plot(time_array,float_array_er/qs_intrep_counts,color=colors_list[1],marker='^',markersize=2,linewidth=0.5,label='QS ratio')
axs.plot(time_array,(float_array/counts_2_interp),color=colors_list[2],marker='s',markersize=2,linewidth=0.5,label='Ca II h/Mg II k ')

enhancement=((float_array_er/qs_intrep_counts)-(float_array/counts_2_interp))/(float_array_er/qs_intrep_counts)
ax2.plot(time_array,enhancement,label='Enhancment')

axs.set_ylabel("count Ratio ")
#ax2.set_ylabel('Enhancement ')
axs.set_xlabel('Time')
axis_title='Total count'
img_nm='light_curve.png'

plt.xlabel('Time',fontsize=13)
plt.title('Light Curves')

plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()