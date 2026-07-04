

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


#------------Input params------------
flt1='nb04'
flt2='nb03'
fig_title='Mg II h / Mg II k'
img_suffix='c7'

m_cls=datetime.fromisoformat('2024-11-01T02:05:00.000')
m_cls_p=datetime.fromisoformat('2024-11-01T02:16:00.000')

data=(np.loadtxt(f'{flt1}_contours.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB08_Light_curve_data.dat'
NB3_data=(np.loadtxt(f'{flt2}_contours.csv',delimiter=',',skiprows=1,dtype='str')).transpose()

#--------------------------------------


date_array=data[0]
NB3_date_array=NB3_data[0]
area=int(data[3,0])
area_array=(np.array(NB3_data[3])).astype(float)
qs_area=(np.array(NB3_data[4])).astype(float)
print('area:',area)

colors_list = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e', '#000000', '#7f7f7f']

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


float_array = np.array([float(string)/area for string in data[1]])
float_array_er =np.array([float(string) for string in data[2]])


nb3float_array = np.array([float(string)/area for string in NB3_data[1]])
nb3float_array_er =np.array( [float(string) for string in NB3_data[2]])



interp_func = interp1d(t2_stamps, nb3float_array, kind='linear',fill_value="extrapolate")
interp_func_qs = interp1d(t2_stamps,NB3_data[2] , kind='linear',fill_value="extrapolate")
counts_2_interp = interp_func(t1_stamps)
qs_intrep_counts=interp_func_qs(t1_stamps)

fig,axs=plt.subplots(1,1, figsize=(11,5))
#fig.subplots_adjust(right=0.83)
#ax2 = axs.twinx()



qs_sigma_ca=np.sqrt(float_array_er/qs_area[0])
qs_sigma_mg=np.sqrt(qs_intrep_counts/qs_area[0])
ar_sigma_ca=np.sqrt(float_array/qs_area[0])
ar_sigma_mg=np.sqrt(counts_2_interp/qs_area[0])

qs_sigma=(float_array_er/qs_intrep_counts)*np.sqrt((qs_sigma_ca / float_array_er)**2 + (qs_sigma_mg / qs_intrep_counts)**2)
ar_sigma=(float_array/counts_2_interp)*np.sqrt((ar_sigma_ca / float_array)**2 + (ar_sigma_mg / counts_2_interp)**2)

counts_2_interp=np.array(counts_2_interp)
float_array=np.array(float_array)
print(float_array_er[0],qs_intrep_counts[0],float_array_er[0]/qs_intrep_counts[0])
axs.errorbar(time_array,float_array_er/qs_intrep_counts,yerr=qs_sigma,color=colors_list[1],marker='^',markersize=2,linewidth=0.5,label='QS ratio')
axs.errorbar(time_array,(float_array/counts_2_interp),yerr=ar_sigma,color=colors_list[2],marker='s',markersize=2,linewidth=0.5,label='AR ratio ')
img_nm=f'{img_suffix}_ratio_lc_{flt1}_{flt2}.png'
axs.set_ylabel("Count ratio ")



plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')

plt.xlabel('Time',fontsize=13)
plt.title(f'{fig_title} ratio Light Curves')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()

fig2,ax2=plt.subplots(1,1, figsize=(11,5))
enhancement=((float_array_er/qs_intrep_counts)-(float_array/counts_2_interp))/(float_array_er/qs_intrep_counts)
ax2.plot(time_array,enhancement,label='Enhancment')


ax2.set_ylabel('Enhancement ')
axs.set_xlabel('Time')
axis_title='Total count'
img_nm2=f'{img_suffix}_enhnace_lc_{flt1}_{flt2}.png'

plt.xlabel('Time',fontsize=13)
plt.title(f'{fig_title} ratio Light Curves')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm2,dpi=300)
plt.close()

np.savetxt('enhancement.csv',np.c_[date_array,enhancement],delimiter=',',header='date,enhancement', fmt='%s')

#----------------------------AR ratio alone---------------

fig3,ax3=plt.subplots(1,1, figsize=(11,5))
enhancement=((float_array_er/qs_intrep_counts)-(float_array/counts_2_interp))/(float_array_er/qs_intrep_counts)
ax3.errorbar(time_array,(float_array/counts_2_interp),yerr=ar_sigma,color=colors_list[2],marker='s',markersize=2,linewidth=0.5,label='AR ratio ')
ax3.set_ylabel('Ratio')
ax3.set_xlabel('Time')
axis_title='Total count'
img_nm3=f'{img_suffix}_AR_ratio_lc_{flt1}_{flt2}.png'

plt.xlabel('Time',fontsize=13)
plt.title(f'{fig_title} ratio Light Curves')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm3,dpi=300)
plt.close()

#----------------------------QS ratio alone------------------

fig4,ax4=plt.subplots(1,1, figsize=(11,5))
#enhancement=((float_array_er/qs_intrep_counts)-(float_array/counts_2_interp))/(float_array_er/qs_intrep_counts)
ax4.errorbar(time_array,(float_array_er/qs_intrep_counts),yerr=qs_sigma,color=colors_list[1],marker='s',markersize=2,linewidth=0.5,label='QS ratio ')
ax4.set_ylabel('Ratio')
ax4.set_xlabel('Time')
axis_title='Total count'
img_nm4=f'{img_suffix}_QS_ratio_lc_{flt1}_{flt2}.png'

plt.xlabel('Time',fontsize=13)
plt.title(f'{fig_title} ratio Light Curves')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm4,dpi=300)
plt.close()
