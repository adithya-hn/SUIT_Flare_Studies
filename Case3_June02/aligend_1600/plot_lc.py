

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import *
from astropy.io import fits
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
data1=(np.loadtxt(f'NB03_c3_lc_data.csv',delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
data2=(np.loadtxt(f'NB08_c3_lc_data.csv',delimiter=',',dtype='str')).transpose() #'NB08_Light_curve_data.dat'



m_cls=datetime.fromisoformat('2024-06-02T08:50:00.000')
m_cls_p=datetime.fromisoformat('2024-06-02T08:56:00.000')


date_array2=data2[0]
date_array1=data1[0]

time_array1=[]

for i in range(len(date_array1)):
    parsed_time = datetime.fromisoformat(date_array1[i])
    time_array1.append(parsed_time)

date=time_array1[0].strftime('%Y-%m-%d')

fig,axs=plt.subplots(1,1, figsize=(10,5))
float_array1 = [float(string)  for string in data1[1]]
float_array_er1 = [float(string)  for string in data1[2]]
float_array_er1_=np.sqrt(int(data1[3,0]))

axs.errorbar(time_array1,list(map(int,float_array1)),yerr=float_array_er1_,color='tab:blue',capsize=2,markersize=3,fmt='o', label='Mg II k light curve')

img_nm='nb3_lc.png'
plt.ylabel('Total counts',fontsize=13)

plt.xlabel(f"Time ",fontsize=13)
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
#print(time_array1[0].strftime('%y-%m-%d'))

plt.title(f'Mg II k light curve ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.4, 0.5))

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
axs.errorbar(time_array2,list(map(int,float_array2)),yerr=float_array_er2_,fmt='ko',capsize=2,markersize=2,linewidth=0.5, label='Ca II h light curve')

img_nm='nb8_lc.png'
plt.ylabel('Total counts',fontsize=13)
plt.xlabel(f"Time ",fontsize=13)
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.title(f'Ca II h light curve ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()

#--------------------------------------------------- Nb03-08----------------------------------
fig3,ax3=plt.subplots(1,1, figsize=(10,5))
ax31 = ax3.twinx()
ax3.errorbar(time_array2,list(map(int,float_array2)),yerr=float_array_er2_,fmt='ko',capsize=2,markersize=2,linewidth=0.5, label='Ca II h light curve')
ax31.errorbar(time_array1,list(map(int,float_array1)),yerr=float_array_er1_,fmt='co',capsize=2,markersize=2,linewidth=0.5, label='Mg II k light curve')
ax31.set_ylabel('Mg II k total counts',fontsize=13)
ax3.set_ylabel('Ca II H total counts',fontsize=13)
ax3.set_xlabel('Time',fontsize=13)
img_nm='mg2Ca2_lc.png'
plt.axvline(m_cls,color='orange',linestyle='--')#,label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-')#,label='GOES Flare peak time')
plt.title(f'Light curves ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close() 



#------------------------------------Over plotting all the light curves----------------------------------#
fig4,ax4=plt.subplots(1,1, figsize=(10,5))
ax41 = ax4.twinx()
#ax42 = ax4.twinx()
#ax45 = ax4.twinx()
#ax44 = ax4.twinx()

ax4.errorbar(time_array2,list(map(int,float_array2)),yerr=float_array_er2_,fmt='ko',capsize=2,markersize=2,linewidth=0.5, label='Ca II h light curve')
#ax41.errorbar(time_array1,list(map(int,float_array1)),yerr=float_array_er1_,color='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Mg II k light curve')
ax41.errorbar(time_array1,list(map(int,float_array1)),yerr=float_array_er1_,fmt='co',capsize=2,markersize=2,linewidth=0.5, label='Mg II k light curve')
#ax45.errorbar(time_array4,list(map(int,float_array4)),yerr=float_array_er4_,fmt='bo-',capsize=2,markersize=2,linewidth=0.5, label='Mg II k wing light curve')
#ax44.errorbar(time_array5,list(map(int,float_array5)),yerr=float_array_er5_,fmt='mo-',capsize=2,markersize=2,linewidth=0.5, label='Mg II h wing light curve')
#ax42.spines.right.set_position(("axes", 1.08))

#ax44.spines.right.set_position(("axes", 1.12))
img_nm='all_lc.png'


#------------------------------------Helios----------------------------------------------------#

Helios=(np.load("cdte_data_flare_3.npy", allow_pickle=True)).transpose()
print('Helios data',Helios.shape)
cdte1=Helios[1]+Helios[2]
helios_er=Helios[3]+Helios[4]
cdte2=Helios[2]
#print(np.array(Helios[1]))
cdte1_er=np.sqrt(np.array(Helios[1], dtype=np.float64))
datetime_objects = pd.to_datetime(Helios[0])
helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]
ax43 = ax4.twinx()
ax43.errorbar(helio_time_array,cdte1,yerr=helios_er, fmt='ro',capsize=2,markersize=2,linewidth=0.5,label="HEL1OS-(CdTe1+CdTe2)",alpha=0.5)
ax43.set_yscale('log')
ax43.set_ylabel('HEL1OS counts',fontsize=13)
ax43.spines.right.set_position(("axes", 1.1))

ax41.set_ylabel('Mg II k total counts',fontsize=13)
ax4.set_ylabel('Ca II H total counts',fontsize=13)
ax4.set_xlabel('Time',fontsize=13)

#ax43.set_axis_off()
#ax42.set_axis_off()

plt.axvline(m_cls,color='orange',linestyle='--')#,label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-')#,label='GOES Flare peak time')
plt.title(f'Light curves ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close() 