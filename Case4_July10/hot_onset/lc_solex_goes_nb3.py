

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import *
from astropy.io import fits
import math as mt
from datetime import datetime, timedelta
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
data1=(np.loadtxt(f'c4_NB04_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'
Solexs=(np.loadtxt(f'fit_results_AL1_SOLEXS_20240710_SDD2_L1_2407100330_2407100630_TEMP_EM.txt',skiprows=1,dtype='str')).transpose() #'NB08_Light_curve_data.dat'
goes=(np.loadtxt("goes_emperature_emission.csv",skiprows=1,delimiter=',',dtype='str')).transpose()

m_cls_p=datetime.fromisoformat('2024-07-10T05:59:00.000')
m_cls=datetime.fromisoformat('2024-07-10T05:44:00.000')



time_array1=np.array(data1[0], dtype='datetime64')


date=str(time_array1[0])[:10]

fig,axs=plt.subplots(1,1, figsize=(10,5))
float_array1 = [float(string)  for string in data1[1]]
float_array_er1_=np.sqrt(float_array1)
axs.errorbar(time_array1,list(map(int,float_array1)),yerr=float_array_er1_,color='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Mg II k light curve')

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

#---------------------------------------------------------%%%%%%%SOLEXS%%%%%%------------------------------------------#

time_array2=[]

sl_temp=[float(tp) for tp in Solexs[1]]
sl_temp_er=[float(tpe) for tpe in Solexs[2]]
sl_Em=[float(em) for em in Solexs[3]]
sl_Em_er=[float(eme) for eme in Solexs[4]]
sltime=Solexs[0]

time_array2=[datetime.strptime(str(ts)[:19], "%Y-%m-%dT%H:%M:%S") for ts in sltime]


fig,axs=plt.subplots(1,1, figsize=(10,5))
axs1=axs.twinx()
axs.errorbar(time_array2,list(map(int,sl_Em)),yerr=sl_Em_er,fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Emission Measure')
axs1.errorbar(time_array2,list(map(int,sl_temp)),yerr=sl_temp_er,fmt='g-', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Temperature')


img_nm='SOLEXs_EM_temp_lc.png'
axs.set_ylabel('Emission measure',fontsize=13)
axs1.set_ylabel('Temperature (in MK)',fontsize=13)
plt.xlabel(f"Time ",fontsize=13)
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.title(f'SoLEXS EM-Temp light curves ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()

#------------------------------------Over plotting all the light curves----------------------------------#
fig4,ax4=plt.subplots(1,1, figsize=(12,5))


ax42 = ax4.twinx()
ax43 = ax4.twinx()
fig4.subplots_adjust(right=0.85)

goes_time=np.array(goes[0],dtype='datetime64')
goes_temp=np.array(goes[1],dtype='float')

ax4.errorbar(time_array1,list(map(int,float_array1)),yerr=float_array_er1_,color='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Mg II k light curve')
ax4.set_ylabel('Mg II k total counts',fontsize=13)
ax4.set_xlabel('Time',fontsize=13)

bright=np.array(data1[3],dtype='float')
ax43.plot(time_array1,bright,color='y',markersize=2,linewidth=1, label='Brightnings')
ax42.errorbar(time_array2,list(map(int,sl_temp)),yerr=sl_temp_er,fmt='g-', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Temperature')

ax42.set_ylabel('Temperature (in MK)',fontsize=13)
#ax42.set_yscale('log')

ax42.spines.right.set_position(("axes", 1.1))

ax42.plot(goes_time,goes_temp,markersize=2,linewidth=0.5,label="GOES-Temperature",alpha=0.5)
#ax43.set_ylabel('Temperature',fontsize=13)
#ax43.spines.right.set_position(("axes", 1.2))
img_nm='all_lc.png'
plt.axvline(m_cls,color='orange',linestyle='--')#,label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-')#,label='GOES Flare peak time')
plt.title(f'Light curves ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
hst=datetime.fromisoformat('2024-07-10T05:40:00.000')
# hed=datetime.fromisoformat('2024-07-10T05:44:00.000')
plt.xlim(hst,m_cls_p)
plt.savefig(img_nm,dpi=300)
plt.close() 