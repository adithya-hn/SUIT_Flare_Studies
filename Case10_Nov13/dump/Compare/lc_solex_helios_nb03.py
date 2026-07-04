

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import *
from astropy.io import fits
import scipy.misc
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
data1=(np.loadtxt(f'NB04_c10_lc_data.csv',delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
Solexs=(np.loadtxt(f'fit_results_AL1_SOLEXS_20241113_SDD2_L1_2411131500_2411131730_TEMP_EM.txt',skiprows=1,dtype='str')).transpose() #'NB08_Light_curve_data.dat'
Helios=(np.load("cdte_data_flare_10.npy", allow_pickle=True)).transpose()


m_cls=datetime.fromisoformat('2024-11-13T16:57:00.000')
m_cls_p=datetime.fromisoformat('2024-11-13T17:08:00.000')


date_array1=data1[0]

time_array1=[]

for i in range(len(date_array1)):
    parsed_time = datetime.fromisoformat(date_array1[i])
    time_array1.append(parsed_time)

date=time_array1[0].strftime('%Y-%m-%d')

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
#axs.errorbar(time_array2,list(map(int,sl_Em)),yerr=sl_Em_er,fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Emission Measure')
axs.errorbar(time_array2,list(map(int,sl_temp)),yerr=sl_temp_er,fmt='go',linewidth=0.5, capsize=2,markersize=2, label='Temperature')


img_nm='SOLEXs_EM_temp_lc.png'
#axs.set_ylabel('Emission measure',fontsize=13)
axs.set_ylabel('Temperature (in MK)',fontsize=13)
plt.xlabel(f"Time ",fontsize=13)
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.title(f'SoLEXS Temp light curves ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()

#------------------------------------Over plotting all the light curves----------------------------------#
fig4,ax4=plt.subplots(1,1, figsize=(12,5))


ax42 = ax4.twinx()
#ax43 = ax4.twinx()
fig4.subplots_adjust(right=0.85)

print(Helios.shape)
cdte1=Helios[1]+Helios[2]

#print(np.array(Helios[1]))
cdte1_er=np.sqrt(np.array(Helios[1], dtype=np.float64))
datetime_objects = pd.to_datetime(Helios[0])
helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]

#ax4.errorbar(time_array1,list(map(int,float_array1)),yerr=float_array_er1_,color='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Mg II k light curve')
ax4.set_ylabel('Mg II k total counts',fontsize=13)
ax4.set_xlabel('Time',fontsize=13)

'''ax41 = ax4.twinx()
ax41.errorbar(time_array2,list(map(int,sl_Em)),yerr=sl_Em_er,fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Emission Measure')
ax41.set_ylabel('Emission measure (in $10^{46} cm^{-3}$)',fontsize=13)
ax41.set_yscale('log')'''

ax42.errorbar(time_array2,list(map(int,sl_temp)),yerr=sl_temp_er,fmt='g-', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Temperature')
ax42.set_ylabel('Temperature (in MK)',fontsize=13)
#ax42.spines.right.set_position(("axes", 1.1))
#ax42.set_yscale('log')

ax4.errorbar(helio_time_array,cdte1,yerr=cdte1_er, fmt='ro-',capsize=2,markersize=2,linewidth=0.5,label="Helios-CdTe1",alpha=0.5)
ax4.set_ylabel('HEL1OS counts',fontsize=13)
ax4.set_yscale('log')
#ax43.spines.right.set_position(("axes", 1.1))

img_nm='all_lc.png'
plt.axvline(m_cls,color='orange',linestyle='--')#,label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-')#,label='GOES Flare peak time')
plt.title(f'Light curves ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close() 