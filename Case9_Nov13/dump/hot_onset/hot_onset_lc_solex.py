

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
import scienceplots
plt.style.use('science')
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 



data1=(np.loadtxt(f'c9_NB04_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'
Solexs1=fits.open('AL1_SOLEXS_20241112_SDD2_L1_puc_tb_fit_results_T_EM.fits')#[1].data
Solexs2=fits.open('AL1_SOLEXS_20241113_SDD2_L1_puc_tb_fit_results_T_EM.fits')#[1].data

m_cls=datetime.fromisoformat('2024-11-13T00:10:00.000')
m_cls_p=datetime.fromisoformat('2024-11-13T00:22:00.000')
hst=datetime.fromisoformat('2024-11-13T00:05:00.000')
time_array1=np.array(data1[0], dtype='datetime64')
date=str(time_array1[0])[:10]
nb4_count =np.array(data1[1],dtype=float)
nb4_count_er=np.sqrt(nb4_count)
#---------------------------------------------------------%%%%%%%SOLEXS%%%%%%------------------------------------------#

time_array2=[]

sl_temp1 =Solexs1[1].data['temperature']
sl_temp1_er= Solexs1[1].data['TEMPERATURE_ERR']
time1 = Solexs1[1].data['TIME']
s_em1 = Solexs1[1].data['em']
s_em1_er = Solexs1[1].data['em_err']

sl_temp2 = Solexs2[1].data['temperature']
sl_temp2_er= Solexs2[1].data['TEMPERATURE_ERR']
time2 = Solexs2[1].data['TIME']
s_em2 = Solexs2[1].data['em']
s_em2_er = Solexs2[1].data['em_err']

time_array21=[]
time_array22=[]
for t in time1:
    time_array21.append(datetime.fromtimestamp(t)-timedelta(hours=5,minutes=30))
for t in time2:
    time_array22.append(datetime.fromtimestamp(t)-timedelta(hours=5,minutes=30))

time_array2=np.concatenate((time_array21,time_array22))
sl_Em=np.concatenate((s_em1,s_em2))
sl_temp=np.concatenate((sl_temp1,sl_temp2))
sl_temp_er=np.concatenate((sl_temp1_er,sl_temp2_er))
sl_Em_er=np.concatenate((s_em1_er,s_em2_er))
#------------------------------------Over plotting all the light curves----------------------------------#
fig4,ax=plt.subplots(2,1, figsize=(7,5),sharex=True, gridspec_kw={'hspace': 0})
for i in range(len(ax)):  # all but bottom panel
    ax[i].ticklabel_format(style='plain', axis='y', scilimits=(0,0))
    if i ==0:
        ax[i].tick_params(axis="x", which="both", bottom=True, top=True) 
    else:
        ax[i].tick_params(axis="x", which="both", bottom=True, top=False) 
    ax[i].grid(True, which='major', linestyle='--', alpha=0.6)
ax1 = ax[0].twinx()

fig4.subplots_adjust(right=0.85)
ax[1].errorbar(time_array1,nb4_count/1e8,yerr=nb4_count_er/1e8,color='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5, label=r'Mg II h light curve ')
ax[1].set_ylabel(r'Mg II h total counts ($\times 10^8$)')
ax[-1].set_xlabel('Time')
ax[0].errorbar(time_array2,list(map(int,sl_temp)),yerr=sl_temp_er,fmt='g-', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Temperature')
ax1.errorbar(time_array2,list(map(int,sl_Em)),yerr=sl_Em_er,fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Emission Measure')
ax[0].set_ylabel('Temperature (in MK)')
ax1.set_ylabel('Emission measure')
ax1.set_yscale('log')

ax1.set_xlim(hst,m_cls_p)
ax[1].set_xlim(hst,m_cls_p)

img_nm='solexs_hot_onset_lc.png'
plt.axvline(m_cls,color='orange',linestyle='--')
plt.axvline(m_cls_p,color='orange',linestyle='-')
plt.title(f'Light curves ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
# plt.legend()
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close() 