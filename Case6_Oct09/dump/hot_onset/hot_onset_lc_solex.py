

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


data1=(np.loadtxt(f'c6_NB04_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'
Solexs=fits.open('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/light_curve/csv_files/AL1_SOLEXS_20241009_SDD2_L1_puc_tb_fit_results_T_EM.fits') #'NB08_Light_curve_data.dat'




m_cls=datetime.fromisoformat('2024-10-09T01:25:00.000')
m_cls_p=datetime.fromisoformat('2024-10-09T01:56:00.000')
hst=datetime.fromisoformat('2024-10-09T01:20:00.000')
time_array1=np.array(data1[0], dtype='datetime64')
date=str(time_array1[0])[:10]
nb4_count =np.array(data1[1],dtype=float)
nb4_count_er=np.sqrt(nb4_count)
#---------------------------------------------------------%%%%%%%SOLEXS%%%%%%------------------------------------------#
time_array2=[]
sl_temp = Solexs[1].data['temperature']
sl_temp_er= Solexs[1].data['TEMPERATURE_ERR']
time = Solexs[1].data['TIME']
sl_Em = Solexs[1].data['em']
sl_Em_er = Solexs[1].data['em_err']
for t in time:
    time_array2.append(datetime.fromtimestamp(t)-timedelta(hours=5,minutes=30))


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
ax1.errorbar(time_array2,sl_Em,yerr=sl_Em_er,fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Emission Measure')
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