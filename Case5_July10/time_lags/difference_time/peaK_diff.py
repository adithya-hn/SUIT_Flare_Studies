
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
from scipy.signal import argrelextrema
from scipy.signal import find_peaks

Helios = (np.loadtxt('helios_CdTe_c5.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
nb4_suit=(np.loadtxt(f'Diff_img_data_NB04.csv',delimiter=',',skiprows=1,dtype='str')).transpose()

cdte=np.array(Helios[1],dtype=float)
cdte_er=np.array(Helios[2],dtype=float)
helios_dt=np.array(Helios[0],dtype='datetime64[us]')


suit_dt=np.array(nb4_suit[0],dtype='datetime64')
nb4_count=np.array(nb4_suit[5],dtype=float)
nb4_er=np.array(nb4_suit[4],dtype=float)

t_start =  np.datetime64("2024-07-10T13:37:00")
t_end   =  np.datetime64("2024-07-10T15:25:00")

mask = (helios_dt >= t_start) & (helios_dt <= t_end)
t_cut = helios_dt[mask]
img_count_cut  = cdte[mask]


fig4,ax43=plt.subplots(1,1, figsize=(12,5))
fig4.subplots_adjust(right=0.85)
ax5=ax43.twinx()
ax43.errorbar(helios_dt,cdte,yerr=cdte_er, fmt='ro-',capsize=2,markersize=2,linewidth=0.5,alpha=0.5)
ax5.errorbar(suit_dt,nb4_count,yerr=nb4_er,fmt='ko-',capsize=2,markersize=2,linewidth=0.5,alpha=0.5)
ax43.set_ylabel('HEL1OS (CdTe 10-40 keV) counts',fontsize=13)
ax43.set_yscale('log')
ax43.set_xlabel("Time",fontsize=13)

ts = pd.to_datetime(helios_dt)
img_nm='helios_CdTe.png'

threshold1=1e0
peaks=find_peaks(img_count_cut,prominence=0.2)
print('Helios peaks: ', len(peaks))
dt=np.array(t_cut[peaks[0]])
pk=np.array(img_count_cut[peaks[0]])
for i in range(len(peaks[0])):
        plt.axvline(t_cut[peaks[0][i]], color='b',alpha=0.2)
#---------------------------------------------------
threshold=2e7
s_mask = (suit_dt >= t_start) & (suit_dt <= t_end)
suit_cut = suit_dt[s_mask]
suit_img_count_cut  = nb4_count[s_mask]
s_peaks=find_peaks(suit_img_count_cut, height=threshold)

nb4_peak_dt=np.array(suit_cut[s_peaks[0]])
nb4_peak_val=np.array(suit_img_count_cut[s_peaks[0]])
helios_peak_dt=np.array(t_cut[peaks[0]])
helios_peak_val=np.array(img_count_cut[peaks[0]])


for i in range(len(s_peaks[0])):
        ax5.axvline(suit_cut[s_peaks[0][i]], color='g',alpha=0.2)
#------------------------------------------------------
plt.title(f'HEL1OS (CdTe 1+CdTe 2) Light curve)')
#plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show() 

np.savetxt('nb4_peaks.csv',np.c_[nb4_peak_dt.astype(str),nb4_peak_val],delimiter=',',header='peak_time,peak_value',comments='',fmt='%s')
np.savetxt('helios_peaks.csv',np.c_[helios_peak_dt.astype(str),helios_peak_val],delimiter=',',header='peak_time,peak_value',comments='',fmt='%s')