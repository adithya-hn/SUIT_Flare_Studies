

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
data1=(np.loadtxt(f'NB03_c3_lc_data.csv',delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
data2=(np.loadtxt(f'NB08_c3_lc_data.csv',delimiter=',',dtype='str')).transpose() #'NB08_Light_curve_data.dat'
data3=(np.loadtxt(f'NB04_c3_lc_data.csv',delimiter=',',dtype='str')).transpose() #'NB04_Light_curve_data.dat'
Solexs=(np.loadtxt(f'AL1_SOLEXS_20240601_SDD2_L1_puc_tb_fit_results_TEMP_EM.txt',delimiter=' ',dtype='str')).transpose()


m_cls=datetime.fromisoformat('2024-06-01T08:25:00.000')
m_cls_p=datetime.fromisoformat('2024-06-01T08:49:00.000')


date_array3=data3[0]
date_array2=data2[0]
date_array1=data1[0]

time_array1=[]

for i in range(len(date_array1)):
    parsed_time = datetime.fromisoformat(date_array1[i])
    time_array1.append(parsed_time)

date=time_array1[0].strftime('%Y-%m-%d')


fig,axs=plt.subplots(1,1, figsize=(10,5))
float_array1 = np.array(data1[1], dtype=np.float64)/np.array(data1[3], dtype=np.float64) # mean counts
float_array_er1 = np.array(data1[2]) # mean QS
qs_area=np.array(data1[4],dtype= float)
axs.errorbar(time_array1,list(map(int,float_array1)),yerr=np.sqrt(float_array1/qs_area),color='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Mg II k light curve')

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
float_array2 = np.array(data2[1],dtype=float)
float_array_er2 = [float(string)  for string in data2[2]]
float_array_er2_=np.std(float_array_er2)*3*np.sqrt(int(data2[3,0]))
axs.errorbar(time_array2,float_array2,yerr=float_array_er2_,fmt='ko-',capsize=2,markersize=2,linewidth=0.5, label='Ca II h light curve')

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

#---------------------------------------------------------%%%%%%-NB04-%%%%%%%------------------------------------------#

time_array3=[]

for i in range(len(date_array3)):
    parsed_time = datetime.fromisoformat(date_array3[i])
    time_array3.append(parsed_time)

fig,axs=plt.subplots(1,1, figsize=(10,5))
float_array3 = np.array(data3[1],dtype=float)  #[float(string)  for string in data3[1]]
float_array_er3 = [float(string)  for string in data3[2]]
float_array_er3_=np.std(float_array_er3)*3*np.sqrt(int(data3[3,0]))
axs.errorbar(time_array3,list(map(int,float_array3)),yerr=float_array_er3_,fmt='co-',capsize=2,markersize=2,linewidth=0.5, label='Mg II h light curve')

img_nm='nb4_lc.png'
plt.ylabel('Total counts',fontsize=13)
plt.xlabel(f"Time ",fontsize=13)
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.title(f'Mg II h light curve ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()


#---------------------------------------------------------%%%%%%-NB02-%%%%%%%------------------------------------------#


#------------------------------------Over plotting all the light curves----------------------------------#
fig4,ax4=plt.subplots(1,1, figsize=(10,5))
ax41 = ax4.twinx()
#ax42 = ax4.twinx()
#ax43 = ax4.twinx()
#ax44 = ax4.twinx()

ax4.errorbar(time_array2,list(map(int,float_array2)),yerr=float_array_er2_,fmt='ko-',capsize=2,markersize=2,linewidth=0.5, label='Ca II h light curve')
ax41.errorbar(time_array1,float_array1,yerr=float_array1/np.sqrt(float_array1),color='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Mg II k light curve')
#ax42.errorbar(time_array3,list(map(int,float_array3)),yerr=float_array_er3_,fmt='co-',capsize=2,markersize=2,linewidth=0.5, label='Mg II h light curve')
##ax43.errorbar(time_array4,list(map(int,float_array4)),yerr=float_array_er4_,fmt='bo-',capsize=2,markersize=2,linewidth=0.5, label='Mg II k wing light curve')
#ax44.errorbar(time_array5,list(map(int,float_array5)),yerr=float_array_er5_,fmt='mo-',capsize=2,markersize=2,linewidth=0.5, label='Mg II h wing light curve')
#ax42.spines.right.set_position(("axes", 1.08))
#ax44.spines.right.set_position(("axes", 1.12))
img_nm='all_lc.png'


#------------------------------------Helios----------------------------------------------------#

Helios=(np.load("cdte_data_flare_1.npy", allow_pickle=True)).transpose()
print(Helios.shape)
cdte1=Helios[1][19:]
cdte2=Helios[2]
#print(np.array(Helios[1]))
cdte1_er=np.sqrt(np.array(Helios[1][19:], dtype=np.float64))
datetime_objects = pd.to_datetime(Helios[0])
helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]
ax43 = ax4.twinx()

#ax43.spines.right.set_position(("axes", 1.15))
ax43.errorbar(helio_time_array[19:],cdte1,yerr=cdte1_er, fmt='ro-',capsize=2,markersize=2,linewidth=0.5,label="Helios-CdTe1",alpha=0.5)
#ax3.plot(helio_time_array,cdte2, label="Helios")

ax41.set_ylabel('Mg II k total counts',fontsize=13)
ax43.set_yscale('log')
ax43.spines.right.set_position(("axes", 1.1))

ax4.set_ylabel('Ca II h total counts',fontsize=13)
ax43.set_ylabel('HEL1OS counts',fontsize=13)
ax4.set_xlabel('Time',fontsize=13)

plt.axvline(m_cls,color='orange',linestyle='--')#,label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-')#,label='GOES Flare peak time')
plt.title(f'Light curves ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close() 


#---------------------------------------------------------%%%%%%%SOLEXS%%%%%%------------------------------------------#

time_array4=[]

sl_temp=[float(tp) for tp in Solexs[1]]
sl_temp_er=[float(tpe) for tpe in Solexs[2]]
sl_Em=[float(em) for em in Solexs[3]]
sl_Em_er=[float(eme) for eme in Solexs[4]]

sltime=np.array([float(tp) for tp in Solexs[0]])
#sl_time=[datetime.strptime(str(ts)[:19], "%Y-%m-%dT%H:%M:%S") for ts in sltime]
base_time = datetime(2024, 6, 1, 7, 0, 0)  # Jun 1, 2025 07:00:00 UTC
time_seconds = sltime-sltime[0]  # Convert string times to float seconds
# Convert seconds to datetime
time_array4 = [base_time + timedelta(seconds=int(t)) for t in time_seconds]


fig,axs=plt.subplots(1,1, figsize=(10,5))
axs1=axs.twinx()
axs.errorbar(time_array4,list(map(int,sl_Em)),yerr=sl_Em_er,fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Emission Measure')
axs1.errorbar(time_array4,list(map(int,sl_temp)),yerr=sl_temp_er,fmt='g-', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Temperature')


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


#-------------------

#fig, axs = plt.subplots(5, 1, sharex=True, figsize=(12,10))
fig, axs = plt.subplots(5, 1, sharex=True, figsize=(10,12),
                        gridspec_kw={'hspace': 0})  # no vertical spacing
for i in range(len(axs)):  # all but bottom panel
    axs[i].ticklabel_format(style='plain', axis='y', scilimits=(0,0))
    #axs[i].yaxis.get_offset_text().set_fontsize(8)  # smaller font
    #axs[i].yaxis.get_offset_text().set_y(-0.35)
    
    if i ==0:
        axs[i].tick_params(axis="x", which="both", bottom=True, top=True) 
    else:
        #axs[i].label_outer()                     # hide x-labels
        axs[i].tick_params(axis="x", which="both", bottom=True, top=False) 
    #axs[i].yaxis.offsetText.set_position((-0.04,-0.1))  # adjust X,Y offset
    axs[i].grid(True, which='major', linestyle='--', alpha=0.6)

axs1_=axs[1].twinx()
axs[0].errorbar(time_array1, float_array1/10e3,yerr=sl_Em_er,fmt='black', marker="o",capsize=2,markersize=2,linewidth=0.5, label="SUIT Ca II H"); axs[0].legend()
axs[1].errorbar(time_array2, float_array2/10e6,yerr=sl_Em_er,fmt='black', marker="o",capsize=2,markersize=2,linewidth=0.5, label="SUIT Mg II h"); axs[1].legend()
axs1_.errorbar(time_array3, float_array3/10e6,yerr=sl_Em_er,fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label="SUIT Mg II k"); axs[1].legend()
axs[2].errorbar(helio_time_array[19:], cdte1,yerr=sl_Em_er,fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label="HEL1OS (CdTe1+ CdTe2)"); axs[2].legend()
#axs[4].plot(time_array2, float_array2, label="SoLEXS Temp")
#axs[4].plot(time_array2, float_array2, label="SoLEXS EM")
axs[3].errorbar(time_array4,list(map(int,sl_Em)),yerr=sl_Em_er,fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Emission Measure'); axs[3].legend()
axs[4].errorbar(time_array4,list(map(int,sl_temp)),yerr=sl_temp_er,fmt='g-', marker="o",capsize=2,markersize=2,linewidth=0.5, label='Temperature')
axs[4].legend()

axs[0].set_ylabel('Ca II H counts (x$10^{3}$)')
axs[1].set_ylabel('Mg II h counts (x$10^{6}$)')
axs1_.set_ylabel('Mg II k counts (x$10^{6}$)')
axs[2].set_ylabel('HEL1OS counts')
axs[3].set_ylabel('EM(x$10^{43}cm^{-3}$)')
axs[4].set_ylabel('Temperature (MK)')

axs[2].set_yscale('log')
axs[3].set_yscale('log')
axs1_.ticklabel_format(style='plain', axis='y', scilimits=(0,0))
axs[-1].set_xlabel("Time [HH:MM]") # Shared x-label


time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)

# Add global title
fig.suptitle(f"Light Curves of Flare ({date})", fontsize=14, weight='bold')
# Adjust layout so title doesn’t overlap
plt.subplots_adjust(top=0.95)
plt.savefig('test.png',dpi=300)
plt.close()
