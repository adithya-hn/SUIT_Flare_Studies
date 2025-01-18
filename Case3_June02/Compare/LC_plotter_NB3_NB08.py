

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
from subprocess import call
from matplotlib import colors


Filters=['NB08']
param1='magnetogram'
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'NB08_M2.1_Light_curve_data.csv',delimiter=',',dtype='str')).transpose() #'NB03_Light_curve_data.dat'
date_array=data[0]

NB3_data=(np.loadtxt(f'NB03_M2.1_Light_curve_data.csv',delimiter=',',dtype='str')).transpose()
NB3_date_array=NB3_data[0]

Helios=(np.load("cdte_data_flare_3.npy", allow_pickle=True)).transpose()

cdte1=Helios[1]
cdte2=Helios[2]
cdte1_er=Helios[3]
datetime_objects = pd.to_datetime(Helios[0])
helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]

#print(helio_time_array)
Solexs=np.loadtxt('fit_results_AL1_SOLEXS_20240602_SDD2_L1_2406020630_2406020915_TEMP_EM.txt',dtype='str',skiprows=1).transpose()

sl_temp=[float(tp) for tp in Solexs[1]]
sl_temp_er=[float(tpe) for tpe in Solexs[2]]
sl_Em=[float(em) for em in Solexs[3]]
sl_Em_er=[float(eme) for eme in Solexs[4]]

sltime=Solexs[0]
sl_time=[datetime.strptime(str(ts)[:19], "%Y-%m-%dT%H:%M:%S") for ts in sltime]


time_array=[]
#print(len(data[1]))
for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time)

nb3_time_array=[]
for i in range(len(NB3_date_array)):
    parsed_time = datetime.fromisoformat(NB3_date_array[i])
    nb3_time_array.append(parsed_time)

rc('axes', linewidth=1.2)
plt.rcParams["xtick.major.size"] = 10
fig,axs=plt.subplots(1,1, figsize=(11,5))
fig.subplots_adjust(right=0.85)
ax2 = axs.twinx()
#ax4 = axs.twinx()
#ax5 = axs.twinx()
ax3 = axs.twinx()

ax3.spines.right.set_position(("axes", 1.1))
ax3.errorbar(helio_time_array,cdte1,yerr=cdte1_er, fmt='r',capsize=2,markersize=2,linewidth=0.5,label="Helios-CdTe1")
#ax3.plot(helio_time_array,cdte2, label="Helios")
ax3.set_ylabel('Helios',fontsize=13)
ax3.set_yscale('log')
'''

ax4.spines.right.set_position(("axes", 1.33))
ax4.errorbar(sl_time,sl_temp,yerr=sl_temp_er, fmt='g',capsize=2,markersize=2,linewidth=0.5,label="Temperature-SoLExs",alpha=0.5)
ax4.set_ylabel('Temperature',fontsize=13)
ax4.set_yscale('log')

#ax5.spines.right.set_position(("axes", 1.48))
#ax5.errorbar(sl_time,sl_Em,yerr=sl_Em_er, fmt='gray',capsize=2,markersize=2,linewidth=0.5,label="EM-SoLExs",alpha=0.5)
ax5.plot(sl_time,sl_Em,'gray',linewidth=0.5,label="EM-SoLExs",alpha=0.8)
ax5.set_ylabel('Emision measure',fontsize=13)
ax5.set_yscale('log')'''


axs.xaxis.set_tick_params(size=0.5)
axs.yaxis.set_tick_params(size=0.5)
axs.tick_params( axis='both', direction='in', length=6, width=1)
axs.tick_params(which='minor', direction='in', length=3, width=1)
axs.yaxis.set_ticks_position('both')
axs.xaxis.set_ticks_position('both')
axs.minorticks_on()

float_array = [float(string) for string in data[1]]
float_array_er = [float(string) for string in data[2]]

nb3float_array = [float(string) for string in NB3_data[1]]
nb3float_array_er = [float(string) for string in NB3_data[2]]

y_er=np.std(float_array_er)
nb3_y_er=np.std(float_array_er)

axs.errorbar(time_array,list(map(int,float_array)),yerr=y_er,fmt='ko-',capsize=2,markersize=2,linewidth=0.5,label='NB08')
ax2.errorbar(nb3_time_array,list(map(int,nb3float_array)),yerr=nb3_y_er,fmt='bo-',capsize=2,markersize=2,linewidth=0.5,label='NB03')
#ax2.plot(nb3_time_array,hmi_data,'bo--',markersize=2,linewidth=0.5)
ax2.set_ylabel("NB03 Total count ")
axs.set_ylabel('NB08 Total count ')
#ax2.set_yscale('log')
#axs.legend(loc='upper right')
#ax2.legend()

m_cls=datetime.fromisoformat('2024-06-02T08:50:00.000')
#m_cls=datetime.fromisoformat('2024-06-01T08:29:00.000')
m_cls_p=datetime.fromisoformat('2024-06-02T08:56:00.000')

#axs2[0,0].plot(AR_I,AR_M,'ko',markersize=1.5)
Flt=param
axis_title='Total count'
img_nm='Light_curve.png'

#plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
#plt.axvline(m_cls,color='r',label='M class Flare start time',linestyle='dotted')
plt.axvline(m_cls,color='b',linestyle='dotted',label='GOES Flare start time')
#plt.axvline(m_cls_p,color='r',linestyle='-',label='M class Flare peak time')
plt.axvline(m_cls_p,color='b',linestyle='-',label='GOES Flare peak time')
#plt.axhline(2.58e8,color='g',linestyle='dotted')
plt.title('Light Curve')
#plt.ylim(57e4,68e4)#(66e4,700000)
#plt.legend(loc='best')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))

#mpld3.save_html(fig, '12th_June_ROI_CRval.html')
plt.savefig(img_nm,dpi=300)
plt.show() #close()