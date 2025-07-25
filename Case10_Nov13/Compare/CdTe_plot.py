

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


from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()


Helios=(np.load("cdte_data_flare_10.npy", allow_pickle=True)).transpose()

fig4,ax43=plt.subplots(1,1, figsize=(12,5))
fig4.subplots_adjust(right=0.85)

print(Helios.shape)
cdte=Helios[1]+Helios[2]

peaks=np.array(argrelextrema(cdte,np.greater)[0])
#print(cdte[peaks],peaks)
peaks_sorted=np.argsort(cdte[peaks])[::-1]
print(peaks_sorted)
plot_pts='no'

top_peaks=peaks[peaks_sorted[:25]]
pks=[1,4,3,10,9,8]
selected_peaks=[top_peaks[i] for i in pks]
cdte_er=np.sqrt(np.array(cdte, dtype=np.float64))
datetime_objects = pd.to_datetime(Helios[0])
helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]



ax43.errorbar(helio_time_array,cdte,yerr=cdte_er, fmt='ro-',capsize=2,markersize=2,linewidth=0.5,alpha=0.5)
ax43.set_ylabel('HEL1OS (CdTe 10-40 keV) counts',fontsize=13)
ax43.set_yscale('log')
ax43.set_xlabel("Time",fontsize=13)

#ax43.spines.right.set_position(("axes", 1.2))

m_cls=datetime.fromisoformat('2024-11-13T16:57:00.000')
m_cls_p=datetime.fromisoformat('2024-11-13T17:08:00.000')

date=helio_time_array[0].strftime('%Y-%m-%d')

img_nm='helios_CdTe.png'
if plot_pts=='yes':
    for p in range(len(top_peaks)):
        plt.text(helio_time_array[top_peaks[p]],cdte[top_peaks[p]],p,color='k')

spaceing=-10
for pt in selected_peaks:
    spaceing+=30
    ax43.annotate(
    f'Time: {helio_time_array[pt].hour:02d}:{helio_time_array[pt].minute:02d} \nCounts:{cdte[pt]:.2f}',
    xy=(helio_time_array[pt], cdte[pt]), xycoords='data',
    xytext=(0, spaceing), textcoords='offset points',
    size=10,
    bbox=dict(boxstyle="round,pad=.5", fc="0.8",alpha=0.2),
    arrowprops=dict(arrowstyle="->",
                    connectionstyle="angle,angleA=0,angleB=-90,rad=10"),alpha=0.6)
plt.axvline(m_cls,color='orange',linestyle='--')#,label='GOES Flare start time')
plt.axvline(m_cls_p,color='orange',linestyle='-')#,label='GOES Flare peak time')
plt.title(f'HEL1OS (CdTe 1+CdTe 2) Light curve  ({date})')
#plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show() 