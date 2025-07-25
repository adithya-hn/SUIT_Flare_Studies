

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

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

#pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 

data1=(np.loadtxt(f'131_roi_lc.csv',delimiter=',',dtype='str')).transpose()
data2=(np.loadtxt(f'171_roi_lc.csv',delimiter=',',dtype='str')).transpose()
#data3=(np.loadtxt(f'193_roi_lc.csv',delimiter=',',dtype='str')).transpose()
data4=(np.loadtxt(f'94_roi_lc.csv',delimiter=',',dtype='str')).transpose()
data5=(np.loadtxt(f'1600_roi_lc.csv',delimiter=',',dtype='str')).transpose()
#data4=(np.loadtxt(f'NB08_c2_lc_data.csv',delimiter=',',dtype='str')).transpose()

csv_filename ="goes_xray_lightcurve_20240601_850.csv"
df = pd.read_csv(csv_filename, parse_dates=['Time'], index_col='Time')

#data5=(np.loadtxt(f'171_lc.csv',delimiter=',',dtype='str')).transpose()
#data6=(np.loadtxt(f'171_lc.csv',delimiter=',',dtype='str')).transpose()

Helios=(np.load("cdte_data_flare_1.npy", allow_pickle=True)).transpose()

dt1=[datetime.fromisoformat(ts) for ts in data1[0]]
dt2=[datetime.fromisoformat(ts) for ts in data2[0]]
#dt3=[datetime.fromisoformat(ts) for ts in data3[0]]
dt4=[datetime.fromisoformat(ts) for ts in data4[0]]
dt5=[datetime.fromisoformat(ts) for ts in data5[0]]

count1=data1[1].astype(float)
count2=data2[1].astype(float)
#count3=data3[1].astype(float)
count4=data4[1].astype(float)
count5=data5[1].astype(float)


cdte1=Helios[1].astype(float)
cdte2=Helios[2]
datetime_objects = pd.to_datetime(Helios[0])
helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]

rc('axes', linewidth=1.2)
plt.rcParams["xtick.major.size"] = 10
fig,axs=plt.subplots(1,1, figsize=(11,5))
fig.subplots_adjust(right=0.85)
ax2 = axs.twinx()
ax3 = axs.twinx()
ax4 = axs.twinx()
ax5 = axs.twinx()
#ax6 = axs.twinx()
area=9412922.3144


ax3.errorbar(dt1[0:int(len(dt1)*0.75)],count1[0:int(len(dt1)*0.75)],yerr=np.sqrt(count1[0:int(len(dt1)*0.75)])/area,fmt='tab:blue',capsize=2,markersize=1,linewidth=0.5,label='131 roi')
axs.errorbar(dt2[0:int(len(dt2)*0.75)],count2[0:int(len(dt2)*0.75)],yerr=np.sqrt(count2[0:int(len(dt2)*0.75)]),capsize=2,fmt='tab:green',markersize=1,linewidth=0.5,label='171 roi')
ax2.errorbar(dt4[0:int(len(dt4)*0.75)],count4[0:int(len(dt4)*0.75)],yerr=np.sqrt(count4[0:int(len(dt4)*0.75)]),capsize=2,fmt='tab:purple',markersize=1,linewidth=0.5,label='94 roi')
ax5.errorbar(dt5[0:int(len(dt5)*0.75)],count5[0:int(len(dt5)*0.75)],yerr=np.sqrt(count5[0:int(len(dt5)*0.75)]),fmt='ko-',capsize=2,markersize=1,linewidth=0.5,label='1600 roi')
ax4.errorbar(helio_time_array,cdte1,yerr=np.sqrt(cdte1),fmt='ro-',capsize=2,markersize=1,linewidth=0.5,label='HEL1OS')
#ax5.plot(df.index, df['xrsa'], label='0.5–4 Å (short)')
#ax6.plot(df.index, df['xrsb'], label='1–8 Å (long)', color='blue')

ax4.set_yscale('log')
ax3.set_yscale('log')
axs.set_yscale('log')
ax2.set_yscale('log')
axs.set_yscale('log')
#ax6.set_yscale('log')

m_cls=datetime.fromisoformat('2024-06-01T08:25:00.000')
m_cls_p=datetime.fromisoformat('2024-06-01T08:49:00.000')

axis_title='Total count'
img_nm='chop_comb_lc.png'
plt.title('ROI Light Curve')
#
plt.xlabel('Time',fontsize=13)
plt.axvline(m_cls,color='b',linestyle='dotted',label='GOES Flare start time')
plt.axvline(m_cls_p,color='b',linestyle='-',label='GOES Flare peak time')
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))

plt.savefig(img_nm,dpi=300)
plt.show()