
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
import matplotlib.dates as mdates
import seaborn as sns
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

#palette = sns.color_palette("deep")

pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 

data1=(np.loadtxt(f'MgII_2k_lc.csv',delimiter=',',dtype='str')).transpose() #'NB03 fd'
data2=(np.loadtxt(f'NB03_c6_lc_data.csv',delimiter=',',dtype='str')).transpose() #'NB03 roi'
Helios=(np.load("cdte_data_flare_8.npy", allow_pickle=True)).transpose()

date_array1=data1[0]
date_array2=data2[0]

time_array1=[]
time_array2=[]
for i in range(len(date_array1)):
    parsed_time = datetime.fromisoformat(date_array1[i])
    time_array1.append(parsed_time)

for i in range(len(date_array2)):
    parsed_time = datetime.fromisoformat(date_array2[i])
    time_array2.append(parsed_time)

date=time_array1[0].strftime('%Y-%m-%d')

float_array1 = [float(string)  for string in data1[1]]
#float_array_er1 = [float(string)  for string in data1[2]]
#float_array_er1_=np.std(float_array_er1)*3*np.sqrt(int(data1[3,0]))

float_array2 = [float(string)  for string in data2[1]]
float_array_er2 = [float(string)  for string in data2[2]]
float_array_er2_=np.std(float_array_er2)*3*np.sqrt(int(data2[3,0]))

fig4,ax4=plt.subplots(1,1, figsize=(10,5))
ax41 = ax4.twinx()
#ax42 = ax4.twinx()
#ax45 = ax4.twinx()
#ax44 = ax4.twinx()

#ax4.errorbar(time_array2,list(map(int,float_array2)),yerr=float_array_er2_,fmt='ko-',capsize=2,markersize=2,linewidth=0.5, label='Mg II k ROI light curve')
ax41.plot(time_array1,list(map(int,float_array1)),color='tab:blue', marker="o",markersize=2,linewidth=0.5, label='Mg II k FD light curve')

print(Helios.shape)
cdte1=Helios[1]
cdte2=Helios[2]
#print(np.array(Helios[1]))
cdte1_er=np.sqrt(np.array(Helios[1], dtype=np.float64))
datetime_objects = pd.to_datetime(Helios[0])
helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]
ax43 = ax4.twinx()

img_nm='fd_comp.png'

#ax43.spines.right.set_position(("axes", 1.15))
ax43.errorbar(helio_time_array,cdte1,yerr=cdte1_er, fmt='ro-',capsize=2,markersize=2,linewidth=0.5,label="Helios-CdTe1",alpha=0.5)
#ax3.plot(helio_time_array,cdte2, label="Helios")
ax43.set_ylabel('Helios',fontsize=13)
ax41.set_ylabel('Mg II k total counts',fontsize=13)
ax43.set_yscale('log')
ax43.spines.right.set_position(("axes", 1.1))
plt.title(f'Light curves ({date})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show()