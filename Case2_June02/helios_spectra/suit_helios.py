#

'''
Created on 31 Dec 2025
@author: adithya-hn
purpose: combined plot of suit and helios to see matched ligt curves

'''
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.dates as mdates
import pandas as pd
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
import numpy.ma as ma
import seaborn as sns
import glob
import os
from scipy.optimize import curve_fit
from scipy.ndimage import median_filter

scol =sns.color_palette("colorblind")


suit_fl="Diff_img_data_NB04.csv"
helios_fl="helios_CdTe_c2.csv"

all_times = []
all_ints = []


data = np.genfromtxt(suit_fl, delimiter=',', names=True, dtype=None)
helios=np.genfromtxt(helios_fl, delimiter=',', names=True, dtype=None, encoding='utf-8')

# Convert time to datetime64
t = np.array(data['Date'], dtype='datetime64[ms]')
ht= np.array(helios['Time'], dtype='datetime64[ms]')
h_count=np.array(helios['Total'],dtype=float)
h_er=np.array(helios['CdTeEr'])#np.sqrt(h_count)
area=np.array(data['area'],dtype=float)
#x = np.array(data['diff_count'], dtype=float)#/area
x = np.array(data['diff_count'], dtype=float)#/area

case_id=os.path.basename(suit_fl)[0:3]
flt=os.path.basename(suit_fl)[-8:-4]

t_start =  np.datetime64("2024-06-02T02:50:00")
t_end   =  np.datetime64("2024-06-02T04:41:00")
#norm_v=norm[0]


mask = (t >= t_start) & (t <= t_end)
t_cut = t[mask]
img_count_cut  = x[mask]#/norm_v

window = 15  # choose based on timescale over which you expect NO flares
bkg = median_filter(h_count, size=window)

# Avoid zeros in background to prevent division by zero
bkg_safe = np.clip(bkg, 1e-6, None)
sigma = np.sqrt(bkg_safe)
z = (h_count - bkg_safe) / (sigma/np.sqrt(window))    # "Gaussian-equivalent" significance per bin
sigma_thresh = 1  # or 4, 5, etc.
min_separation_bins = 2  # minimum distance between distinct peaks

peaks, props = find_peaks(z,height=sigma_thresh, distance=min_separation_bins)


peak_times_ = (ht[peaks].astype('datetime64[s]')).astype(str)
peak_times = ht[peaks]
peak_counts = h_count[peaks]

fig,ax=plt.subplots(figsize=(14,6))
ax2=plt.twinx()
threshold=1.0e5 #5661*4 #2.5e6 #
peaks2=find_peaks(img_count_cut)
dt2=(t_cut[peaks2[0]].astype('datetime64[s]')).astype(str)
pk2=np.array(img_count_cut[peaks2[0]])

dt=np.array(ht[peaks[0]])
pk=np.array(h_count[peaks[0]])

all_times.append(dt)
all_ints.append(pk)
data_v2 = np.genfromtxt('v2_Diff_img_data_NB04.csv', delimiter=',', names=True, dtype=None)
print('Comuns headers:', data_v2.dtype.names)
v2_count=(np.array(data_v2['diff_count'], dtype=float))
v2_dt=np.array(data_v2['date_time'], dtype='datetime64[ms]')
ax2.plot(v2_dt,v2_count,'g--')
ax.errorbar(ht,h_count,yerr=h_er)
ax2.plot(t_cut,img_count_cut,'ro-',markersize=2)

for i in range(len(peaks)):
    plt.axvline(ht[peaks[i]], color='b',alpha=0.2)

for i in range(len(peaks2[0])):
    plt.axvline(t_cut[peaks2[0][i]],color='r',alpha=0.2)

np.savetxt('helios_peaks.csv',np.c_[peak_times_,peak_counts],header='date_time,helio_count',comments='',delimiter=',',fmt='%s')
np.savetxt('suit_diff_peaks.csv',np.c_[dt2,pk2],header='date_time,suit_coun',comments='',delimiter=',',fmt='%s')



plt.title(f'{case_id}')
ax.set_yscale('log')
ax2.set_yscale('log')
plt.savefig(f'{case_id}_{flt}_local_peaks.png',dpi=300)
plt.show()


'''

elif case_id=='c05':
    t_start =  np.datetime64("2024-07-10T13:37:00")
    t_end   =  np.datetime64("2024-07-10T15:25:00")
    norm_v=norm[1]

elif case_id=='c06':
    t_start =  np.datetime64("2024-10-08T23:56:00")
    t_end   =  np.datetime64("2024-10-09T01:25:00")
    norm_v=norm[2]

elif case_id=='c07':
    t_start =  np.datetime64("2024-11-01T00:16:00")
    t_end   =  np.datetime64("2024-11-01T02:05:00")
    norm_v=norm[3]

elif case_id=='c08':
    t_start =  np.datetime64("2024-11-01T12:31:00")
    t_end   =  np.datetime64("2024-11-01T14:18:00")
    norm_v=norm[4]

elif case_id=='c09':
    t_start =  np.datetime64("2024-11-12T22:22:00")
    t_end   =  np.datetime64("2024-11-13T00:10:00")
    norm_v=norm[5]

elif case_id=='c10':
    t_start =  np.datetime64("2024-11-13T15:08:00")
    t_end   =  np.datetime64("2024-11-13T16:57:00")
    norm_v=norm[6]
'''