import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.dates as mdates
import pandas as pd
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
import numpy.ma as ma
import seaborn as sns
scol =sns.color_palette("colorblind")


data=(np.loadtxt('Diff_img_data_NB04.csv',delimiter=',',skiprows=1,dtype='str')).transpose()

t_start =  np.datetime64("2024-11-01T12:31:00")
t_end   =  np.datetime64("2024-11-01T14:18:00")

dt=np.array(data[0], dtype='datetime64')
spikes=np.array(data[3],dtype=float)
img_count=np.array(data[5],dtype=float)
mask = (dt >= t_start) & (dt <= t_end)
dt_cut = dt[mask]
img_count_cut  = img_count[mask]
mid_cut=np.median(img_count_cut)

fig,ax=plt.subplots(figsize=(14,6))
peaks=np.array(argrelextrema(img_count_cut,np.greater)[0])

t=np.array(dt_cut[peaks])
pk=np.array(img_count_cut[peaks])

#np.savetxt(f'c6_lc_data.csv',np.c_[date_array, fltr_count,smoothed_count,area,exposure_times,meas_exp_times],delimiter=',',header='Time,Total,Area,Exposure,Meas_exposure',comments='' ,fmt='%s')
print(t[0],pk[0])
df = pd.DataFrame({"Peak time": t,  "Peak total count": pk})
df.to_csv("c8_peaks.csv", index=False)
ax.plot(dt_cut,img_count_cut)

for i in range(len(peaks)):
    plt.axvline(dt_cut[[peaks[i]]], color='b',alpha=0.2)
#plt.yscale('log')
plt.axhline(mid_cut)
plt.savefig('local_peaks.png',dpi=300)
plt.show()
