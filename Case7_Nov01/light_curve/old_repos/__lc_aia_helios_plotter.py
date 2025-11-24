import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.dates as mdates
import pandas as pd
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
import numpy.ma as ma
from astropy.convolution import Box1DKernel, convolve


from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

# ----------------------------
# CONFIGURE
# ----------------------------
channel = '131'
rgn_lc = f'csv_files/aia131_region_lc.csv'
output_plot = f'helios_aia_{channel}_fd_roi_lc.png'
Helios=(np.loadtxt("csv_files/helios_CdTe_c7.csv",skiprows=1,delimiter=',',dtype='str')).transpose()
csv_file = f'csv_files/{channel}_lc.csv'
# ----------------------------
# LOAD CSV
# ----------------------------
# Load without header
data = np.genfromtxt(csv_file, delimiter=',', dtype=str, skip_header=1)
rgn_lc_data=np.genfromtxt(rgn_lc, delimiter=',', dtype=str, skip_header=1)

# Extract columns
exposure=rgn_lc_data[:, 1].astype(float)
dates_str = data[:, 0]
int_full = data[:, 1].astype(float)
int_roi = rgn_lc_data[:, 3].astype(float)/exposure #region 4
int_outside = data[:, 3].astype(float)

cdte1=np.array(Helios[1],dtype=float)
cdte_er=np.array(Helios[2],dtype=float)
helio_time_array=np.array(Helios[0],dtype='datetime64')

# datetime_objects = pd.to_datetime(Helios[0])
# helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]

# Convert date strings to datetime
times = [datetime.strptime(d, "%Y-%m-%d %H:%M:%S.%f") if '.' in d else datetime.strptime(d, "%Y-%m-%dT%H:%M:%S") for d in dates_str]

# ----------------------------
# PLOT
# ----------------------------
fig=plt.figure(figsize=(10, 5))
ax1=fig.add_subplot(111)
ax2=ax1.twinx()
ax3=ax1.twinx()
#ax4=ax1.twinx()

ax3.plot(times, int_full, label='Full Disk', marker='o', markersize=1,color='orange',linestyle='-.')
ax2.plot(times, int_roi, label='ROI Only', marker='s', markersize=1,color='green',linestyle='--')
ax2.plot(np.nan, np.nan, label='AIA 131 peaks', markersize=1,color='b',alpha=0.2)
#ax2.plot(times, int_outside, label='Outside ROI', marker='x', markersize=1)
ax1.errorbar(helio_time_array,cdte1,yerr=cdte_er, fmt='ro-',capsize=2,markersize=2,linewidth=0.5,label="HEL1OS (CdTe1+CdTe2)",alpha=0.5)
ax1.set_ylabel('HEL1OS (counts/min)',fontsize=13)
ax1.set_yscale('log')
#ax4.spines.right.set_position(("axes", 1.2))

#------
from scipy.signal import savgol_filter
flux_smooth = savgol_filter(int_roi, window_length=15, polyorder=3)
int_roi_ =np.where(flux_smooth>2*np.median(flux_smooth),2*np.median(flux_smooth),flux_smooth)#np.ma.masked_array(int_roi, mask=mask)

#peaks=np.array(argrelextrema(int_roi_,np.greater)[0])
#peaks_sort=np.argsort(int_roi_[peaks])[::-1]

#peaks=peaks[peaks_sort[:8]]

peaks, _ = find_peaks(int_roi_, prominence=1, width=10)#, height=1.03*np.median(int_roi))

dt=np.array(times)

for i in range(len(peaks)):
    plt.axvline(dt[[peaks[i]]], color='b',alpha=0.2)


#ax2.axhline(3*np.median(int_roi))
date=times[0].strftime('%Y-%m-%d')
ax3.spines.right.set_position(("axes", 1.15))
plt.gcf().autofmt_xdate()
ax1.set_xlabel("Time (in UT)")
ax3.set_ylabel("FD Intensity (DN/s)")
ax2.set_ylabel("ROI Intensity (DN/s)")
plt.title(f"AIA {channel} Å Light Curves ({date})")
ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')
#plt.grid(True)
#plt.legend()
plt.tight_layout()
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.6, 0.5))

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(output_plot, dpi=300)
plt.close()

from astropy.convolution import Box1DKernel, convolve

smoothed=convolve(int_roi, kernel=Box1DKernel(100), boundary='extend')


spiks=np.where((int_roi-smoothed)<100000,0,(int_roi-smoothed))
peaks, _ = find_peaks(spiks ,distance=10)
# peaks=np.array(argrelextrema(spiks,np.greater)[0])
# peaks_sort=np.argsort(spiks[peaks])[::-1]
# peaks=peaks[peaks_sort[:10]]
#flux_smooth = savgol_filter(int_roi, window_length=15, polyorder=3)

plt.plot(dt,smoothed)
plt.plot(dt,int_roi)
plt.plot(dt,spiks)
for i in range(len(peaks)):
    plt.axvline(dt[[peaks[i]]], color='b',alpha=0.2)
plt.show()
