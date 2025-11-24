import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.dates as mdates
import pandas as pd
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
import numpy.ma as ma

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

# ----------------------------
# CONFIGURE
# ----------------------------
channel = '131'
rgn_lc = f'csv_files/aia131_region_lc.csv'
output_plot = f'helios_aia_{channel}_fd_roi_lc.eps'
Helios=(np.loadtxt("csv_files/helios_CdTe_c4.csv",skiprows=1,delimiter=',',dtype='str')).transpose()
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
int_roi = rgn_lc_data[:, 2].astype(float)/exposure #region 4
int_outside = data[:, 3].astype(float)

cdte1=np.array(Helios[1],dtype=float)
cdte_er=np.array(Helios[2],dtype=float)
helio_time_array=np.array(Helios[0],dtype='datetime64')

# Convert date strings to datetime
times = [datetime.strptime(d, "%Y-%m-%d %H:%M:%S.%f") if '.' in d else datetime.strptime(d, "%Y-%m-%dT%H:%M:%S") for d in dates_str]

hel1=60
hel2=-31

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
ax1.errorbar(helio_time_array[hel1:hel2],cdte1[hel1:hel2],yerr=cdte_er[hel1:hel2], fmt='ro-',capsize=2,markersize=2,linewidth=0.5,label="HEL1OS (CdTe1+CdTe2)",alpha=0.5)
ax1.set_ylabel('HEL1OS (counts/min)',fontsize=13)
ax1.set_yscale('log')
#ax4.spines.right.set_position(("axes", 1.2))

#------
# int_roi_ = np.where(int_roi>2*np.median(int_roi),2*np.median(int_roi),int_roi)#np.ma.masked_array(int_roi, mask=mask)
# peaks=np.array(argrelextrema(int_roi_,np.greater)[0])
# peaks_sort=np.argsort(int_roi_[peaks])[::-1]
# peaks=peaks[peaks_sort[:5]]

# dt=np.array(times)
# for i in range(len(peaks)):
#     plt.axvline(dt[[peaks[i]]], color='b',alpha=0.2)



date=times[0].strftime('%Y-%m-%d')
ax3.spines.right.set_position(("axes", 1.15))
#plt.gcf().autofmt_xdate()
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
plt.show()




