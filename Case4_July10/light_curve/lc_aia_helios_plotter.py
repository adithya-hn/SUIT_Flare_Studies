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
from matplotlib.ticker import FixedLocator, FuncFormatter

scol =sns.color_palette("colorblind")
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
#set_pub_style()

# ----------------------------
# CONFIGURE
# ----------------------------
channel = '131'
rgn_lc = f'csv_files/aia131_region_lc.csv'
Helios=(np.loadtxt("csv_files/helios_CdTe_c4.csv",skiprows=1,delimiter=',',dtype='str')).transpose()
csv_file = f'csv_files/{channel}_lc.csv'
nb4_csv  ='csv_files/Diff_img_data_NB04.csv'
data1 =(np.loadtxt(f'csv_files/c4_NB04_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose()
output_plot = f'c4_helios_aia_{channel}_fd_roi_lc.pdf'

# ----------------------------
# Load without header
data = np.genfromtxt(csv_file, delimiter=',', dtype=str, skip_header=1)
rgn_lc_data=np.genfromtxt(rgn_lc, delimiter=',', dtype=str, skip_header=1)
spikes_nb4=(np.loadtxt(nb4_csv,delimiter=',',skiprows=1,dtype='str')).transpose() 

# Extract columns
exposure=rgn_lc_data[:, 1].astype(float)
dates_str = data[:, 0]
int_full = data[:, 1].astype(float)
int_roi = rgn_lc_data[:, 2].astype(float)/exposure #region 4
#int_outside = data[:, 3].astype(float)

cdte1=np.array(Helios[1],dtype=float)
cdte_er=np.array(Helios[2],dtype=float)
helio_time_array=np.array(Helios[0],dtype='datetime64')
times = [datetime.strptime(d, "%Y-%m-%d %H:%M:%S.%f") if '.' in d else datetime.strptime(d, "%Y-%m-%dT%H:%M:%S") for d in dates_str]

hel1=60
hel2=-31

suit_exp= np.array(data1[3],dtype=float)
time_array1=np.array(spikes_nb4[0], dtype='datetime64')
nb4_counts=np.array(spikes_nb4[3],dtype=float)
nb4_counts_er=np.sqrt(nb4_counts*(suit_exp/1000))/(suit_exp/1000)
time_array1_sec = np.array(spikes_nb4[0], dtype='datetime64[s]')
dt = np.diff(time_array1_sec).astype('timedelta64[s]').astype(float)
gap_threshold = 120  # seconds
gap_indices = np.where(dt > gap_threshold)[0]


# ----------------------------
# PLOT
# ----------------------------
fig=plt.figure(figsize=(18, 8))
plt.rcParams["font.size"]=30
plt.rcParams["axes.labelsize"]=30
plt.rcParams["xtick.labelsize"]=30
plt.rcParams["ytick.labelsize"]=30
plt.rcParams["legend.fontsize"]=21
plt.rcParams["figure.titlesize"]=30
plt.rcParams["axes.titlesize"]=30
ax1=fig.add_subplot(111)
ax2=ax1.twinx()
ax3=ax1.twinx()
ax4=ax1.twinx()

ax3.plot(times, np.log10(int_full), label='AIA 131 Full-disk', color=scol[3],linestyle='dotted',linewidth=2)
#ax2.plot(np.nan, np.nan, label='AIA 131 peaks', markersize=1,color='b',alpha=0.2)
ax2.plot(times, np.log10(int_roi), label='AIA 131 ROI only',color='k',linestyle='--',linewidth=2)
start=0
for idx in gap_indices:
    ax4.errorbar(time_array1[start:idx+1], nb4_counts[start:idx+1],yerr=nb4_counts_er[start:idx+1],color=scol[0], marker="o",capsize=2,markersize=1,linewidth=1)
    start=idx+1
ax4.errorbar(time_array1[start:], nb4_counts[start:],yerr=nb4_counts_er[start:],marker="o",capsize=2,markersize=1,linewidth=1, label='SUIT difference image intensity',color=scol[0])
ax1.errorbar(helio_time_array,cdte1,yerr=cdte_er, color=scol[2],marker='o',capsize=2,markersize=2,linewidth=0.5,label="HEL1OS (CdTe1+CdTe2)",alpha=0.5)
ax2.ticklabel_format(style='plain', axis='y')


#ax4.spines.right.set_position(("axes", 1.2))

#------
# int_roi_ = np.where(int_roi>2*np.median(int_roi),2*np.median(int_roi),int_roi)#np.ma.masked_array(int_roi, mask=mask)
# peaks=np.array(argrelextrema(int_roi_,np.greater)[0])
# peaks_sort=np.argsort(int_roi_[peaks])[::-1]
# peaks=peaks[peaks_sort[:5]]

# dt=np.array(times)
# for i in range(len(peaks)):
#     plt.axvline(dt[[peaks[i]]], color='b',alpha=0.2)

# def fmt(y, pos):
#     return f"{y:.0f}"
# ax2.yaxis.set_major_formatter(FuncFormatter(fmt))
# ax3.yaxis.set_major_formatter(FuncFormatter(fmt))


date=times[0].strftime('%Y-%m-%d')
plt.title(f"AIA {channel} Å Light Curves ({date})")

ax3.spines.right.set_position(("axes", 1.12))
ax4.spines.right.set_position(("axes", 1.24))
ax1.set_xlabel("Time (UT)")
ax3.set_ylabel(r"$\log ~\mathrm{FD~ Intensity ~(DN/s)}$",color=scol[3])
ax3.tick_params(axis='y', colors=scol[3])
ax3.tick_params(axis='y', which='minor', color=scol[3])
ax3.spines['right'].set_color(scol[3])

ax2.set_ylabel(r"$\log~ \mathrm{ROI~ intensity ~(DN/s)}$",color='k')    

ax4.set_ylabel("Difference image intensity (DN/s)",color=scol[0])
ax4.tick_params(axis='y', colors=scol[0])
ax4.tick_params(axis='y', which='minor', color=scol[0])
ax4.spines['right'].set_color(scol[0])
ax1.set_ylabel('HEL1OS (counts/min)',color=scol[2])
ax1.set_yscale('log')
ax4.set_yscale('log')

# ticks = range(2, 14, 4)   # 2,4,6,8,10,12
# ax2.yaxis.set_major_locator(FixedLocator(ticks))
# ax2.yaxis.set_major_formatter(FuncFormatter(lambda val, pos: f"{val:g}"))


plt.tight_layout()
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.46, 0.5155))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(output_plot, dpi=300)
plt.show()




