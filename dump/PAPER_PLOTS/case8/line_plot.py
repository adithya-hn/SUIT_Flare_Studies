
#-----------------------------------------------------------------------

#Date created: 2025-07-10
#Author: @adithya-hn
#Description: This script plots the AIA full-disk and region-of-interest (ROI) light curves from a CSV file.
#---------Removes short exposure images---------------------------------


#-----------------------------------------------------------------------

import timeit
from scipy import stats as S
import scipy as sp
import pathlib
import pandas as pd
from subprocess import call
from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import *
from astropy.io import fits
import math as mt
from datetime import datetime
import matplotlib.dates as mdates
import os

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

# ----------------------------
# CONFIGURE
# ----------------------------
channel = '131'
csv_file = f'csv_files/aia_{channel}_all_lc.csv'
output_plot = f'aia_{channel}_fd_roi_lc.png'
            #csv_files/aia_131_all_lc.csv

# ----------------------------
# LOAD CSV
# ----------------------------
# Load without header
data = np.genfromtxt(csv_file, delimiter=',', dtype=str, skip_header=1)

# Extract columns
dates_str = data[:, 0]
int_full = data[:, 1].astype(float)
int_roi = data[:, 2].astype(float)
int_outside = data[:, 3].astype(float)
exposure = np.array(data[:,4].astype(float))

# Convert date strings to datetime
times = np.array([datetime.strptime(d, "%Y-%m-%dT%H:%M:%S.%f") if '.' in d else datetime.strptime(d, "%Y-%m-%dT%H:%M:%S") for d in dates_str])



# ----------------------------
# PLOT
# ----------------------------

fig=plt.figure(figsize=(11, 5))
ax1=fig.add_subplot(111)
ax2=ax1.twinx()
ax3=ax1.twinx()

clip_val=2.9

times=times[exposure>2.9]
int_full=int_full[exposure>2.9]
int_roi=int_roi[exposure>2.9]
int_outside=int_outside[exposure>2.9]


ax1.plot(times, int_full, label='Full Disk', marker='o', markersize=1,color='orange')
ax2.plot(times, int_roi, label='ROI Only', marker='s', markersize=1,color='green')
ax3.plot(times, int_outside, label='Outside ROI', marker='x', markersize=1)


date=times[0].strftime('%Y-%m-%d')
ax3.spines.right.set_position(("axes", 1.12))
plt.gcf().autofmt_xdate()
ax1.set_xlabel("Time")
ax1.set_ylabel("FD Intensity (DN/s)")
ax3.set_ylabel("Intensity Outside ROI (DN/s)")
ax2.set_ylabel("ROI Intensity (DN/s)")
plt.title(f"AIA {channel} Å Light Curves ({date})")
ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')

#plt.grid(True)
#plt.legend()
plt.tight_layout()
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(output_plot, dpi=300)
plt.close()