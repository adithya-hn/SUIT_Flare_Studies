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


scol =sns.color_palette("colorblind")

files = sorted(glob.glob("csv_files/*.csv"))
norm =np.loadtxt('norm_vals.csv',delimiter=',')[1]


all_times = []
all_ints = []

for f in files:
    data = np.genfromtxt(f, delimiter=',', names=True, dtype=None, encoding='utf-8')
    
    # Convert time to datetime64
    t = np.array(data['Date'], dtype='datetime64[ms]')
    area=np.array(data['area'],dtype=float)
    #x = np.array(data['diff_count'], dtype=float)#/area
    x = np.array(data['Img_count'], dtype=float)#/area
    
    case_id=os.path.basename(f)[0:3]
    flt=os.path.basename(f)[-8:-4]
    #print(flt)
    
    if case_id=='c04':
        t_start =  np.datetime64("2024-07-10T03:59:00")
        t_end   =  np.datetime64("2024-07-10T05:44:00")
        norm_v=norm[0]

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
    
    mask = (t >= t_start) & (t <= t_end)
    t_cut = t[mask]
    img_count_cut  = x[mask]/norm_v
    mid_cut=np.median(img_count_cut)

    fig,ax=plt.subplots(figsize=(14,6))

    threshold=1.0e5 #5661*4 #2.5e6 #
    peaks=find_peaks(img_count_cut, height=threshold)
    
    dt=np.array(t_cut[peaks[0]])
    pk=np.array(img_count_cut[peaks[0]])

    all_times.append(dt)
    all_ints.append(pk)
    


    df = pd.DataFrame({"Peak time": dt,  "Peak total count": pk})
    df.to_csv(f"pks/{case_id}_{flt}_peaks.csv", index=False)
    ax.plot(t_cut,img_count_cut)

    for i in range(len(peaks[0])):
        plt.axvline(t_cut[peaks[0][i]], color='b',alpha=0.2)
    #plt.yscale('log')
    plt.title(f'{case_id}')
    plt.axhline(mid_cut,color='r',linestyle='-',label='Median')
    plt.axhline(threshold,color='g', linestyle=':',label='Threshold')
    plt.legend()
    plt.savefig(f'pks/{case_id}_{flt}_local_peaks.png',dpi=300)
    plt.close()

all_c = np.concatenate(all_ints)
'''
print(x.min(),x.max())
bins=np.arange(x.min(),x.max(),10000000)

print(len(all_c))
plt.figure(figsize=(12,5))
counts, edges = np.histogram(x, bins=bins)
bin_centers = 0.5 * (edges[1:] + edges[:-1])
mask = counts > 0
y = counts[mask]
x_fit = bin_centers[mask]
logx = np.log10(x_fit)
logy = np.log10(y)

# slope, intercept = np.polyfit(logx, logy, 1)

# alpha = -slope
# A = 10**intercept

def expo_func(x,a,c):
    f=c*np.exp



# histogram
plt.hist(x, bins=bins, alpha=0.5, label='data')

# fitted curve
# x_plot = np.logspace(np.log10(min(x_fit)), np.log10(max(x_fit)), 200)
# y_plot = A * x_plot**(-alpha)

# plt.plot(x_plot, y_plot, label=f'Power law fit: alpha={alpha:.2f}')
plt.show()

'''
# from scipy.stats import linregress

# # intensities = your list of intensities

intensities=all_c
print('No of peaks: ',len(all_c))
# intensities = np.array(intensities)

# # log-spaced bins
# bins = np.logspace(np.log10(intensities.min()), np.log10(intensities.max()), 1000)

# hist, bin_edges = np.histogram(intensities, bins=bins)
# bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])

# # remove empty bins
# mask = hist > 0
# x1 = np.log10(bin_centers[mask])
# y1 = np.log10(hist[mask])

# slope, intercept, r, p, stderr = linregress(x1, y1)
# alpha = -slope

# # Plot
# plt.figure(figsize=(7, 5))

# plt.scatter(bin_centers[mask], hist[mask], label="Data (binned)")
# #plt.plot(bin_centers[mask],10**(intercept + slope * x),label=f"Fit: α = {alpha:.2f}",  linewidth=2)

# plt.xscale("log")
# plt.yscale("log")
# plt.xlabel("Flare Intensity")
# plt.ylabel("Number of Flares")
# plt.legend()
# plt.tight_layout()
# plt.show()

import numpy as np
import matplotlib.pyplot as plt

I = np.array(intensities)
I = np.sort(I)

# CCDF: P(I' >= I)
N = len(I)
ccdf = 1 - np.arange(1, N+1) / N

plt.loglog(I, ccdf, 'o', label='CCDF')

# Overplot fitted power law
alpha = 1 + N / np.sum(np.log(I / I[0]))
C = ccdf[0] * I[0]**(alpha-1)
plt.loglog(I, C * I**(-(alpha-1)), label=f"Fit α={alpha:.2f}")

plt.xlabel("Intensity")
plt.ylabel("P(I' ≥ I)")
plt.legend()
plt.tight_layout()
plt.savefig('powerlaw.png',dpi=300)
plt.show()


I = np.array(intensities)
I = I[I > 0]

Imin = I.min()   # or your completeness threshold

alpha = 1 + len(I) / np.sum(np.log(I / Imin))
print("Power-law index α =", alpha)
