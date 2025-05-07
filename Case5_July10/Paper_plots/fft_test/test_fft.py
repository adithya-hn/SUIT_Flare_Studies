import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.signal import detrend
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
import matplotlib.dates as mdates
from scipy.interpolate import interp1d

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()


colors_list = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e', '#000000', '#7f7f7f']
Filters=['NB08']
param1='magnetogram'
param=Filters[0]
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data=(np.loadtxt(f'core_nb04_lc_data.csv',delimiter=',',dtype='str')).transpose() #'NB08_Light_curve_data.dat'
NB3_data=(np.loadtxt(f'core_nb03_lc_data.csv',delimiter=',',dtype='str')).transpose()

date_array=data[0]
NB3_date_array=NB3_data[0]



time_array=[]
t1_stamps=[]
t2_stamps=[]
#print(len(data[1]))
for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    t1_stamps.append(parsed_time.timestamp())
    time_array.append(parsed_time)

nb3_time_array=[]
for i in range(len(NB3_date_array)):
    parsed_time = datetime.fromisoformat(NB3_date_array[i])
    t2_stamps.append(parsed_time.timestamp())
    nb3_time_array.append(parsed_time)
#___________________________________________________________

# Interpolating Light Curve 2 to match timestamps of Light Curve 1
print('-------->',len(data[0]),len(data[1]))
interp_func = interp1d(t2_stamps,NB3_data[1] , kind='linear',fill_value="extrapolate")
interp_func_qs = interp1d(t2_stamps,NB3_data[2] , kind='linear',fill_value="extrapolate")
counts_2_interp = interp_func(t1_stamps)
qs_intrep_counts=interp_func_qs(t1_stamps)


#fig.subplots_adjust(right=0.83)
#ax2 = axs.twinx()

float_array = np.array([float(string) for string in data[1]])
float_array_er =np.array([float(string) for string in data[2]])


nb3float_array = np.array([float(string) for string in NB3_data[1]])
nb3float_array_er =np.array( [float(string) for string in NB3_data[2]])

qs_sigma_ca=np.sqrt(float_array_er/20736)
qs_sigma_mg=np.sqrt(qs_intrep_counts/20736)
ar_sigma_ca=np.sqrt(float_array/20736)
ar_sigma_mg=np.sqrt(counts_2_interp/20736)

qs_sigma=(float_array_er/qs_intrep_counts)*np.sqrt((qs_sigma_ca / float_array_er)**2 + (qs_sigma_mg / qs_intrep_counts)**2)
ar_sigma=(float_array/counts_2_interp)*np.sqrt((ar_sigma_ca / float_array)**2 + (ar_sigma_mg / counts_2_interp)**2)

counts_2_interp=np.array(counts_2_interp)
float_array=np.array(float_array)
print(float_array_er[0],qs_intrep_counts[0],float_array_er[0]/qs_intrep_counts[0])
#axs.errorbar(time_array,float_array_er/qs_intrep_counts,yerr=qs_sigma,color=colors_list[1],marker='^',markersize=2,linewidth=0.5,label='QS ratio')
#axs.errorbar(time_array,(float_array/counts_2_interp),yerr=ar_sigma,color=colors_list[2],marker='s',markersize=2,linewidth=0.5,label='AR ratio ')

# Load your light curve: Replace with your actual file
# Assuming 'time' and 'intensity' are numpy arrays


time = t1_stamps #np.loadtxt('your_lc.csv', delimiter=',', usecols=0)
intensity = (float_array/counts_2_interp)#np.loadtxt('your_lc.csv', delimiter=',', usecols=1)

# Detrend the signal to focus on periodic variations
detrended = detrend(intensity)

# Sampling parameters
dt = time[1] - time[0]  # time step
N = len(detrended)

# FFT
yf = fft(detrended)
xf = fftfreq(N, dt)[:N//2]  # positive frequencies

# Power spectrum
power = np.abs(yf[:N//2])**2

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(xf, power, color='red', label='Power Spectrum')
#plt.plot(xf, yf[:N//2], color='blue', label='FFT')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power')
plt.title('FFT Power Spectrum')
plt.grid(True)
plt.savefig(f'Figures/FFT_Power_Spectrum_{param1}_{param}.png', dpi=300)
plt.show()
