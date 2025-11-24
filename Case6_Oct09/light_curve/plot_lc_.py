

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta                                                                 
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

#----------------------Input-parameters------------------

C_n=6 #case number
data1=(np.loadtxt(f'csv_files/c{C_n}_NB04_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'
#data2=(np.loadtxt(f'csv_files/NB08_c{C_n}_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
#data3=(np.loadtxt(f'csv_files/NB04_c{C_n}_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
Solexs=(np.loadtxt(f'csv_files/AL1_SOLEXS_20241009_SDD2_L1_puc_tb_fit_results_TEMP_EM.txt',skiprows=1,dtype='str')).transpose()
Helios=(np.load(f"csv_files/cdte_data_flare_{C_n}.npy", allow_pickle=True)).transpose()
spikes_nb3=(np.loadtxt(f'csv_files/Diff_img_data_NB04.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
goes = (np.loadtxt('csv_files/goes_xray_lightcurve.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 

# Plotting
# plt.figure(figsize=(10, 6))
# plt.plot(df.index, df['xrsa'], label='0.5–4 Å (short)', color='red')
# plt.plot(df.index, df['xrsb'], label='1–8 Å (long)', color='blue')

goes_dt=np.array(goes[0],dtype='datetime64')
xrs_a=np.array(goes[1],dtype=float)
xrs_b=np.array(goes[2],dtype=float)


# plt.yscale('log')
# plt.xlabel('Time (UTC)')
# plt.ylabel('Flux (W/m²)')
# plt.title('GOES Soft X-ray Light Curve')

# plt.grid(True)
# plt.tight_layout()

nb3_c1=36 #Data gap indicies
nb8_c1=37
nb8_c2=88

#-------------------------------------------------------
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 

time_array1=np.array(data1[0], dtype='datetime64')
time_array5=np.array(spikes_nb3[0], dtype='datetime64')
nb3_counts=np.array(spikes_nb3[3],dtype=float)

#time_array2=np.array(data2[0], dtype='datetime64')
#time_array3=np.array(data3[0], dtype='datetime64')
#time_array6=np.array(spikes_nb8[0], dtype='datetime64')
#nb8_counts=np.array(spikes_nb8[1],dtype=float)




date=str(time_array1[0])#[:10] #time_array1[0].strftime('%Y-%m-%d')


lc1_mean = np.array(data1[1],dtype=float)
lc1_mean_er= np.array(data1[2],dtype=float)

# lc2_mean = np.array(data2[1],dtype=float)
# lc2_mean_er= np.array(data2[2],dtype=float)

cdte=Helios[1]+Helios[2]
cdte1_er=(np.array(Helios[3], dtype=np.float64))**2
cdte2_er=(np.array(Helios[4], dtype=np.float64))**2
cdte_er=np.sqrt(cdte1_er+cdte2_er)

datetime_objects = pd.to_datetime(Helios[0])
helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]

#---------------------------------------------------------%%%%%%%SOLEXS%%%%%%------------------------------------------#

time_array4=[]

'''
sl_temp=np.array(Solexs[1],dtype=float)   #[float(tp) for tp in Solexs[1]]
sl_temp_er=np.array(Solexs[2],dtype=float)#[float(tpe) for tpe in Solexs[2]]
sl_Em=np.array(Solexs[3],dtype=float)     #[float(em) for em in Solexs[3]]
sl_Em_er=np.array(Solexs[4],dtype=float)  #[float(eme) for eme in Solexs[4]]
'''
sl_temp=[float(tp) for tp in Solexs[1]]
sl_temp_er=[float(tpe) for tpe in Solexs[2]]
sl_Em=[float(em) for em in Solexs[3]]
sl_Em_er=[float(eme) for eme in Solexs[4]]
slt=Solexs[0]
sltime=np.array([float(tp) for tp in Solexs[0]])
base_time = datetime(2024, 10, 9, 1, 0, 0)  # Jun 1, 2025 07:00:00 UTC
slt = sltime-sltime[0]
time_seconds = sltime-sltime[0]  # Convert string times to float seconds
# Convert seconds to datetime
time_array4 = [base_time + timedelta(seconds=int(t)) for t in time_seconds]



#fig, axs = plt.subplots(5, 1, sharex=True, figsize=(12,10))
fig, axs = plt.subplots(4, 1, sharex=True, figsize=(10,12),
                         gridspec_kw={'hspace': 0})  # no vertical spacing 
for i in range(len(axs)):  # all but bottom panel
    axs[i].ticklabel_format(style='plain', axis='y', scilimits=(0,0))
    #axs[i].yaxis.get_offset_text().set_fontsize(8)  # smaller font
    #axs[i].yaxis.get_offset_text().set_y(-0.35)
    
    if i ==0:
        axs[i].tick_params(axis="x", which="both", bottom=True, top=True) 
    else:
        #axs[i].label_outer()                     # hide x-labels
        axs[i].tick_params(axis="x", which="both", bottom=True, top=False) 
    #axs[i].yaxis.offsetText.set_position((-0.04,-0.1))  # adjust X,Y offset
    axs[i].grid(True, which='major', linestyle='--', alpha=0.6)

soLen=len(time_array4)
hel1=25
hel2=-28
sole1=0 
sole2=-75

# Sort and calculate time difference
time_array1_sec = np.array(data1[0], dtype='datetime64[s]')
dt = np.diff(time_array1_sec).astype('timedelta64[s]').astype(float)
gap_threshold = 120  # seconds
gap_indices = np.where(dt > gap_threshold)[0]
print("Gap indices:", gap_indices)
axs1=axs[0].twinx()
start = 0
for idx in gap_indices:
    axs[0].errorbar(time_array1[start:idx+1], (lc1_mean)[start:idx+1]/10e8,yerr=lc1_mean_er[start:idx+1]/10e8,fmt='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5)
    axs1.plot(time_array5[start:idx+1],nb3_counts[start:idx+1],'y')
    start = idx + 1
axs[0].errorbar(time_array1[start:], (lc1_mean)[start:]/10e8,yerr=lc1_mean_er[start:]/10e8,fmt='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5, label="SUIT Mg II h")
axs1.plot(time_array5[start:],nb3_counts[start:],'y',label=r'Difference image intensity above 4 $\sigma$')

axs[1].errorbar(helio_time_array[hel1:hel2], cdte[hel1:hel2],yerr=cdte_er[hel1:hel2],fmt='tab:red', marker="o",capsize=2,markersize=2,linewidth=0.5, label="HEL1OS (CdTe1+CdTe2)"); axs[1].legend(loc='upper center')
axs[2].errorbar(time_array4[sole1:sole2],sl_temp[sole1:sole2],yerr=sl_temp_er[sole1:sole2],fmt='g-', marker="o",capsize=2,markersize=2,linewidth=0.5, label='SoLEXS Temperature'); axs[2].legend(loc='upper center')
#axs[3].errorbar(time_array4[sole1:sole2],sl_Em[sole1:sole2],yerr=sl_Em_er[sole1:sole2],fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label='SoLEXS Emission Measure'); axs[3].legend(loc='upper center')
axs[3].plot(goes_dt,xrs_a,color='b',markersize=2,linewidth=1, label='GOES: 0.5–4 Å'); axs[3].legend(loc='center')
axs[3].plot(goes_dt,xrs_b,color='r',markersize=2,linewidth=1, label='GOES: 1–8 Å'); axs[3].legend(loc='center')



axs1.set_yscale('log')

# ask matplotlib for the plotted objects and their labels
lines, labels = axs[0].get_legend_handles_labels()
lines2, labels2 = axs1.get_legend_handles_labels()
axs[0].legend(lines + lines2, labels + labels2, loc='upper center')

#axs[0].set_ylim(1.865,1.933)
ul=np.power(10,7.7)
ll=np.power(10,5.4)
#axs1.set_ylim(ll,ul)
axs1.set_ylabel('Difference image counts ')
axs[0].set_ylabel(r'Mg II k counts ($\times 10^8$ DN) ')
axs[1].set_ylabel('HEL1OS counts/min')
#axs[3].set_ylabel(r'EM($\mathrm{\times10^{43}cm^{-3}}$)')
axs[3].set_ylabel('X-ray flux (W/m²)')
axs[2].set_ylabel('Temperature (MK)')

axs[1].set_yscale('log')
axs[3].set_yscale('log')
axs[-1].set_xlabel(f"Time (in UTC)") # Shared x-label

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)

# Add global title
fig.suptitle(f"Flare Light Curves ({date[:10]}) ", fontsize=20, weight='bold')
# Adjust layout so title doesn’t overlap
plt.subplots_adjust(top=0.95)
# Add panel labels (a), b), c), ...)
panel_labels = ['a)', 'b)', 'c)', 'd)']
for i, ax in enumerate(axs):
    ax.text(0.02, 0.9, panel_labels[i],
            transform=ax.transAxes,
            fontsize=16, fontweight='bold',
            va='top', ha='left')

plt.savefig(f'case{C_n}_lc.eps',dpi=300)
plt.close()


