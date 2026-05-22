

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
scol =sns.color_palette("colorblind")

#palette = sns.color_palette("deep")

#----------------------Input-parameters------------------

C_n=5 #case number
data1 =(np.loadtxt(f'csv_files/c{C_n}_NB04_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'
Solexs=(np.loadtxt(f'csv_files/fit_results_AL1_SOLEXS_20240710_SDD2_L1_2407101330_2407101600_TEMP_EM.txt',skiprows=1,dtype='str')).transpose()
Helios=(np.loadtxt(f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/light_curve/csv_files/Helios_LightCurve.csv', delimiter=',',skiprows=1,dtype='str')).transpose()
spikes_nb3=(np.loadtxt(f'csv_files/Diff_img_data_NB04.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
goes = (np.loadtxt('csv_files/goes_xray_lightcurve.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 

peaks_pos=(np.loadtxt('csv_files/helios_peaks.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
peaks_dt =np.array(peaks_pos[0],dtype='datetime64')
suit_pks_pos=(np.loadtxt('csv_files/suit_diff_peaks.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
suit_pks_dt =np.array(suit_pks_pos[0],dtype='datetime64')
rgn_lc = f'csv_files/aia131_region_lc.csv'
csv_file = f'csv_files/{131}_lc.csv'
#-------------------------------------------------------
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True)

goes_dt=np.array(goes[0],dtype='datetime64')
xrs_a=np.array(goes[1],dtype=float)
xrs_b=np.array(goes[2],dtype=float)

exposure= np.array(data1[3],dtype=float)

time_array1=np.array(data1[0], dtype='datetime64')
time_array5=np.array(spikes_nb3[0], dtype='datetime64')
nb3_counts=np.array(spikes_nb3[3],dtype=float)
nb3_counts_er=np.sqrt(nb3_counts*(exposure/1000))/(exposure/1000)

date=str(time_array1[0])#[:10] #time_array1[0].strftime('%Y-%m-%d')


lc1_tot = np.array(data1[1],dtype=float)
lc1_mean_er= np.sqrt(lc1_tot*(exposure/1000))/(exposure/1000)


cdte=np.array(Helios[1],dtype=float)
cdte_er=np.array(np.sqrt(cdte*60)/60)
helio_time_array=np.array(Helios[0],dtype='datetime64[us]')
#--------------------------131--------------------

data = np.genfromtxt(csv_file, delimiter=',', dtype=str, skip_header=1)
rgn_lc_data=np.genfromtxt(rgn_lc, delimiter=',', dtype=str, skip_header=1)

# Extract columns
exposure=rgn_lc_data[:, 1].astype(float)
dates_str = data[:, 0]
int_full = data[:, 1].astype(float)
int_roi = rgn_lc_data[:, 2].astype(float)/exposure #region 1
times = [datetime.strptime(d, "%Y-%m-%d %H:%M:%S.%f") if '.' in d else datetime.strptime(d, "%Y-%m-%dT%H:%M:%S") for d in dates_str]


#---------------------------------------------------------%%%%%%%SOLEXS%%%%%%------------------------------------------#

time_array4=[]

sl_temp=[float(tp) for tp in Solexs[1]]
sl_temp_er=[float(tpe) for tpe in Solexs[2]]
sl_Em=[float(em) for em in Solexs[3]]
sl_Em_er=[float(eme) for eme in Solexs[4]]
slt=Solexs[0]
time_array4=[datetime.strptime(str(ts)[:19], "%Y-%m-%dT%H:%M:%S") for ts in slt]


fig, axs = plt.subplots(4, 1, sharex=True, figsize=(12,16), gridspec_kw={'hspace': 0})  # no vertical spacing 
plt.rcParams["font.size"]=22
plt.rcParams["axes.labelsize"]=22
plt.rcParams["xtick.labelsize"]=22
plt.rcParams["ytick.labelsize"]=22
plt.rcParams["legend.fontsize"]=18
plt.rcParams["figure.titlesize"]=22
plt.rcParams["axes.titlesize"]=22
for i in range(len(axs)):  # all but bottom panel
    axs[i].ticklabel_format(style='plain', axis='y', scilimits=(0,0))
    
    if i ==0:
        axs[i].tick_params(axis="x", which="both", bottom=True, top=True) 
    else:
        #axs[i].label_outer()                     # hide x-labels
        axs[i].tick_params(axis="x", which="both", bottom=True, top=False) 
    #axs[i].yaxis.offsetText.set_position((-0.04,-0.1))  # adjust X,Y offset
    axs[i].grid(True, which='major', linestyle='--', alpha=0.6)

soLen=len(time_array4)

# Sort and calculate time difference
time_array1_sec = np.array(data1[0], dtype='datetime64[s]')
dt = np.diff(time_array1_sec).astype('timedelta64[s]').astype(float)
gap_threshold = 120  # seconds
gap_indices = np.where(dt > gap_threshold)[0]
print("Gap indices:", gap_indices)
axs1=axs[0].twinx()
axs2=axs[2].twinx()
axs3=axs[3].twinx()
start = 0
for idx in gap_indices:
    axs[0].errorbar(time_array1[start:idx+1], (lc1_tot)[start:idx+1]/10e8,yerr=lc1_mean_er[start:idx+1]/10e5,fmt='k', marker="o",capsize=2,markersize=2,linewidth=0.5)
    axs1.errorbar(time_array5[start:idx+1],nb3_counts[start:idx+1],yerr=nb3_counts_er[start:idx+1]*1e3,color=scol[0], marker="o",capsize=2,markersize=2,linewidth=0.5)
    start = idx + 1
axs[0].errorbar(time_array1[start:], (lc1_tot)[start:]/10e8,yerr=lc1_mean_er[start:]/10e5,fmt='k', marker="o",capsize=2,markersize=2,linewidth=0.5, label=r"SUIT Mg II h (errors multiplied by $ 10^3$)")
axs1.errorbar(time_array5[start:],nb3_counts[start:],yerr=nb3_counts_er[start:]*1e3,color=scol[0], marker="o",capsize=2,markersize=2,linewidth=0.5,label=r'Excess intensity (errors multiplied by $ 10^3$)')
axs[1].errorbar(helio_time_array, cdte,yerr=cdte_er,color=scol[2], marker="o",capsize=2,markersize=2,linewidth=0.5, label="HEL1OS (CdTe1+CdTe2)"); axs[1].legend(loc='upper center')
axs[2].errorbar(time_array4,sl_temp,yerr=sl_temp_er,color=scol[3], marker="o",capsize=2,markersize=2,linewidth=0.5, label='SoLEXS Temperature'); axs[2].legend(loc='upper center')
axs2.plot(goes_dt,xrs_b/1e-6,color=scol[9],markersize=2,linewidth=1, label='GOES: 1–8 Å')
axs1.set_yscale('log')
axs[3].plot(times, np.log10(int_full), label='AIA 131 Full-disk', color='tab:green',linestyle='dotted',linewidth=2); axs[3].legend(loc='upper center')
axs3.plot(times, np.log10(int_roi), label='AIA 131 ROI only',color='k',linestyle='--',linewidth=2); axs[3].legend(loc='upper center')
# axs3.set_ylim(6.775,6.794)
print(time_array1[-1],time_array4[-1],time_array5[-1],times[-1])
for pk in peaks_dt:
    axs[0].axvline(pk,alpha=0.2,color='r')
    axs[1].axvline(pk,alpha=0.2,color='r')
    axs[2].axvline(pk,alpha=0.2,color='r')
    axs[3].axvline(pk,alpha=0.2,color='r')
for pk in suit_pks_dt:
    axs[0].axvline(pk,alpha=0.6,color='tab:purple')
axs[1].axvline(pk,alpha=0.0,color='r',label='HEL1OS peaks')

# ask matplotlib for the plotted objects and their labels
lines, labels = axs[0].get_legend_handles_labels()
lines2, labels2 = axs1.get_legend_handles_labels()
axs[0].legend(lines + lines2, labels + labels2, loc='upper center',bbox_to_anchor=(0.4, 0.99))

ln1, lbl1 = axs[2].get_legend_handles_labels()
ln2, lbl2 = axs2.get_legend_handles_labels()
axs[2].legend(ln1 + ln2, lbl1 + lbl2, loc='upper center')

ln3, lbl3 = axs[3].get_legend_handles_labels()
ln4, lbl4 = axs3.get_legend_handles_labels()
axs[3].legend(ln3 + ln4, lbl3 + lbl4, loc='upper center')

ul=np.power(10,7.7)
ll=np.power(10,5.4)
#axs1.set_ylim(ll,ul)
axs1.set_ylabel('Excess intensity (DN/s) ', fontsize=20,color=scol[0])
axs[0].set_ylabel(r'Mg II h counts ($\times 10^8$ DN/s) ',color='k')
axs[1].set_ylabel('HEL1OS counts/s',color=scol[2])
#axs[3].set_ylabel(r'EM($\mathrm{\times10^{43}cm^{-3}}$)')
axs2.set_ylabel(r'X-ray flux ($\times 10^-6$  W/m²)',color=scol[9])
axs[2].set_ylabel('Temperature (MK)',color=scol[3])
axs[3].set_ylabel(r"$\log ~\mathrm{FD~ intensity ~(DN/s)}$",color='tab:green')  
axs3.set_ylabel(r"$\log~ \mathrm{ROI~ intensity ~(DN/s)}$",color='k') 

axs2.ticklabel_format(style='plain', axis='y')
axs2.yaxis.get_major_formatter().set_scientific(False)
axs2.yaxis.get_offset_text().set_visible(False)
axs[1].set_yscale('log')
#axs2.set_yscale('log')
axs[-1].set_xlabel(f"Time (UT)") # Shared x-label

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)

# Add global title
fig.suptitle(f"Flare Light Curves (SOL2024-07-10T15:37) ", fontsize=20, weight='bold')
# Adjust layout so title doesn’t overlap
plt.subplots_adjust(top=0.95)
# Add panel labels (a), b), c), ...)
panel_labels = ['a)', 'b)', 'c)', 'd)']
for i, ax in enumerate(axs):
    ax.text(0.02, 0.9, panel_labels[i],
            transform=ax.transAxes,
            fontsize=20, fontweight='bold',
            va='top', ha='left')

plt.savefig(f'case{C_n}_lc.pdf',dpi=300)
plt.close()
