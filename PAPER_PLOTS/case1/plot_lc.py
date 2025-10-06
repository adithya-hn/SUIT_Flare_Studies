

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

pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 
data1=(np.loadtxt(f'csv_files/NB03_c1_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'
data2=(np.loadtxt(f'csv_files/NB08_c1_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
data3=(np.loadtxt(f'csv_files/NB04_c1_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
Solexs=(np.loadtxt(f'csv_files/AL1_SOLEXS_20240601_SDD2_L1_puc_tb_fit_results_TEMP_EM.txt',delimiter=' ',dtype='str')).transpose()
Helios=(np.load("csv_files/cdte_data_flare_1.npy", allow_pickle=True)).transpose()

m_cls=datetime.fromisoformat('2024-06-01T08:25:00.000')
m_cls_p=datetime.fromisoformat('2024-06-01T08:49:00.000')

time_array1=np.array(data1[0], dtype='datetime64')
time_array2=np.array(data2[0], dtype='datetime64')
time_array3=np.array(data3[0], dtype='datetime64')
date=str(time_array1[0])#[:10] #time_array1[0].strftime('%Y-%m-%d')

lc1_mean = np.array(data1[1],dtype=float)/np.array(data1[5],dtype=float) # total/area
qs1_mean = np.array(data1[3],dtype=float)/np.array(data1[6],dtype=float)
lc1_mean_er= np.array(data1[2],dtype=float)/np.array(data1[5],dtype=float)
qs1_mean_er= np.array(data1[4],dtype=float)/np.array(data1[6],dtype=float)
n_lc1_er = (lc1_mean / qs1_mean) * np.sqrt( (lc1_mean_er / lc1_mean)**2 + (qs1_mean_er / qs1_mean)**2 )

lc2_mean = np.array(data2[1],dtype=float)/np.array(data2[5],dtype=float) # total/area
qs2_mean = np.array(data2[3],dtype=float)/np.array(data2[6],dtype=float)
lc2_mean_er= np.array(data2[2],dtype=float)/np.array(data2[5],dtype=float)
qs2_mean_er= np.array(data2[4],dtype=float)/np.array(data2[6],dtype=float)
n_lc2_er = (lc2_mean / qs2_mean) * np.sqrt( (lc2_mean_er / lc2_mean)**2 + (qs2_mean_er / qs2_mean)**2 )



lc3_mean = np.array(data3[1],dtype=float)/np.array(data3[5],dtype=float) # total/area
qs3_mean = np.array(data3[3],dtype=float)/np.array(data3[6],dtype=float)

lc3_mean_er= np.array(data3[2],dtype=float)/np.array(data3[5],dtype=float)
qs3_mean_er= np.array(data3[4],dtype=float)/np.array(data3[6],dtype=float)


cdte=Helios[1][19:]+Helios[2][19:]
cdte1_er=np.sqrt(np.array(Helios[3][19:], dtype=np.float64))
cdte2_er=np.sqrt(np.array(Helios[4][19:], dtype=np.float64))
cdte_er=cdte1_er+cdte2_er

datetime_objects = pd.to_datetime(Helios[0])
helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]

#---------------------------------------------------------%%%%%%%SOLEXS%%%%%%------------------------------------------#

time_array4=[]

sl_temp=np.array(Solexs[1],dtype=float)   #[float(tp) for tp in Solexs[1]]
sl_temp_er=np.array(Solexs[2],dtype=float)#[float(tpe) for tpe in Solexs[2]]
sl_Em=np.array(Solexs[3],dtype=float)     #[float(em) for em in Solexs[3]]
sl_Em_er=np.array(Solexs[4],dtype=float)  #[float(eme) for eme in Solexs[4]]

sltime=np.array([float(tp) for tp in Solexs[0]])
base_time = datetime(2024, 6, 1, 7, 0, 0)  # Jun 1, 2024 07:00:00 UTC
time_seconds = sltime-sltime[0]  # Convert string times to float seconds
time_array4 = [base_time + timedelta(seconds=int(t)) for t in time_seconds]
qs_fig=plt.figure(figsize=(8,5))
mean_line=np.mean(qs2_mean)
plt.axhline(mean_line,color='r',linestyle='--',label=f'mean QS  {mean_line:.2f}')
plt.plot(time_array2,qs2_mean)
plt.savefig('Figures/qs_lc.png',dpi=300)
plt.legend()
plt.show()


#-------------------

#fig, axs = plt.subplots(5, 1, sharex=True, figsize=(12,10))
fig, axs = plt.subplots(5, 1, sharex=True, figsize=(10,14),
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
#axs1_=axs[1].twinx()
axs[0].errorbar(time_array1, lc1_mean,yerr=n_lc1_er,fmt='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5, label="SUIT Mg II k"); axs[0].legend()
axs[1].errorbar(time_array2, lc2_mean,yerr=n_lc2_er,fmt='black', marker="o",capsize=2,markersize=2,linewidth=0.5, label="SUIT Ca II H"); axs[1].legend()
#axs[1].errorbar(time_array3, float_array3/10e6,yerr=float_array_er3/10e6,fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label="SUIT Mg II h"); axs[1].legend()
axs[2].errorbar(helio_time_array[19:], cdte,yerr=cdte_er,fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label="HEL1OS (CdTe1+CdTe2)"); axs[2].legend()
axs[3].errorbar(time_array4[4:soLen-8],sl_Em[4:soLen-8],yerr=sl_Em_er[4:soLen-8],fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label='SoLEXS Emission Measure'); axs[3].legend()
axs[4].errorbar(time_array4[4:soLen-8],sl_temp[4:soLen-8],yerr=sl_temp_er[4:soLen-8],fmt='g-', marker="o",capsize=2,markersize=2,linewidth=0.5, label='SoLEXS Temperature'); axs[4].legend()


axs[0].set_ylabel('Mg II k counts ')
axs[1].set_ylabel('Ca II H counts')
#axs1_.set_ylabel('Mg II k counts (x$10^{6}$)')
axs[2].set_ylabel('HEL1OS counts/min')
axs[3].set_ylabel('EM(x$10^{43}cm^{-3}$)')
axs[4].set_ylabel('Temperature (MK)')

#axs[1].set_yscale('log')

axs[2].set_yscale('log')
axs[3].set_yscale('log')
#axs1_.ticklabel_format(style='plain', axis='y', scilimits=(0,0))

axs[-1].set_xlabel(f"Start Time ({date})") # Shared x-label

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)

# Add global title
fig.suptitle(f"Flare Light Curves ", fontsize=20, weight='bold')
# Adjust layout so title doesn’t overlap
plt.subplots_adjust(top=0.95)
# Add panel labels (a), b), c), ...)
panel_labels = ['a)', 'b)', 'c)', 'd)', 'e)']
for i, ax in enumerate(axs):
    ax.text(0.02, 0.9, panel_labels[i],
            transform=ax.transAxes,
            fontsize=14, fontweight='bold',
            va='top', ha='left')
plt.savefig('Figures/case1_lc.png',dpi=300)
plt.close()
