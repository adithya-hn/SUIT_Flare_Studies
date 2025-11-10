

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

C_n=8 #case number
data1=(np.loadtxt(f'csv_files/NB04_c{C_n}_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() #'NB03_Light_curve_data.dat'
#data2=(np.loadtxt(f'csv_files/NB08_c{C_n}_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
#data3=(np.loadtxt(f'csv_files/NB04_c{C_n}_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
Solexs=(np.loadtxt(f'csv_files/fit_results_AL1_SOLEXS_20241101_SDD2_L1_2411011100_2411011500_TEMP_EM.txt',skiprows=1,dtype='str')).transpose()
Helios=(np.load(f"csv_files/cdte_data_flare_{C_n}.npy", allow_pickle=True)).transpose()
spikes_nb3=(np.loadtxt(f'csv_files/NB04_brightenings.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
spikes_nb8=(np.loadtxt(f'csv_files/NB08_brightenings.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 



m_cls=datetime.fromisoformat('2024-06-02T04:41:00.000')
m_cls_p=datetime.fromisoformat('2024-06-02T04:50:00.000')

nb3_c1=36 #Data gap indicies
nb8_c1=37
nb8_c2=88

#-------------------------------------------------------
pathlib.Path("Figures").mkdir(parents=True, exist_ok=True) 

time_array1=np.array(data1[0], dtype='datetime64')
#time_array2=np.array(data2[0], dtype='datetime64')
#time_array3=np.array(data3[0], dtype='datetime64')

time_array5=np.array(spikes_nb3[0], dtype='datetime64')
time_array6=np.array(spikes_nb8[0], dtype='datetime64')
nb3_counts=np.array(spikes_nb3[1],dtype=float)
nb8_counts=np.array(spikes_nb8[1],dtype=float)




date=str(time_array1[0])#[:10] #time_array1[0].strftime('%Y-%m-%d')


lc1_mean = np.array(data1[1],dtype=float)
lc1_mean_er= np.array(data1[2],dtype=float)

# lc2_mean = np.array(data2[1],dtype=float)
# lc2_mean_er= np.array(data2[2],dtype=float)

cdte=Helios[1]+Helios[2]
cdte1_er=(np.array(Helios[1], dtype=np.float64))**2
cdte2_er=(np.array(Helios[2], dtype=np.float64))**2
cdte_er=np.sqrt(cdte1_er+cdte2_er)

datetime_objects = pd.to_datetime(Helios[0])
helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]

#---------------------------------------------------------%%%%%%%SOLEXS%%%%%%------------------------------------------#

time_array4=[]

'''sl_temp=np.array(Solexs[1],dtype=float)   #[float(tp) for tp in Solexs[1]]
sl_temp_er=np.array(Solexs[2],dtype=float)#[float(tpe) for tpe in Solexs[2]]
sl_Em=np.array(Solexs[3],dtype=float)     #[float(em) for em in Solexs[3]]
sl_Em_er=np.array(Solexs[4],dtype=float)  #[float(eme) for eme in Solexs[4]]
'''
sl_temp=[float(tp) for tp in Solexs[1]]
sl_temp_er=[float(tpe) for tpe in Solexs[2]]
sl_Em=[float(em) for em in Solexs[3]]
sl_Em_er=[float(eme) for eme in Solexs[4]]
slt=Solexs[0]
#sltime=np.array([float(tp) for tp in Solexs[0]])
time_array4=[datetime.strptime(str(ts)[:19], "%Y-%m-%dT%H:%M:%S") for ts in slt]

'''base_time = datetime(2024, 6, 1, 7, 0, 0)  # Jun 1, 2024 07:00:00 UTC
time_seconds = sltime-sltime[0]  # Convert string times to float seconds
time_array4 = [base_time + timedelta(seconds=int(t)) for t in time_seconds]

'''


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
hel1=75
sole1=150
#axs1_=axs[1].twinx()
axs[0].errorbar(time_array1[:nb3_c1], (lc1_mean)[:nb3_c1],yerr=lc1_mean_er[:nb3_c1],fmt='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5, label="SUIT Mg II h"); axs[0].legend(loc='lower right')
axs[0].errorbar(time_array1[nb3_c1:], (lc1_mean)[nb3_c1:],yerr=lc1_mean_er[nb3_c1:],fmt='tab:blue', marker="o",capsize=2,markersize=2,linewidth=0.5); axs[0].legend(loc='lower right')
#axs[1].errorbar(time_array2[:nb8_c1], (lc2_mean)[:nb8_c1],yerr=lc2_mean_er[:nb8_c1],fmt='black', marker="o",capsize=2,markersize=2,linewidth=0.5, label="SUIT Ca II H"); axs[1].legend(loc='lower right')
#axs[1].errorbar(time_array2[nb8_c1:nb8_c2], (lc2_mean)[nb8_c1:nb8_c2],yerr=lc2_mean_er[nb8_c1:nb8_c2],fmt='black', marker="o",capsize=2,markersize=2,linewidth=0.5); axs[1].legend(loc='lower right')
#axs[1].errorbar(time_array2[nb8_c2:], (lc2_mean)[nb8_c2:],yerr=lc2_mean_er[nb8_c2:],fmt='black', marker="o",capsize=2,markersize=2,linewidth=0.5); axs[1].legend(loc='lower right')
#axs[1].errorbar(time_array3, float_array3/10e6,yerr=float_array_er3/10e6,fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label="SUIT Mg II h"); axs[1].legend()
axs[2].errorbar(helio_time_array[hel1:], cdte[hel1:],yerr=cdte_er[hel1:],fmt='tab:red', marker="o",capsize=2,markersize=2,linewidth=0.5, label="HEL1OS (CdTe1+CdTe2)"); axs[2].legend(loc='lower right')
axs[4].errorbar(time_array4[sole1:],sl_Em[sole1:],yerr=sl_Em_er[sole1:],fmt='gray', marker="o",capsize=2,markersize=2,linewidth=0.5, label='SoLEXS Emission Measure'); axs[4].legend(loc='lower right')
axs[3].errorbar(time_array4[sole1:],sl_temp[sole1:],yerr=sl_temp_er[sole1:],fmt='g-', marker="o",capsize=2,markersize=2,linewidth=0.5, label='SoLEXS Temperature'); axs[3].legend(loc='lower right')
axs1=axs[0].twinx()
axs2=axs[1].twinx()
axs1.plot(time_array5[:nb3_c1],nb3_counts[:nb3_c1],'y')
axs1.plot(time_array5[nb3_c1:],nb3_counts[nb3_c1:],'y')
axs2.plot(time_array6[:nb8_c1],nb8_counts[:nb8_c1],'y')
axs2.plot(time_array6[nb8_c1:nb8_c2],nb8_counts[nb8_c1:nb8_c2],'y')
axs2.plot(time_array6[nb8_c2:],nb8_counts[nb8_c2:],'y')


axs[0].set_ylabel('Mg II k counts ')
axs[1].set_ylabel('Ca II H counts')
#axs1_.set_ylabel('Mg II k counts (x$10^{6}$)')
axs[2].set_ylabel('HEL1OS counts/min')
axs[4].set_ylabel('EM(x$10^{43}cm^{-3}$)')
axs[3].set_ylabel('Temperature (MK)')

#axs[1].set_yscale('log')

axs[2].set_yscale('log')
axs[4].set_yscale('log')
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
plt.savefig(f'case{C_n}_lc.png',dpi=300)
plt.close()
