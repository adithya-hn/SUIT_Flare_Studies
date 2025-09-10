

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os
import pandas as pd
import matplotlib.dates as mdates
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
#set_pub_style()

#-----------------------------------

data_file='NB03_threshold_count.csv'
time_col=0
plot_col=2
er_col=3
param1=1
out_file='Threshold_plot.png'
plot_title='Threshold light curve'
plt_x_label='Time'
plt_y_label= 'Threshold_total_count'
plot_lable='Mg II k'
flare_start='2024-06-02T08:50:00.000'
flare_peak ='2024-06-02T08:56:00.000'

#-----------------------------------
data_file2='cdte_data_flare_3.npy'
time_col2=0
plot_col2=1
plot_col3=2
er_col2=3
er_col2=4

helios=np.load(data_file2, allow_pickle=True).transpose()
helios_time_array = [datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in helios[time_col2]]
helios_float_array1 = helios[plot_col2]+ helios[plot_col3]

#-----------------------------------


data=(np.loadtxt(data_file,delimiter=',',dtype='str')).transpose() 
date_array=data[time_col] 
area_array=np.array(data[param1],dtype=float)

time_array=[]
print(len(date_array))
for i in range(len(date_array)):
    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time)


fig,axs=plt.subplots(1,1, figsize=(10,5))

ax2=axs.twinx()
float_array = np.array(data[plot_col],dtype=float) 
float_array_er = np.array(data[er_col],dtype=float)
y_er=np.std(float_array_er)

axs.plot(helios_time_array, helios_float_array1, color='tab:orange', marker="o", markersize=2, linewidth=0.5, label='HEL1OS (CdTe 10-40 keV) counts')
ax2.plot(time_array,float_array/area_array,markersize=2,linewidth=0.5,label=plot_lable)
#ax2.plot(time_array,list(map(int,float_array_er)),color='k',markersize=2,linewidth=0.5,label=plot_lable)
#axs.errorbar(time_array,list(map(int,float_array)),yerr=y_er,fmt='ko',capsize=2,markersize=2,linewidth=0.5,label=Filters[0])

axs.set_yscale('log')
m_cls=datetime.fromisoformat(flare_start)
m_cls_p=datetime.fromisoformat(flare_peak)


axis_title='Total count'
img_nm=out_file

plt.ylabel(axis_title,fontsize=13)
plt.xlabel('Time',fontsize=13)
plt.axvline(m_cls_p,color='orange',linestyle='-',label='GOES Flare peak time')
plt.axvline(m_cls,color='orange',linestyle='--',label='GOES Flare start time')
plt.title('Mg II h Light Curve and HELIOS')
plt.legend(loc='best')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.show()
