import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.dates as mdates
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

data=(np.loadtxt('hmi_lc.csv',delimiter=',',skiprows=1,dtype='str')).transpose()
suit=(np.loadtxt('Diff_img_data_NB04.csv',delimiter=',',skiprows=1,dtype='str')).transpose()
suit_dt=np.array(suit[0],dtype='datetime64')
hmi_dt=np.array(data[0],dtype='datetime64')
suit_data=np.array(suit[1:,:],dtype=float)
nb4_bright_ct=suit_data[2]
bright_ct_er=suit_data[3]

sns_cl3=sns.color_palette('bright')
fig,ax=plt.subplots(figsize=(16,10))

print(data.shape)
data = np.array(data[1:,:], dtype=float)
print(data.shape)
#ax.plot(hmi_dt,np.array(data[7],dtype=float))

for i in range(data.shape[0]):
    ax=plt.plot(hmi_dt,data[i]/np.max(data[i]),color=sns_cl3[i],label=f'Box {i}')#
    #ax.set_yscale('log')

plt.yscale('log')
plt.legend()
plt.savefig('all_hmi_boxes.png',dpi=300)
plt.close()

for i in range(data.shape[0]):
    fig,ax=plt.subplots(figsize=(16,10))
    ax2=ax.twinx()
    
    lc=np.array(data[i])
    
    ax.plot(hmi_dt,lc,color=sns_cl3[i],label=f'Box {i}')# 
    ax2.errorbar(suit_dt,nb4_bright_ct,yerr=bright_ct_er, fmt='o-',color='k',capsize=2,capthick=1,markersize=1.5)
    
    ax2.set_yscale('log')
    ax.set_yscale('log')
    ax.set_ylabel('Total LOS Magnetic flux (in G)')
    ax2.set_ylabel('SUIT enhancemnt count')
    ax.legend()
    time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
    plt.gca().xaxis.set_major_formatter(time_formatter)
    plt.savefig(f'hmi_boxes_{i}.png',dpi=300)
    plt.close()