import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

data=(np.loadtxt('QS_data.csv',delimiter=',',dtype='str',skiprows=1)).transpose()

qs=np.array(data[1],dtype=float)
ar=np.array(data[2],dtype=float)
date=np.array(data[0], dtype='datetime64')

fig,axs=plt.subplots(1,1, figsize=(10,5))
axs1=axs.twinx()
axs.plot(date,qs,'co-',markersize=2,linewidth=0.5, label='QS mode')
axs1.plot(date,ar,'bo-',markersize=2,linewidth=0.5, label='AR Sum')
img_nm='qs_ar_lightcurve_with_hmi.png'
axs.set_ylabel('QS Mean counts',fontsize=13)
axs1.set_ylabel('AR total counts',fontsize=13)
plt.xlabel(f"Time ",fontsize=13)
plt.title(f'QS and AR Light curves ({date[0]})')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(img_nm,dpi=300)
plt.close()
