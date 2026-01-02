from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from sunpy.time import parse_time
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
#set_pub_style()



data=np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/nb08_degrade.csv',delimiter=',')

dt_dy=data[:,0]
# dt =parse_time(Time(dt_dy, format='decimalyear').iso)
dt_ =Time(dt_dy, format='decimalyear').to_datetime()

idx = np.argsort(dt_)
dt = dt_[idx]

val=np.array(data[:,1][idx],dtype=float)
plt.figure(figsize=(14,8))
plt.title('NB04 degradation')
plt.scatter(dt,val,marker='+',color='r')
plt.xlabel('Time (UT)')
plt.ylabel('Exposure normalized relative intensity')
#time_formatter = mdates.DateFormatter('%Y')  # Format as HH:MM
#plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig('degradation.png',dpi=300)
plt.show()
