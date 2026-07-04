

"""
@Author      : Adithya H N
@Created On  : 2026-04-05
@Last Updated: 2026-04-05
@Project     : Project Name
@Version     : 1.0

@Description
-----------
Brief description
"""

import numpy as np
import matplotlib.pyplot as plt
import  pandas as pd



thermal_data=np.genfromtxt('7_Nov01_stix_timeResolved_nth_fit.csv',delimiter=',', names=True, dtype=None, encoding='utf-8')
print(thermal_data.dtype.names)
dt1=np.array(thermal_data['time_start'],dtype='datetime64')
T1=thermal_data['T']
T1_er=thermal_data['T_er1']
EM1=thermal_data['EM']
EM1_er=thermal_data['EM_er1']
L1=thermal_data['L']

non_thermal=np.genfromtxt('thermal/7_Nov01_stix_timeResolved_th_fit.csv',delimiter=',', names=True, dtype=None, encoding='utf-8')
dt2=np.array(non_thermal['time_start'],dtype='datetime64')
T2=non_thermal['T']
T2_er=non_thermal['T_er1']
EM2=non_thermal['EM']
EM2_er=non_thermal['EM_er1']
L2=non_thermal['L']


fig,ax=plt.subplots(3, 1, sharex=True, figsize=(12,16), gridspec_kw={'hspace': 0})
ax[0].errorbar(dt1,T1,yerr=T1_er,label='1 model temeperature')
ax[0].errorbar(dt2,T2,yerr=T2_er,label='2 model temeperature')
ax[1].errorbar(dt1,EM1,yerr=EM1_er,label='1 model EM')
ax[1].set_yscale('log')
# ax[0].set_yscale('log')
ax[0].set_ylim(4,30)
ax[0].axhline(10)
ax[1].errorbar(dt2,EM2,yerr=EM2_er,label='2 model EM')
# ax[2].errorbar(dt,,yerr=EM2_er,label='2 model EM')
ax[2].plot(dt1,L1,label='1 model Likelihood')
ax[2].plot(dt2,L2,label='2 model Liklihood')
ax[2].set_ylim(-100,5)

plt.legend()
plt.savefig('Hot_onset_condition.png',dpi=300)
plt.show()



