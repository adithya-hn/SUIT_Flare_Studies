import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

Solexs=(np.loadtxt(f'fit_results_AL1_SOLEXS_20240602_SDD2_L1_2406020630_2406020915_TEMP_EM.txt',skiprows=1,dtype='str')).transpose()
goes=np.genfromtxt(f'goes_temp_em.csv', delimiter=',',names=True, dtype=None, encoding='utf-8')
goes_2=np.genfromtxt(f'tn_woods_goes_temp_em.csv', delimiter=',',names=True, dtype=None, encoding='utf-8')

#---------------------------------------------------------%%%%%%%SOLEXS%%%%%%------------------------------------------#
time_array4=[]
sl_temp=[float(tp) for tp in Solexs[1]]
sl_temp_er=[float(tpe) for tpe in Solexs[2]]
sl_Em=[float(em)*1e46 for em in Solexs[3]]
sl_Em_er=[float(eme) for eme in Solexs[4]]
slt=Solexs[0]
time_array4=[datetime.strptime(str(ts)[:19], "%Y-%m-%dT%H:%M:%S") for ts in slt]
#---------------------

goes_dt=np.array(goes['time'],dtype='datetime64')
goes_temp=np.array(goes['temperature'])

goes2_dt=np.array(goes_2['date_time'],dtype='datetime64')
goes2_temp=np.array(goes_2['temperature'])

# float_list = [float(num_str) for num_str in goes_temp]
# print(float_list)
# print(goes_temp)
fig,ax=plt.subplots(figsize=(14,8))
ax2=ax.twinx()
ax.plot(time_array4,sl_temp)
ax.plot(goes_dt,goes_temp)
ax.plot(goes2_dt,goes2_temp)
# ax2.plot(time_array4,sl_Em)
# ax2.plot(goes_dt,goes['emission_measure'])
#ax2.plot(goes2_dt,goes_2['em'])
plt.show()
