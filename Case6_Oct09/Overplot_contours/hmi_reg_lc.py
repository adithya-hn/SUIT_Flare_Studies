import numpy as np
import matplotlib.pyplot as plt


data=(np.loadtxt('hmi_lc.csv',delimiter=',',skiprows=1,dtype='str')).transpose()
hmi_dt=np.array(data[0],dtype='datetime64')

fig,ax=plt.subplots(figsize=(14,8))
print(data.shape)
# ax.plot(hmi_dt,np.array(data[3],dtype=float))
for i in range(6):
    ax=plt.plot(hmi_dt,np.array(data[i+1],dtype=float),label=f'Box {i}')#
    #ax.set_yscale('log')
plt.yscale('log')
plt.legend()
plt.show()