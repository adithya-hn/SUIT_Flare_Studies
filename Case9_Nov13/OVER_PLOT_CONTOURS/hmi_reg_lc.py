import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

data=(np.loadtxt('hmi_lc.csv',delimiter=',',skiprows=1,dtype='str')).transpose()
hmi_dt=np.array(data[0],dtype='datetime64')
sns_cl3=sns.color_palette('bright')
fig,ax=plt.subplots(figsize=(14,8))
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
plt.show()

for i in range(data.shape[0]):
    fig,ax=plt.subplots(figsize=(14,8))
    ax=plt.plot(hmi_dt,data[i],color=sns_cl3[i],label=f'Box {i}')#
    #ax.set_yscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig(f'hmi_boxes_{i}.png',dpi=300)
    plt.show()