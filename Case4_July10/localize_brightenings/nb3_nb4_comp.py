import numpy as np
import matplotlib.pyplot as plt


nb3_data=np.loadtxt('180_Diff_img_data_NB03.csv',skiprows=1,delimiter=',',dtype='str').transpose()
nb4_data=np.loadtxt('Diff_img_data_NB03.csv',skiprows=1,delimiter=',',dtype='str').transpose()

print(nb3_data.shape)

nb3_dt=np.array(nb3_data[0],dtype='datetime64[s]')
nb4_dt=np.array(nb4_data[0],dtype='datetime64[s]')


plt.figure(figsize=(10, 6))
plt.plot(nb3_dt, np.array(nb3_data[3],dtype=float),marker='o', label='NB3', color='blue')
plt.plot(nb4_dt, np.array(nb4_data[3],dtype=float), marker='s', label='NB4', color='orange')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('NB3 vs NB4 Data Comparison') 
plt.legend()
plt.grid()
plt.savefig('nb3_vs_nb4_comparison.png')
plt.show()
