# read diffrence image intensity from csv file

import numpy as np
import matplotlib.pyplot as plt
import glob

files = sorted(glob.glob("csv_files/*NB04.csv"))
o_files = sorted(glob.glob("csv_files/prev_csv_fils/*NB04.csv"))
norm =np.loadtxt('norm_vals.csv',delimiter=',')[1]




nb3_diff_img=np.genfromtxt(files[2],delimiter=',',names=True,dtype=None, encoding='utf-8')
nb4_diff_img=np.genfromtxt(o_files[2],delimiter=',',names=True,dtype=None, encoding='utf-8')

nb3_date=np.array(nb3_diff_img['Date'],dtype='datetime64')
nb4_date=np.array(nb4_diff_img['Date'],dtype='datetime64')

nb3_count=np.array(nb3_diff_img['diff_count'],dtype=float)
nb4_count=np.array(nb4_diff_img['diff_count'],dtype=float)

plt.figure(figsize=(14,6))
plt.plot(nb3_date,nb3_count,'bo-',markersize=1.5,label='new files')
plt.plot(nb4_date,nb4_count,'ro-',markersize=1.5,label='old files')
plt.legend()
plt.show()
