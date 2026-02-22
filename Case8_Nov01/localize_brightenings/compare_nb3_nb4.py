# read diffrence image intensity from csv file

import numpy as np
import matplotlib.pyplot as plt




nb3_diff_img=np.genfromtxt('Diff_img_data_NB03.csv',delimiter=',',names=True,dtype=None, encoding='utf-8')
nb4_diff_img=np.genfromtxt('Diff_img_data_NB04.csv',delimiter=',',names=True,dtype=None, encoding='utf-8')

nb3_date=np.array(nb3_diff_img['Date'],dtype='datetime64')
nb4_date=np.array(nb4_diff_img['Date'],dtype='datetime64')

nb3_count=np.array(nb3_diff_img['diff_count'],dtype=float)
nb4_count=np.array(nb4_diff_img['diff_count'],dtype=float)


plt.plot(nb3_date,nb3_count,'bo-',markersize=1.5)
plt.plot(nb4_date,nb4_count,'ro-',markersize=1.5)
plt.show()
