

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import *
from astropy.io import fits
import scipy.misc
import math as mt
from datetime import datetime
import os
import timeit
from scipy import stats as S
import scipy as sp
import pathlib
import pandas as pd
import matplotlib.dates as mdates
from sunpy.time import parse_time

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()



iris_data=(np.loadtxt('iris_sji_2796_lc.csv',delimiter=',',skiprows=1,dtype='str')).transpose()
nb3_data=(np.loadtxt('nb03_lc.csv',delimiter=',',dtype='str',skiprows=1)).transpose()
iris_preflare=(np.loadtxt('pre_flare_iris_sji_2796_lc.csv',delimiter=',',skiprows=1,dtype='str')).transpose()


nb3_date_array=parse_time(nb3_data[0]).datetime
date_array=parse_time(iris_data[0]).datetime
pf_iris_date=parse_time(iris_preflare[0]).datetime

nb3_count=nb3_data[1]
iris_count=iris_data[1]
iris_qs1=iris_data[2]
iris_qs2=iris_data[3]

pf_iris_count=iris_preflare[1].astype('float')
pf_iris_qs1=iris_preflare[2].astype('float')
pf_iris_qs2=iris_preflare[3].astype('float')

print('iris count:',len(iris_qs2),len(iris_qs1))

nb3_count=np.array(nb3_count, dtype=float)
iris_count=iris_count.astype('float')
iris_qs1=iris_qs1.astype('float')
iris_qs2=iris_qs2.astype('float')
#print('IRIS count:',nb3_count,iris_count)

fig, ax = plt.subplots(figsize=(12, 6))
ax2=ax.twinx()

ax.scatter(nb3_date_array,nb3_count, label='NB3', color='blue', marker='o', s=1,linewidth=0.5)
ax.plot(date_array,iris_qs1, label='IRIS_qs1', color='blue', marker='o', markersize=1)

ax.plot(date_array,iris_qs2, label='IRIS_qs2', color='green', marker='o', markersize=1)
ax2.scatter(date_array,iris_count, label='IRIS', color='red', marker='o', s=1,linewidth=0.5)
ax2.scatter(pf_iris_date,pf_iris_count, label='IRIS-preflare', color='red', marker='o', s=1,linewidth=0.5)
#ax.set_ylabel('NB3 Counts')
ax.set_ylabel('IRIS Counts')
ax.set_xlabel('Date')
plt.title('IRIS data')

plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig('../results/iris_nb3_flare_lc.png', dpi=300)
plt.show()
plt.close()