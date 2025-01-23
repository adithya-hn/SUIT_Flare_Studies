

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
from subprocess import call
from matplotlib import colors
import matplotlib.dates as mdates



Helios=(np.load("cdte_data_flare_3.npy", allow_pickle=True)).transpose()

cdte1=Helios[1]
cdte2=Helios[2]
cdte1_er=Helios[3]
datetime_objects = pd.to_datetime(Helios[0])
helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]
helio_sec_array=[(datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f").timestamp()) for ts in datetime_objects]

time_array=[]

rc('axes', linewidth=1.2)
plt.rcParams["xtick.major.size"] = 10
fig,ax3=plt.subplots(1,1, figsize=(11,5))


ax3.spines.right.set_position(("axes", 1.1))
ax3.errorbar(helio_time_array,cdte1,yerr=cdte1_er, fmt='ro',capsize=2,markersize=2,label="Helios-CdTe1")
#ax3.plot(helio_time_array,cdte2, label="Helios")
ax3.set_ylabel('Helios',fontsize=13)
ax3.set_yscale('log')
plt.show()
np.savetxt('helios_timeseries.csv',np.c_[helio_time_array,helio_sec_array,cdte1,cdte1_er],delimiter=',',fmt='%s')