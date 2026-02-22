import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.dates as mdates
import pandas as pd
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
import numpy.ma as ma
import seaborn as sns
import glob
import os
from scipy.optimize import curve_fit


scol =sns.color_palette("colorblind")

files = sorted(glob.glob("csv_files/*.csv"))
# norm =np.loadtxt('norm_vals.csv',delimiter=',')[1]
suit_count=[]
helios_count=[]
for f in files:
    data = np.loadtxt(f, delimiter=',')
    print(data[:,0])
    suit_count.append(data[:,0])
    helios_count.append(data[:,1])
print(helios_count)
suit_count=np.concatenate(suit_count)
helios_count=np.concatenate(helios_count)
plt.scatter(suit_count/1e6,helios_count)
plt.show()
    