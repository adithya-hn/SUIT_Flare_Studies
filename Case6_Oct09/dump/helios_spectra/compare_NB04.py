


'''
Created on 4 Jan 2026
@author: adithya-hn
purpose: compare diff count img count

'''
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
from scipy.ndimage import median_filter

suit_fl="Diff_img_data_NB04.csv"
data = np.genfromtxt(suit_fl, delimiter=',', names=True, dtype=None)

t = np.array(data['Date'], dtype='datetime64[ms]')
area=np.array(data['area'],dtype=float)
diff_c = np.array(data['diff_count'], dtype=float)#/area
img_c  = np.array(data['Img_count'], dtype=float)#/area

fig,ax=plt.subplots(figsize=(14,6))
ax2=plt.twinx()

ax.plot(t,diff_c,'bo-',label='Difference image count',markersize=4)
ax2.plot(t,img_c,'ro-',label='Image contour count',markersize=4)

plt.title(f'Comapre diff count vs image count')
ax.set_yscale('log')
ax2.set_yscale('log')
plt.savefig(f'compare_NB04.png',dpi=300)
plt.show()
