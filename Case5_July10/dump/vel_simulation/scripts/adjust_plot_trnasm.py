
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
from scipy.interpolate import interp1d

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

data1 = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/raw/NB8_1_inband.txt")

col2_file1 = data1[:, 1].astype(float)
wavelength=data1[:,0]
interp_func = interp1d(wavelength, col2_file1, kind='linear' )
wavln=np.arange(3950000,3985000,1)/10000
print(wavln)
intre_inten=interp_func(wavln)
shifted_inten=intre_inten[926:]
print(len(intre_inten)-len(shifted_inten))
pad_zeros=np.zeros(926)
intensity=np.concatenate((shifted_inten,pad_zeros))
print(len(intre_inten),len(intensity))
# Plot for visualization
fig,axs=plt.subplots(1,1, figsize=(11,5))
axs.plot(wavln, intensity, label="In band", linestyle="dashed",linewidth=0.5)

axs.set_ylabel('Intensity')
axs.set_xlabel("Wavelength")
plt.title('Theoritical spectra and Filter Transmission')
plt.show()
print(wavln[argmax(intensity)])#,np.where(wavln==396.8469),max(intre_inten),np.where(intre_inten==0.15220557367168844))

np.savetxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/processed/Shifted_transm_profile.csv',np.c_[wavln,intensity],delimiter=',')