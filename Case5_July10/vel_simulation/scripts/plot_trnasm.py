
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
from plots_styl import set_pub_style
set_pub_style()

data1 = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/raw/NB8_1_inband.txt")
data2 = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/raw/NB8_1_oob.txt")

col2_file1 = data1[:, 1].astype(float)
col2_file2 = data2[:, 1].astype(float)


# Set negative values to zero
col2_file1[col2_file1 < 0] = 0
col2_file2[col2_file2 < 0] = 0

wavelength=data1[:,0]
# Plot for visualization
fig,axs=plt.subplots(1,1, figsize=(11,5))
axs.plot(wavelength, col2_file1, label="In band", linestyle="dashed",linewidth=0.5)

#ax2=axs.twinx()
#ax2.plot(wavelengths_response,transmission)
#ax2.set_ylabel('$\%$ Transmission')
axs.set_ylabel('Intensity')
axs.set_xlabel("Wavelength")
plt.title('Theoritical spectra and Filter Transmission')
#plt.xlabel("Wavelength")
print(argmax())
plt.legend()
plt.savefig('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/results/in_transm.png')
plt.show()