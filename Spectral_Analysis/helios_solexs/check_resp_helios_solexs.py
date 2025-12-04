
#other packages 
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from datetime import datetime as dt

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt



sdd_arf_data = fits.open('solexs_arf_SDD2_v1.arf')
cdte_arf_hdul = fits.open('hel1os_cdte_arf_v03.fits')
inp_enes = np.mean(np.array([cdte_arf_hdul[1].data['ENERG_LO'], cdte_arf_hdul[1].data['ENERG_HI']]).T, axis=1)
solex_inp_enes = np.mean(np.array([sdd_arf_data[1].data['ENERG_LO'], sdd_arf_data[1].data['ENERG_HI']]).T, axis=1)

fig,ax=plt.subplots(figsize=[9,5])
ax.plot(solex_inp_enes,sdd_arf_data[1].data['SPECRESP'],label='SoLEXS')
ax.plot(inp_enes,cdte_arf_hdul[1].data['SPECRESP'],label='HEL1OS' ,linewidth=1)
ax.set_ylabel(r'Effective Area (cm$^2$)')
ax.set_xlabel('Energy (keV)')
ax.legend()
plt.show()