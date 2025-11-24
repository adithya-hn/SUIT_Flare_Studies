
#other packages 
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from datetime import datetime as dt

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt

file='stx_srm_2410315184.fits'
stixr= fits.open(file) #not sure why with open() as : syntax doesn't work here?
stixr.info()
print(stixr[1].header)
file='stx_spectrum_2410315184.fits'
stixf= fits.open(file) #not sure why with open() as : syntax doesn't work here?
#stixf.info()

cdte_arf_hdul = fits.open('hel1os_cdte_arf_v03.fits')
inp_enes = np.mean(np.array([cdte_arf_hdul[1].data['ENERG_LO'], cdte_arf_hdul[1].data['ENERG_HI']]).T, axis=1)

tresp_lo=stixr[1].data['ENERG_LO']
tresp_hi=stixr[1].data['ENERG_HI']
trmatrix=stixr[1].data['MATRIX']

emin=stixf[2].data['E_MIN']
emax=stixf[2].data['E_MAX']

ylabels=[f"{n:.0f}-{x:.0f}" for n,x in zip(emin,emax)]

fig,ax=plt.subplots(figsize=[9,5])

for i,chan in enumerate(ylabels):
    ax.plot(tresp_lo,trmatrix[:,i],label=f"{chan} keV")
ax.legend(loc=1, prop={'size': 6})
ax.set_title('STIX response matrix')
ax.set_xlabel('Lower Energy Bound (keV)')
#ax.plot(inp_enes,cdte_arf_hdul[1].data['SPECRESP'], linewidth=1)
plt.show()

cdte_srf_hdul=fits.open('hel1os_cdte_srf_v03.fits')

Tresp_lo=cdte_srf_hdul[1].data['ENERG_LO']
Tresp_hi=cdte_srf_hdul[1].data['ENERG_HI']
Trmatrix=cdte_srf_hdul[1].data['MATRIX']

Emin=cdte_srf_hdul[2].data['E_MIN']
Emax=cdte_srf_hdul[2].data['E_MAX']
Ylabels=[f"{n:.0f}-{x:.0f}" for n,x in zip(Emin,Emax)]
#print(Ylabels)
fig,ax=plt.subplots(figsize=[9,5])

for i,chan in enumerate(Ylabels):
    ax.plot(Tresp_lo,Trmatrix[:,i],label=f"{chan} keV")
ax.legend(loc=1, prop={'size': 6})
ax.set_title('HEL1OS response matrix')
ax.set_xlabel('Lower Energy Bound (keV)')
plt.show()



# plt.imshow(np.log10(cdte_srf_hdul[1].data['MATRIX']), extent=[cdte_srf_hdul[2].data['E_MIN'][0],
#  cdte_srf_hdul[2].data['E_MAX'][-1], cdte_srf_hdul[1].data['ENERG_HI'][-1], cdte_srf_hdul[1].data['ENERG_LO'][0]],cmap='coolwarm')

# plt.figure()
# plt.semilogy(inp_enes,cdte_arf_hdul[1].data['SPECRESP'], linewidth=3)
# plt.xlabel('Energy (keV)')
# plt.ylabel('Effective area (cm$^2$)')
# plt.title('CdTe Effective area')
# plt.xlim([0,60])
# plt.grid()
# plt.yscale('log')
# plt.show()