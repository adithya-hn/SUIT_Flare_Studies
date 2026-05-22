
import xspec #run this to make sure everything is working correctly
#other packages 
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from datetime import datetime as dt
import matplotlib.pyplot as plt
#surpress warnings for demo...
import warnings
warnings.filterwarnings("ignore")


file='stx_spectrum_2410315184.fits'
stixf= fits.open(file) 
stixf.info()
# print(stixf[1].header)

#mjd2any returns time in seconds from 1-jan-1979
def mjd2any(timezero,spectime):
    '''MJD in days to MJD in seconds, given reference time (timezero) and spectrum time (spectime)'''
    return ((timezero)*86400.+spectime)/86400.


reftime = stixf[1].header['MJDREF']
timezero = stixf[1].header['TIMEZERO'] #this is 1-jan-1979 in MJD, so no need to add MJDREF to it as in COMMENT

rate=stixf[1].data['RATE']
stat_err=stixf[1].data['STAT_ERR']
exp=stixf[1].data['EXPOSURE']
livetime=stixf[1].data['LIVETIME']
print(len(livetime))
spectime=stixf[1].data['TIME'] #COMMENT absTime[i] = mjd2any(MJDREF + TIMEZERO) + TIME[i]
timedel=stixf[1].data['TIMEDEL'] #seconds
rate.shape,spectime.shape

tt=Time(mjd2any(timezero,spectime)+reftime,format='mjd')
timevec=tt.to_value('datetime')

#for plotting - get the energy bins from ENEBAND data (index 2)
emin=stixf[2].data['E_MIN']
emax=stixf[2].data['E_MAX']
ylabels=[f"{n:.0f}-{x:.0f}" for n,x in zip(emin,emax)]

np.empty(rate.T.shape).shape
#matplotlib for nbviewer...
fig,ax=plt.subplots(figsize=[9,5])
cbar=ax.imshow(np.log10(rate.T),origin='lower')
im_ratio=29./77.
#plt.colorbar(cbar, label='Count Rate',fraction=0.046*im_ratio, pad=0.04)
# ax.set_yticks(np.arange(len(ylabels)))
# ax.set_yticklabels(ylabels)
# ax.set_ylabel('Energy Bin (keV)')
# ax.set_xlabel('Index')
# ax.set_title(f"Spectrogram {timevec[0]:%Y-%m-%d %H:%M:%S}")
# plt.close()


# stixr[1].header
# tresp_lo=stixr[1].data['ENERG_LO']
# tresp_hi=stixr[1].data['ENERG_HI']
# trmatrix=stixr[1].data['MATRIX']
# fig,ax=plt.subplots(figsize=[9,5])
# for i,chan in enumerate(ylabels):
#     ax.plot(tresp_lo,trmatrix[:,i],label=f"{chan} keV")
# ax.legend(loc=1, prop={'size': 6})
# ax.set_title('STIX response matrix')
# ax.set_xlabel('Lower Energy Bound (keV)')
# plt.show()
