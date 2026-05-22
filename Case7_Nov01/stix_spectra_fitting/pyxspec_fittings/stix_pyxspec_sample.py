
"""
@Author      : Adithya H N
@Created On  : 2026-01-25
@Last Updated: 2026-02-21
@Project     : Pre-flare study
@Version     : 1.0

@Description
-----------
Brief description: Spectral fitting for pre-flare transients,  trying to fit thermal(apec/vapec) and non-thermal (thick target) model.
The current file isn use is type II file whic contain multiple pha files in one file
SRM file accounts for both RMF and ARF, 
    
"""


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


# file='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/stix_spectra_fitting/pyxspec_fittings/stix_012230_012300.pha'
# stixf= fits.open(file) 
# stixf.info()
# print(stixf[1].header)

# #mjd2any returns time in seconds from 1-jan-1979
# def mjd2any(timezero,spectime):
#     '''MJD in days to MJD in seconds, given reference time (timezero) and spectrum time (spectime)'''
#     return ((timezero)*86400.+spectime)/86400.


# reftime = stixf[1].header['MJDREF']
# timezero = stixf[1].header['TIMEZERO'] #this is 1-jan-1979 in MJD, so no need to add MJDREF to it as in COMMENT

# rate=stixf[1].data['RATE']
# #stat_err=stixf[1].data['STAT_ERR']
# #exp=stixf[1].data['EXPOSURE']
# #livetime=stixf[1].data['LIVETIME']
# spectime=stixf[1].data['TIME'] #COMMENT absTime[i] = mjd2any(MJDREF + TIMEZERO) + TIME[i]
# timedel=stixf[1].data['TIMEDEL'] #seconds
# rate.shape,spectime.shape

# tt=Time(mjd2any(timezero,spectime)+reftime,format='mjd')
# timevec=tt.to_value('datetime')

# #for plotting - get the energy bins from ENEBAND data (index 2)
# emin=stixf[2].data['E_MIN']
# emax=stixf[2].data['E_MAX']
# ylabels=[f"{n:.0f}-{x:.0f}" for n,x in zip(emin,emax)]

# np.empty(rate.T.shape).shape
# #matplotlib for nbviewer...
# fig,ax=plt.subplots(figsize=[9,5])
# cbar=ax.imshow(np.log10(rate.T),origin='lower')
# im_ratio=29./77.
# #plt.colorbar(cbar, label='Count Rate',fraction=0.046*im_ratio, pad=0.04)
# # ax.set_yticks(np.arange(len(ylabels)))
# # ax.set_yticklabels(ylabels)
# # ax.set_ylabel('Energy Bin (keV)')
# # ax.set_xlabel('Index')
# # ax.set_title(f"Spectrogram {timevec[0]:%Y-%m-%d %H:%M:%S}")
# # plt.close()


# # stixr[1].header
# # tresp_lo=stixr[1].data['ENERG_LO']
# # tresp_hi=stixr[1].data['ENERG_HI']
# # trmatrix=stixr[1].data['MATRIX']
# # fig,ax=plt.subplots(figsize=[9,5])
# # for i,chan in enumerate(ylabels):
# #     ax.plot(tresp_lo,trmatrix[:,i],label=f"{chan} keV")
# # ax.legend(loc=1, prop={'size': 6})
# # ax.set_title('STIX response matrix')
# # ax.set_xlabel('Lower Energy Bound (keV)')
# # plt.show()

print('1')
xspec.AllData.clear()
xspec.AllData("1:1 stix_012230_012300.pha") #the row index (from 1) in curly braces tells Xspec which row to load
xspec.AllData.nGroups, xspec.AllData.nSpectra #number of data groups and spectra

xspec.Plot.setGroup("1") #("1 2 3") if using 3 rows, for example
xspec.Plot.add=False 
xspec.Plot.device = '/null'
xspec.Plot.xAxis = "keV"
xspec.Plot('data') 
x1=xspec.Plot.x()
y1=xspec.Plot.y()
y1err=xspec.Plot.yErr()

fig,ax=plt.subplots()
# ax.errorbar(x1,y=y1,yerr=y1err,label=f'{timevec[12]:%Y-%m-%d %H:%M:%S}',marker='+',linestyle='none')
# ax.set_title(f'STIX spectrum at {timevec[12]:%Y-%m-%d %H:%M:%S}')
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Counts')
ax.set_ylim([.01,20000])
ax.set_yscale('log')
ax.set_xlim([0,100])
plt.show()

m1=xspec.Model("bknpower")#apec
m1.setPars({2:"15.0 -.1,,,18"})
p2=getattr(m1.bknpower,'BreakE') 
p2.frozen=False #unfreeze break energy

#select appropriate fit range
fitstart=8.0 #keV
fitend=30.0 #keV
xspec.AllData.ignore(f"0.-{fitstart} {fitend}-**")
xspec.Fit.statMethod = "chi" #Valid names: 'chi' | 'cstat' | 'lstat' | 'pgstat' | 'pstat' | 'whittle'. 
xspec.Fit.query = "no"
xspec.Fit.renorm()
xspec.Fit.nIterations=1000
xspec.Fit.perform()
#print statistic if you can't find it in the notebook or terminal
xspec.Fit.statistic #chi-squared statistic unless using cstat
xspec.Fit.error("maximum 100 1.0 1-3") #in case of Xspec error about chisq being too high, reset maximum to some value
c1=m1.bknpower
idx1=c1.PhoIndx1.values[0]
idx1_lbound=c1.PhoIndx1.error[0]
idx1_ubound=c1.PhoIndx1.error[1]
breakE=c1.BreakE.values[0]
break_lbound=c1.BreakE.error[0]
break_ubound=c1.BreakE.error[1]
idx2=c1.PhoIndx2.values[0]
idx2_lbound=c1.PhoIndx2.error[0]
idx2_ubound=c1.PhoIndx2.error[1]
c1.BreakE.error,c1.PhoIndx1.error 
xspec.AllData.notice("5.0-50.0")
xspec.Plot.setGroup("1") 
xspec.Plot.add=False
xspec.Plot('data') 
title=f"{timevec[12]:%Y-%m-%d %H:%M:%S}"
fittext=f"power-law index 1: {idx1:.2f} ({idx1_lbound:.2f}-{idx1_ubound:.2f})<br>\
break energy: {breakE:.2f} ({break_lbound:.2f}-{break_ubound:.2f}) keV<br>\
power-law index 2: {idx2:.2f} ({idx2_lbound:.2f}-{idx2_ubound:.2f})"

fig,axs=plt.subplots(2,1,figsize=(5.5,7),gridspec_kw=dict( height_ratios=[4,1],hspace=0.05))

axs[0].errorbar(x1,y1,yerr=y1err,ms=5,color='k',fmt='o')
axs[0].plot(xspec.Plot.x(),xspec.Plot.model(),color='firebrick',drawstyle='steps-mid')
axs[0].set_ylabel('STIX counts')
ylim=[0.1,2*10**4]
xlim=[5,50]
axs[0].set_ylim(ylim)
axs[0].set_yscale('log')
for aa in axs:
    aa.set_xlim(xlim)
    aa.label_outer()
    
axs[0].set_title(title)
    
fiter=[fitstart,fitend]
axs[0].plot([fiter[0],fiter[0]],[ylim[0],10**(np.log10(ylim[1]))],':',color='grey')
axs[0].plot([fiter[1],fiter[1]],[ylim[0],10**(np.log10(ylim[1]))],':',color='grey')

# Calculate and plot the residuals on the bottom plot
resid=(np.subtract(xspec.Plot.y(),xspec.Plot.model()))
axs[1].plot(xspec.Plot.x(),resid,'.',ms=5,color='k')
axs[1].set_ylim([-10,10])
axs[1].set_xlabel('Energy [keV]')
axs[1].plot(fiter,[0,0],'--',color='grey')
axs[1].set_ylabel('Residuals')
axs[0].text(20,10000,fittext[:fittext.find('<br>')])
axs[0].text(20,7000,fittext[fittext.find('<br>')+4:fittext.rfind('<br>')])
axs[0].text(20,5000,fittext[fittext.rfind('<br>')+4:])
plt.show()
