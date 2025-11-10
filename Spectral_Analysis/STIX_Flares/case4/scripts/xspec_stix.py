import xspec
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from datetime import datetime as dt
import os
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
import glob
from stixdcpy.science import Spectrogram, spec_fits_crop, spec_fits_concatenate
from stixdcpy.net import FitsQuery

start_time='2022-03-11T22:00:00'
end_time=  '2022-03-11T23:00:00'

start_dt = Time(start_time).to_datetime()
end_dt = Time(end_time).to_datetime()

#surpress warnings for demo...
import warnings
warnings.filterwarnings("ignore")

project_path= os.path.abspath('..')
#print(f"Project path: {project_path}")
image_list= glob.glob(os.path.join(project_path, 'data/*.fits'))

spec1_file =image_list[0]
spec1 = Spectrogram(spec1_file, None)
spec1.peek()
plt.close()


fl='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/STIX_Flares/case4/data/ql/solo_L1_stix-ql-spectra_20220311_V02.fits'



'''
from plot_stix_spectrogram import *
plot_stix_spec(fl, tickinterval=1000)
fig=plot_stix_spec(fl, mode='scatter')
fig.update_layout(yaxis_title = 'Raw counts', yaxis_range = [0,45000])
cropped_spec = spec_fits_crop(fl, start_time,end_time)
plot_stix_spec(cropped_spec, tickinterval=200)



#stixf= fits.open(file) #not sure why with open() as : syntax doesn't work here?
#stixf.info()

stixf= fits.open(fl) #not sure why with open() as : syntax doesn't work here?
stixf.info()
print(stixf.header)

#mjd2any returns time in seconds from 1-jan-1979
def mjd2any(timezero,spectime):
    MJD in days to MJD in seconds, given reference time (timezero) and spectrum time (spectime)
    return ((timezero)*86400.+spectime)/86400.
#reftime = stixf[1].header['MJDREF']
#timezero = stixf[1].header['TIMEZERO'] #this is 1-jan-1979 in MJD, so no need to add MJDREF to it as in COMMENT
rate=stixf[1].data['RATE']
#stat_err=stixf[1].data['STAT_ERR']
#exp=stixf[1].data['EXPOSURE']
#livetime=stixf[1].data['LIVETIME']
spectime=stixf[1].data['TIME'] #COMMENT absTime[i] = mjd2any(MJDREF + TIMEZERO) + TIME[i]
timedel=stixf[1].data['TIMEDEL'] #seconds
rate.shape,spectime.shape
#tt=Time(mjd2any(timezero,spectime)+reftime,format='mjd')
#timevec=tt.to_value('datetime')'''

xspec.AllData.clear()
xspec.AllData(fl)
#xspec.AllData.nGroups, xspec.AllData.nSpectra #number of data groups and spectra