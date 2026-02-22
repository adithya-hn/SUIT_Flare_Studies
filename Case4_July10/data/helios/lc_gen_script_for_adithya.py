#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.style as style
import sys
import os
import glob
from datetime import datetime
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
from astropy.io import fits


# Easy and lazy  way to generate a bandwise lightcurve : Using existing spectrogram
cdte_spgm_filename='hel1os_cdte_spectra_cdte1.fits';
ff = fits.open(cdte_spgm_filename);     # load the cdte1 or cdte2 spectrogram file
cdte1_spgm = ff[1].data['COUNTS'];    # verify from the fits info that the appropriate HDU is selected

# the mean  relation for CdTe1 channel-to-energy  conversion is : E = 0.172*Ch + 1.55;
# this of course varies with temperature,flux. But for low count rate application this can be used.

# select energy bands to make bandwise lightcurve
band1 = [30,50]; band2 = [50,80]; band3 = [80,170]; band4 = [170,250]; band5 = [60,170]
all_band = [band1,band2,band3,band4,band5];
#they roughly correspond to 5- 10keV, 10- 15keV, 15-30keV,  30-45 keV, 12 - 30 keV
integration_time = ff[1].data['EXPOSURE'];  # the cadence at which spectrogram is generated
    
cdte1_lc = np.zeros((len(cdte1_spgm),len(all_band)));
for b in range(len(all_band)):
    cdte1_lc[:,b] = np.sum(cdte1_spgm[:,all_band[b][0]:all_band[b][-1]],axis=1)/integration_time


cdte1_countrate_bandlc = cdte1_lc #/integration_time;     # now it is in the units of counts/s

plt.plot()


# The actual way for the lc generation requires event data with UT time-stamp.

CHANNELS = np.arange(512)
gain1 = 0.172; offset1 = 1.55;   # these are to be used for CdTe1 only and for low count rates
arf_cdte = fits.open('CdTeResponseReader/hel1os_cdte_arf_v03.fits')

print(arf_cdte[1].ENERG_LO)
# To convert the spectrogram into flux : first select an energy band 
# this generates bandflux in W/m2
cdte1_band_flux =  np.zeros((len(cdte1_spgm)))
e_low = 10; e_high =20; # change it based on requirement
for k in range(len(cdte1_spgm)):
    cdte1_bin_centres = (CHANNELS*gain1) + offset1;      # gives the energy bin centres of the spectrum
    cdte1_binwidth = np.mean(np.diff(cdte1_bin_centres)) # in keV
    cdte1_bin_centres_inJ = cdte1_bin_centres*1.60218E-16;    # in joules
    arf_ene_vec = np.mean([arf_cdte[1].data['ENERG_LO'],arf_cdte[1].data['ENERG_HI']],axis=0);  # if the ARF is in fits format.
    #int_arf_cdte1 = np.interp(cdte1_bin_centres,arf_cdte[0],arf_cdte[1])*1E-4;    # in units of sq cm, to scale it to a single CdTe(ARF includes area of 2 cdtes)
    cdte1_corr_rate = cdte1_spgm[k,:]/integration_time[k];    # spectrum in units of counts/s
    cdte1_power = (cdte1_corr_rate*cdte1_bin_centres_inJ)/(arf_ene_vec); # into J/s.sqm
    cdte1_band_flux[k] = np.sum(cdte1_power[np.logical_and(cdte1_bin_centres_inJ>=e_low,cdte1_bin_centres_inJ<=e_high)]);   # integral of flux in 10 to 20 keV 
   