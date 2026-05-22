#Data created on 14th Oct 2025
#Author : Adithya HN

#Script to create light curves for different filters in a specified ROI and background box
#The ROI and background box coordinates can be modified as per requirement




import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import pandas as pd
import datetime
import timeit
import pathlib
import numpy as np
import glob
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from scipy import stats
from sys import path as sys_path
from astropy.convolution import convolve

sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
#set_pub_style()

start = timeit.default_timer()


# ----------Input Parameters------

fdir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/aligned_crop_fits/'
fol_nm = os.getcwd() + '/lc_images/'
Filters = ['NB03','NB04','NB08']

#----------------------------------

def gaussian_kernel(size, sigma):
    ax = np.linspace(-(size // 2), size // 2, size)
    xx, yy = np.meshgrid(ax, ax)
    kernel = np.exp(-(xx**2 + yy**2) / (2 * sigma**2))
    return kernel / np.sum(kernel)

for fltr in Filters:

    files = sorted(glob.glob(fdir+ '*3'+fltr+'.fits'))
    if not files:
        print(f"No files found for filter {fltr}")
        continue

    print('Processing filter:', fltr)
    
    pathlib.Path(f'{fol_nm}/{fltr}').mkdir(parents=True, exist_ok=True)
    pathlib.Path('csv_files').mkdir(parents=True, exist_ok=True)
    dates = []
    fltr_count = []
    date_array = []
    meas_exp = []
    area=[]
    exposure_times=[]
    meas_exp_times=[]
    smoothed_count=[]

    Sequence = sunpy.map.Map(files, sequence=True)
    aligned_maps = Sequence
    kernel = gaussian_kernel(3, sigma=2) 

    for i, suit_map in enumerate(Sequence):
        exposure = suit_map.meta.get('CMD_EXPT')
        meas_exp = suit_map.meta.get('MEAS_EXP')
        data     = suit_map.data* 1000 / exposure
        sm_data  = convolve(suit_map.data, kernel, boundary='extend')* 1000 / exposure #smoothed data

        fltr_count.append(np.sum(data))
        smoothed_count.append(np.sum(sm_data))
        date_array.append(suit_map.date.datetime)
        area.append(suit_map.data.shape[0] * suit_map.data.shape[1])
        exposure_times.append(exposure)
        meas_exp_times.append(meas_exp)

    # Save light curve
    np.savetxt(f'csv_files/c4_{fltr}_lc_data.csv',
               np.c_[date_array, fltr_count,smoothed_count,area,exposure_times,meas_exp_times],
               delimiter=',',header='Time,Total,Area,Exposure,Meas_exposure',comments='' ,fmt='%s')
    


stop = timeit.default_timer()
print('Run Time: ', (stop - start)/60, 'Mins')

