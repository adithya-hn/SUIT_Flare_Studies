#Data created on 14th Oct 2025
#Author : Adithya HN

#Script to create light curves for different filters in a specified ROI and background box
#The ROI and background box coordinates can be modified as per requirement

#Implementing ROI counter rotation to keep account for submap issues.



import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import pandas as pd
import datetime
import timeit
import pathlib
from astropy.coordinates import SkyCoord
import numpy as np
import glob
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from scipy import stats
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

start = timeit.default_timer()


# ----------Input Parameters------

fdir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/aligned_crop/'
nb3_qs=np.loadtxt('csv_files/NB04QS_data.csv', delimiter=',', skiprows=1, dtype='str').transpose()
nb8_qs=np.loadtxt('csv_files/NB08QS_data.csv', delimiter=',', skiprows=1, dtype='str').transpose()
nb3_thresh=np.array(nb3_qs[3], dtype=float)+3*np.array(nb3_qs[4], dtype=float)
nb8_thresh=np.array(nb8_qs[3], dtype=float)+3*np.array(nb8_qs[4], dtype=float)

fol_nm = os.getcwd() + '/lc_images/'
Filters = ['NB03','NB04','NB08']

for fltr in Filters:

    files = sorted(glob.glob(fdir+ '*3'+fltr+'.fits'))
    if not files:
        print(f"No files found for filter {fltr}")
        continue

    print('Processing filter:', fltr)
    
    dates = []
    box_pth = f'{fol_nm}/{fltr}/Box'
    pathlib.Path(f'{fol_nm}/{fltr}').mkdir(parents=True, exist_ok=True)
    pathlib.Path(box_pth).mkdir(parents=True, exist_ok=True)
    pathlib.Path('csv_files').mkdir(parents=True, exist_ok=True)

    fltr_count = []
    date_array = []
    qs = []
    area=[]
    brightenings=[]
    bright_pixels=[]
    exposure_times=[]

    Sequence = sunpy.map.Map(files, sequence=True)
    aligned_maps = Sequence

    for i, suit_map in enumerate(Sequence):
        fnm = os.path.basename(files[i])[:-5]
        F_name = f'{fol_nm}/{fltr}/{fnm}.jpg'
        exposure = suit_map.meta.get('CMD_EXPT')

        valid = suit_map.data[(suit_map.data > 100)]
        valid_int = valid.astype(int)  
        mode_val = stats.mode(valid_int, keepdims=True).mode[0]*1000/exposure  #normalised mode value    
        #print('Mode value: ', mode_val,'Threshold:', nb3_thresh[i]) 
        threshold = nb3_thresh[i]  # Threshold for brightening detection
        if fltr == 'NB08':
            threshold = nb8_thresh[i] 
            print('NB08 Threshold:', threshold)
        data=suit_map.data* 1000 / exposure
 
        fltr_count.append(np.sum(data))
        qs.append(mode_val)
        date_array.append(suit_map.date.datetime)
        brightenings.append(np.sum(data> threshold))  
        bright_pixels.append(np.sum(data[data > threshold]))
        area.append(suit_map.data.shape[0] * suit_map.data.shape[1])
        exposure_times.append(exposure)
        img=np.where(data> threshold, data, 0)
        plt.imshow(img, origin='lower', cmap='gray', vmin=threshold, vmax=mode_val*2)
        plt.colorbar(label='DN/s')
        plt.title(f'Brightenings in {fltr} at {suit_map.date}')
        plt.savefig(f'{box_pth}/{fnm}_brightenings.jpg', dpi=200)
        plt.close()


    # Save light curve
    np.savetxt(f'csv_files/c6_{fltr}_lc_data.csv',
               np.c_[date_array, fltr_count, qs,brightenings,bright_pixels, area,exposure_times],
               delimiter=',',header='Time,Total,Mode,Bright_count,Bright_area,Area,Exposure',comments='' ,fmt='%s')
    
    # Plot light curve
    plt.rcParams["figure.dpi"] = 120
    fig, ax1 = plt.subplots(1, 1, figsize=(10, 5))
    ax2 = ax1.twinx()
    ax1.plot(date_array, fltr_count, color='C0', label='Total Intensity')
    ax2.plot(date_array, brightenings, color='C1', label='Brightenings Intensity')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Total Intensity (DN)', color='C0')
    ax2.set_ylabel('Brightenings Intensity (DN/s)', color='C1')
    plt.title(f'Light Curve for Filter {fltr}')
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    plt.savefig(f'{fltr}_light_curve.jpg', dpi=200)
    plt.close()

    # Plot light curve
    plt.rcParams["figure.dpi"] = 120
    fig, ax2 = plt.subplots(1, 1, figsize=(10, 5))
    ax2.plot(date_array, qs, color='C1', label='Mode Intensity')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Mode Intensity (DN/s)', color='C1')
    plt.title(f'Light Curve for Filter {fltr}')
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    plt.savefig(f'{fltr}_QS_light_curve.jpg', dpi=200)
    plt.close()
stop = timeit.default_timer()
print('Run Time: ', (stop - start)/60, 'Mins')
