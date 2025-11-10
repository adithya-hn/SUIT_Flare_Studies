'''
Created on 14 oct 2025
Author: Adithya-hn
purpose: create images of fits files

'''



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

fdir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/aligned_crop/'
fol_nm = os.getcwd() + '/lc_images/'
Filters = ['NB08','NB03','NB04']

for fltr in Filters:

    files = sorted(glob.glob(fdir+ '*3'+fltr+'.fits'))
    if not files:
        print(f"No files found for filter {fltr}")
        continue

    print('Processing filter:', fltr)
    
    dates = []
    box_pth = f'{fol_nm}/{fltr}/Box'
    pathlib.Path(f'{fol_nm}/{fltr}').mkdir(parents=True, exist_ok=True)

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
        Imax,Imin=1000,35000
        if suit_map.meta['FTR_NAME']== 'NB08' :
            Imin,Imax=2000,25000
        suit_map.plot(vmin=Imin,vmax=Imax)
        plt.colorbar(label='DN',fraction=0.046, pad=0.04)
        plt.savefig(f'{fol_nm}/{fltr}/{fnm}.png', dpi=300)
        plt.close()