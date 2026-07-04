

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
from astropy.coordinates import SkyCoord
import numpy as np
import glob
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from scipy import stats
from skimage.morphology import disk, closing
from skimage import filters, measure
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
#set_pub_style()

fdir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/raw/'
fol_nm = os.getcwd() + '/lc_images/'
#Filters = ['NB01','NB02','NB03','NB04','NB05','NB06','NB07','NB08']
Filters = ['NB04']

for fltr in Filters:

    files = sorted(glob.glob(fdir+ '*3'+fltr+'.fits'))
    if not files:
        print(f"No files found for filter {fltr}")
        continue

    print('Processing filter:', fltr)
    Sequence = sunpy.map.Map(files, sequence=True)
    aligned_maps = Sequence
    date_array = []
    shut_pos = []

    for i, suit_map in enumerate(Sequence):
        
        shut_pos.append(suit_map.meta.get('SHTR_STR'))
        date_array.append(suit_map.date.datetime)
    plt.plot(date_array,shut_pos)
    plt.savefig('shutter_position_c2.png',dpi=300)
    plt.show()

