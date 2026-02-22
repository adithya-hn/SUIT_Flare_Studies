
'''
Created on 28 Dec 2025
Author: adithya-hn
Purpose: Derotate and crop the AIA frame aligned  SUIT images


'''


import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from sunpy.map import MapSequence
from sunpy.net import Fido
import glob
import datetime
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
from sunkit_image.coalignment import calculate_match_template_shift,apply_shifts
from datetime import timedelta
import timeit
import pathlib
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
#import ImagesToMovie_pkg
import matplotlib.image as mpimg
from PIL import Image
import pandas as pd
from scipy.ndimage import shift
from scipy.ndimage import gaussian_filter
from sunpy.time import parse_time
from astropy.time import Time
from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord
from tqdm import tqdm
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from sunpy.visualization.animator import MapSequenceAnimator
from sunpy.coordinates import Helioprojective, SphericalScreen, propagate_with_solar_surface
from scipy import stats

import warnings
import logging

warnings.simplefilter('ignore')
log_ = logging.getLogger('sunpy')
log_.setLevel('WARNING')
logging.getLogger('reproject').setLevel(logging.WARNING)

#------------------------------------------------------------------------------#


suit_raw_files= '../data/aligned_fits/'
fltr='NB04'
fltr_fl = glob.glob(suit_raw_files + '*3'+f'{fltr}.fits')
fltr_fl=sorted(fltr_fl, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

pathlib.Path('../data/derotated_pngs').mkdir(exist_ok=True,parents=True)
pathlib.Path('../data/aligned_crop_pngs').mkdir(exist_ok=True,parents=True)
pathlib.Path('../data/aligned_crop_fits').mkdir(exist_ok=True,parents=True)

trX=[]
trY=[]
blX=[]
blY=[]
dates=[]
ref_baseMap=Map(fltr_fl[0])
print('Reference derot frame: ',ref_baseMap.date)
maps=[]
raw_maps=[]
print('No. of SUIT images: ',len(fltr_fl))
print('Derotating images')
for i in tqdm(range(len(fltr_fl))): #
    suit_map=Map(fltr_fl[i])
    #-- Derotate aia maps --
    with propagate_with_solar_surface():
        SuitMap=suit_map.reproject_to(ref_baseMap.wcs) #,parallel=True,dask_method='memmap'
    SuitMap.meta['f_name']=suit_map.meta['f_name']
    SuitMap.meta['dhobt_dt'] = suit_map.meta['dhobt_dt']
    SuitMap.meta['grt_dt']= suit_map.meta['grt_dt']
    SuitMap.meta['mfgdate'] = suit_map.meta['mfgdate']
    SuitMap.meta['SHTR_STR'] = suit_map.meta['SHTR_STR']
    SuitMap.meta['SHTR_STO'] = suit_map.meta['SHTR_STO']
    SuitMap.meta['EXPTIME'] = suit_map.meta['EXPTIME']
    SuitMap.meta['WAVELNTH']= suit_map.meta['WAVELNTH']
    SuitMap.meta['crtime'] = suit_map.meta['crtime']
    SuitMap.meta['cmd_expt'] = suit_map.meta['cmd_expt']
    SuitMap.meta['meas_exp'] = suit_map.meta['meas_exp']
    SuitMap.meta['x2'] = suit_map.meta['x2']
    SuitMap.meta['y2'] = suit_map.meta['y2']
    SuitMap.meta['x1'] = suit_map.meta['x1']
    SuitMap.meta['y1'] = suit_map.meta['y1']
    SuitMap.meta['FLAT_CF'] = suit_map.meta['FLAT_CF']
    SuitMap.meta['SCAT_CF'] = suit_map.meta['SCAT_CF']
    SuitMap.meta['MISSION'] = suit_map.meta['MISSION']
    SuitMap.meta['PNAME'] = suit_map.meta['PNAME']
    SuitMap.meta['DESCR'] = suit_map.meta['DESCR']
    SuitMap.meta['WAVEUNIT'] = suit_map.meta['WAVEUNIT']
    SuitMap.meta['DATE-OBS'] = suit_map.meta['DATE-OBS']
    SuitMap.meta['R_SUN'] = suit_map.meta['R_SUN']
    SuitMap.meta['P_ANGLE'] = suit_map.meta['P_ANGLE']
    SuitMap.meta['DATE-AVG'] = suit_map.meta['DATE-OBS']
    data=SuitMap.data
    nonzero = np.argwhere(data > 100)  # Find nonzero pixels
    ymin, xmin = nonzero.min(axis=0) # Get bounding box
    ymax, xmax = nonzero.max(axis=0)
    fig=plt.figure()
    ax=fig.add_subplot(111,projection=ref_baseMap)
    SuitMap.plot(axes=ax)
    plt.scatter([xmin], [ymin], c='red', label='Top-right')
    plt.scatter([xmax], [ymax], c='blue', label='Bottom-right')
    plt.savefig(f'../data/derotated_pngs/{str(SuitMap.meta["F_NAME"])[:-5]}.jpg')
    plt.close()
    
    
    trX.append(xmax)
    trY.append(ymax)
    blX.append(xmin)
    blY.append(ymin)
    dates.append(suit_map.date.datetime)
    maps.append(SuitMap)


fig,ax1=plt.subplots(1,1, figsize=(10,5))
ax2=ax1.twinx()


ax1.plot(dates, trX,color='C0', label='Top Right X')
ax2.plot(dates, trY,color='C1',  label='Top Right Y')
ax1.plot(dates, blX,color='C2',  label='Bottom Left X')
ax2.plot(dates, blY,color='C3',  label='Bottom Left Y')

plt.xlabel('Time')
plt.ylabel('Pixel Coordinates')
plt.title('Stability of Image Alignment Over Time')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
plt.xticks(rotation=45)
plt.savefig('alignment_stability.png', dpi=200)
plt.close()

print(f'Crop coordinates: x1 {max(blX)}+10, y1 {max(blY)}+10, x2 { min(trX)}-10, y2 { min(trY)}-10')

# x1,y1,x2,y2=100,130,610,750 #FOR CONTINUM CHANNELS
x1,y1,x2,y2=max(blX)+30,max(blY)+30,min(trX)-30,min(trY)-30
#os._exit(0) #------------------

print('Cropping out common frames')
for m in maps:
    cropped = m.submap(bottom_left = [x1,y1]*u.pix,top_right = [x2,y2]*u.pix)
    fig=plt.figure()
    ax=fig.add_subplot(111,projection=cropped)
    cropped.plot()
    plt.savefig(f'../data/aligned_crop_pngs/{str(cropped.meta["F_NAME"])[:-4]}jpg',dpi=200)
    plt.close()
    cropped.save(f'../data/aligned_crop_fits/{str(cropped.meta["F_NAME"])}',overwrite=True)

