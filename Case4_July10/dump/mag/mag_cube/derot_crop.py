
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
from astropy.io import fits

import warnings
import logging

warnings.simplefilter('ignore')
log_ = logging.getLogger('sunpy')
log_.setLevel('WARNING')
logging.getLogger('reproject').setLevel(logging.WARNING)

#------------------------------------------------------------------------------#


suit_raw_files= '/media/adithya/Adi_disk4/SUIT_flare_work/case4_jul10/data/HMI/HMI_cutouts/'
fltr_fl = sorted(glob.glob(suit_raw_files + '*.fits')) 
# pathlib.Path('derotated_pngs').mkdir(exist_ok=True,parents=True)
pathlib.Path('aligned_crop_pngs').mkdir(exist_ok=True,parents=True)
pathlib.Path('aligned_crop_fits').mkdir(exist_ok=True,parents=True)

trX=[]
trY=[]
blX=[]
blY=[]
dates=[]
ref_baseMap=Map(fltr_fl[0])
print(str(os.path.basename(fltr_fl[0]))[:-4])
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
   
    data=SuitMap.data
    nonzero = np.argwhere(data > 100)  # Find nonzero pixels
    ymin, xmin = nonzero.min(axis=0) # Get bounding box
    ymax, xmax = nonzero.max(axis=0)
    # fig=plt.figure()
    # ax=fig.add_subplot(111,projection=ref_baseMap)
    # SuitMap.plot(axes=ax)
    # plt.scatter([xmin], [ymin], c='red', label='Top-right')
    # plt.scatter([xmax], [ymax], c='blue', label='Bottom-right')
    # plt.savefig(f'../data/derotated_pngs/{str(SuitMap.meta["F_NAME"])[:-5]}.jpg')
    # plt.close()
    
    
    trX.append(xmax)
    trY.append(ymax)
    blX.append(xmin)
    blY.append(ymin)
    dates.append(suit_map.date.datetime)
    maps.append(SuitMap)


# fig,ax1=plt.subplots(1,1, figsize=(10,5))
# ax2=ax1.twinx()

# ax1.plot(dates, trX,color='C0', label='Top Right X')
# ax2.plot(dates, trY,color='C1',  label='Top Right Y')
# ax1.plot(dates, blX,color='C2',  label='Bottom Left X')
# ax2.plot(dates, blY,color='C3',  label='Bottom Left Y')

# plt.xlabel('Time')
# plt.ylabel('Pixel Coordinates')
# plt.title('Stability of Image Alignment Over Time')
# plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
# plt.xticks(rotation=45)
# plt.savefig('alignment_stability.png', dpi=200)
# plt.close()

print(f'Crop coordinates: x1 {max(blX)}+10, y1 {max(blY)}+10, x2 { min(trX)}-10, y2 { min(trY)}-10')

# x1,y1,x2,y2=100,130,610,750 #FOR CONTINUM CHANNELS
x1,y1,x2,y2=max(blX)+30,max(blY)+30,min(trX)-30,min(trY)-30
#os._exit(0) #------------------
cube=[]
print('Cropping out common frames')
k=0
for m in maps:
    cropped = m.submap(bottom_left = [x1,y1]*u.pix,top_right = [x2,y2]*u.pix)
    fig=plt.figure()
    ax=fig.add_subplot(111,projection=cropped)
    cropped.plot()
    plt.savefig(f'aligned_crop_pngs/{str(os.path.basename(fltr_fl[k]))[:-4]}jpg',dpi=200)
    plt.close()
    cropped.save(f'aligned_crop_fits/{str(os.path.basename(fltr_fl[k]))}',overwrite=True)
    cube.append(cropped.data)
    k+=1

cube = np.array(cube).astype("float32")
hdu = fits.PrimaryHDU(cube)
hdu.writeto("HMI_derot_cube.fits", overwrite=True)