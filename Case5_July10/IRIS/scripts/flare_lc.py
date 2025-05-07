#----------------
# Author: Adithya H.N.
# Date: 2025-04-18
#----------------
# Purpose: Light curve from IRIS Mg II k.
# Description: 


import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from datetime import timedelta
import timeit
import pathlib
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
#from irispy.io import read_files  
from sunpy.time import TimeRange


# Set the path to the directory where the files will be saved
path = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/IRIS/Flare_time'
#
# Create the directory if it doesn't exist
if not os.path.exists(path):
    os.makedirs(path)

file='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/IRIS/data/raw/flare_data/iris_l2_20240710_152851_3620108477_SJI_2796_t000.fits'

# Load the IRIS SJI 2796 map
iris_img = fits.open(file)
#iris_img.info()

sji_imgs=iris_img[0].data

#maps=[sunpy.map.Map(sji_imgs[i].data, iris_img[0].header) for i in range(10)]
#print(iris_img[0].header)
ref_map=sunpy.map.Map(sji_imgs[10], iris_img[0].header)
'''
maps=[]
for i in range(len(sji_imgs)):
    maps.append(sunpy.map.Map(sji_imgs[i], iris_img[0].header))

mgIIk_sq=sunpy.map.Map(maps,sequence=True)
ani = mgIIk_sq.peek()
ani.save(path+'/iris_sji_2796.gif', fps=10, writer='imagemagick')'''

fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(projection=ref_map)
#er_coords = SkyCoord(Tx=(-10, 90) * u.arcsec, Ty=(0, -100) * u.arcsec, frame=ref_map.coordinate_frame)
qs_coords1 = SkyCoord(Tx=(-80, -45) * u.arcsec, Ty=(-115, -140) * u.arcsec, frame=ref_map.coordinate_frame)
qs_coords2 = SkyCoord(Tx=(-70, -30) * u.arcsec, Ty=(-275, -250) * u.arcsec, frame=ref_map.coordinate_frame)
ref_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
ref_map.draw_quadrangle(qs_coords1,axes=ax,edgecolor="blue",linestyle="-",linewidth=1,label='QS sample 1',alpha=0.5)
ref_map.draw_quadrangle(qs_coords2,axes=ax,edgecolor="yellow",linestyle="-",linewidth=1,label='QS sample 2', alpha=0.5)
plt.colorbar()
plt.title('IRIS SJI 2796')
plt.legend()
plt.show()

