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
from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord
from skimage import measure
from matplotlib.path import Path
from astropy.time import Time, TimeDelta
from astropy.table import Table

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()
#----------------

# Set the path to the directory where the files will be saved
from astropy.io import fits

# Path to your SJI FITS file
filename = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/IRIS/data/raw/pre_flare_data/iris_l2_20240710_141828_3620106457_SJI_2796_t000.fits"

# Open the FITS file
hdulist = fits.open(filename)

# Each extension (starting from 1) contains an image and header with timestamp
timestamps = []
hdu =hdulist # skip primary header [0]
table_data = hdu[1].data  # Assuming the table is in HDU 1
print(len(table_data))
# Extract the timestamps from the header
table=Table.read(filename) 

#print(table.columns)

path = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/IRIS/results/pre-flare_iris_lc/'
#
# Create the directory if it doesn't exist
if not os.path.exists(path):
    os.makedirs(path)


# Load the IRIS SJI 2796 map
iris_img = fits.open("../data/raw/pre_flare_data/iris_l2_20240710_141828_3620106457_SJI_2796_t000.fits")

sji_imgs=iris_img[0].data
iris_ref_map=sunpy.map.Map(sji_imgs[21], iris_img[0].header) #SJI flare map
base_time = iris_img[0].header['STARTOBS']
times=hdulist[1].data[:,0]
t0 = Time(base_time, format='isot', scale='utc')
actual_times = t0 + TimeDelta(times, format='sec')
date_array_=[ t.iso for t in actual_times]
date_array = [datetime.datetime.strptime(ts, "%Y-%m-%d %H:%M:%S.%f") for ts in date_array_]

#print('Start time:', date_array)

#SUIT plate scale correction
suit_ref_map = sunpy.map.Map("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/raw/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.35.28.135_0983NB03.fits")
suit_pos = get_horizons_coord(-21, suit_ref_map.date)
suit_ref_map.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))

# QS region for threshold values
qs_coord=SkyCoord(Tx=(-10, 90) * u.arcsec, Ty=(-55, -110) * u.arcsec, frame=suit_ref_map.coordinate_frame)
qs_map=suit_ref_map.submap(qs_coord)
qs_data = qs_map.data * 1000 / qs_map.meta.get('CMD_EXPT')
print(np.median(qs_data),np.mean(qs_data),np.std(qs_data))
Thresh_val=np.median(qs_data)*3
print('Threshold: ',Thresh_val)

qs_coords1 = SkyCoord(Tx=(-80, -45) * u.arcsec, Ty=(-115, -140) * u.arcsec, frame=iris_ref_map.coordinate_frame)
qs_coords2 = SkyCoord(Tx=(-70, -35) * u.arcsec, Ty=(-275, -250) * u.arcsec, frame=iris_ref_map.coordinate_frame)

normalized_data = suit_ref_map.data * 1000 / suit_ref_map.meta.get('CMD_EXPT')
threshold = Thresh_val  # adjust as needed
contours_pix = measure.find_contours(normalized_data, threshold)
largest_contour = max(contours_pix, key=len) #finding the largest contour

# Create a mask for the contour

# Converting to Helip-projective coordinates
hpc_coords2 = suit_ref_map.pixel_to_world(largest_contour[:, 1]*u.pixel, largest_contour[:, 0]*u.pixel)
hpc_coords2= hpc_coords2.transform_to(iris_ref_map.coordinate_frame) # not much chage

x_pix, y_pix = iris_ref_map.world_to_pixel(hpc_coords2)
'''
plt.imshow(iris_ref_map.data, cmap='gray', origin='lower')
plt.plot(x_pix, y_pix, color='red', linewidth=1,label='Flare contour')
plt.colorbar()
plt.show()'''

# Create polygon path from contour in pixel space
contour_path = Path(np.vstack([x_pix.value, y_pix.value]).T)

# Create grid of all pixel positions
ny, nx = iris_ref_map.data.shape
X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
points = np.vstack((X.ravel(), Y.ravel())).T

# Create mask
mask = contour_path.contains_points(points).reshape((ny, nx))

maps=[]
flare_ct=[]
qs_count1=[]
qs_count2=[]
#time_range = TimeRange(iris_ref_map.date - timedelta(minutes=5), iris_ref_map.date + timedelta(minutes=5))
for i in range(111):
    maps.append(sunpy.map.Map(sji_imgs[i], iris_img[0].header))
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(projection=iris_ref_map)
    maps[i].plot(axes=ax, vmin=0,vmax=2500,alpha=1)
    #suit_ref_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent,alpha=0.8)
    #suit_ref_map.draw_quadrangle(qs_coord,axes=ax,edgecolor="blue",linestyle="-",linewidth=1,label='QS region',alpha=0.5)
    qs_map1=maps[i].submap(qs_coords1)
    qs_map2=maps[i].submap(qs_coords2)
    print(qs_map1.meta.get('CMD_EXPT'))
    qs_data1 = qs_map1.data 
    qs_data2 = qs_map2.data 
    qs_count1.append(np.mean(qs_data1))
    qs_count2.append(np.mean(qs_data2))
    iris_ref_map.draw_quadrangle(qs_coords1,axes=ax,edgecolor="blue",linestyle="-",linewidth=1,label='QS sample 1',alpha=0.5)
    iris_ref_map.draw_quadrangle(qs_coords2,axes=ax,edgecolor="yellow",linestyle="-",linewidth=1,label='QS sample 2', alpha=0.5)
    counts = maps[i].data[mask].sum() # Get counts under the mask
    flare_ct.append(counts)
    print(f"Counts in Filter 1 under the contour: {counts}")
    ax.plot_coord(hpc_coords2, color='red', linewidth=1,label='Flare contour')
    plt.colorbar()
    plt.title(f'IRIS SJI 2796: {date_array[i]}')
    plt.legend()
    plt.savefig(f"../results/pre-flare_iris_lc/iris_sji_2796_{i}.png", dpi=300)
    plt.close()

flare_ct = np.array(flare_ct)
date_array = np.array(date_array)
qs_count1 = np.array(qs_count1)
qs_count2 = np.array(qs_count2)
np.savetxt(f"../results/pre_flare_iris_sji_2796_lc.csv", np.c_[date_array[:111], flare_ct,qs_count1,qs_count2], delimiter=',', header='Date,Counts', fmt='%s')
# Save the light curve data to a CSV file
#print(len(flare_ct),date_array[:111])
# Plotting the light curve
plt.figure(figsize=(10, 5))
plt.plot( date_array[:111], flare_ct, marker='o', linestyle='-',markersize=1,linewidth=0.5, color='blue')
plt.xticks(rotation=45)
plt.tight_layout()
plt.xlabel('Time (UTC)')
plt.ylabel('Counts')
plt.title('IRIS SJI 2796 Light Curve')
plt.xticks(rotation=45)
plt.savefig(f"../results/pre-flare_iris_lc/iris_sji_2796_light_curve.png", dpi=300)
plt.show()
#/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/IRIS/results/pre-flare_iris_lc