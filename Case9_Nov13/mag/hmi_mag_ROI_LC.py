
#==========================================
#Differential rotation implimented
#==========================================

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
from datetime import timedelta
import timeit
import pathlib
from astropy.coordinates import SkyCoord
#from colormap import filterColor
import numpy as np
from sunpy.coordinates import RotatedSunFrame
from tqdm import tqdm
import matplotlib.dates as mdates

start = timeit.default_timer()
now = datetime.datetime.now()-timedelta(days=1)
fol_nm='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case9_Nov13/mag'
fol_nm=fol_nm+'/'+'Imgs_fm'
pathlib.Path(fol_nm).mkdir(parents=True, exist_ok=True)

search_fold='/media/adithya/Adi_disk4/SUIT_flare_work/case9_nov13/data/HMI/cut_outs/HMI_cutouts/'

fdir =os.listdir(search_fold)
fdir.sort()
#print(fdir)
key=['magnetogram']
param=key[0]

plot_data=[]
dates=[]


box_pth=fol_nm+'/'+param
pathlib.Path(box_pth+'/Box').mkdir(parents=True, exist_ok=True)
#pathlib.Path(box_pth).mkdir(parents=True, exist_ok=True)
#print(fdir+fltr+'/'+'*3'+fltr+'.fits')
files = sorted(glob.glob(f'{search_fold}*.fits'))
print('Total files:',len(files))
#print(files)
Sequence = sunpy.map.Map(files, sequence=True)
fltr_count=[]
date_array=[]
for i in tqdm(range (len(Sequence))):
    
    hmi_map=Sequence[i]
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(projection=hmi_map)
    hmi_map.plot(axes=ax)
    #hmi_map.peek()diffrot_point2
    #coords = SkyCoord(Tx=(-520, -390) * u.arcsec,Ty=(-323, -260) * u.arcsec,frame=hmi_map.coordinate_frame,)
    #coords = SkyCoord(lat=(-21,-16)  * u.deg,lon=(348, 354)* u.deg,frame=hmi_map.coordinate_frame,)

    ###point1 = SkyCoord(-345*u.arcsec, -260*u.arcsec, frame=hmi_map.coordinate_frame) # ca second box
    #point1 = SkyCoord(-466*u.arcsec, -255*u.arcsec, frame=hmi_map.coordinate_frame) #ca bright box
    #point1 = SkyCoord(-500*u.arcsec, -312*u.arcsec, frame=hmi_map.coordinate_frame) #flare core
    point1 = SkyCoord(-160*u.arcsec, -280*u.arcsec, frame=hmi_map.coordinate_frame) #flull region
    
    #point1 = SkyCoord(-500*u.arcsec, -323*u.arcsec, frame=hmi_map.coordinate_frame) #wide box
    diffrot_point1 = SkyCoord(RotatedSunFrame(base=point1, duration=45*(i+1)*u.second))
    transformed_diffrot_point1 = diffrot_point1.transform_to(hmi_map.coordinate_frame)

    ###point2 = SkyCoord(-300*u.arcsec, -200*u.arcsec, frame=hmi_map.coordinate_frame) #
    #point2 = SkyCoord(-427*u.arcsec, -191*u.arcsec, frame=hmi_map.coordinate_frame) #
    #point2 = SkyCoord(-400*u.arcsec, -260*u.arcsec, frame=hmi_map.coordinate_frame) #small box
    point2 = SkyCoord(120*u.arcsec, -150*u.arcsec, frame=hmi_map.coordinate_frame)

    #point2 = SkyCoord(-375*u.arcsec, -260*u.arcsec, frame=hmi_map.coordinate_frame) #wide box
    diffrot_point2 = SkyCoord(RotatedSunFrame(base=point2, duration=45*(i+1)*u.second))
    transformed_diffrot_point2 = diffrot_point2.transform_to(hmi_map.coordinate_frame)

    #print(transformed_diffrot_point1.Tx.value,transformed_diffrot_point2.Tx.value)
    coords = SkyCoord(Tx=(transformed_diffrot_point1.Tx.value, transformed_diffrot_point2.Tx.value) * u.arcsec,Ty=(transformed_diffrot_point1.Ty.value, transformed_diffrot_point2.Ty.value) * u.arcsec,frame=hmi_map.coordinate_frame,)
    hmi_map.draw_quadrangle(coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='2-element SkyCoord')
    fnm=(os.path.basename(files[i]))[:-5]
    hmi_box=hmi_map.submap(coords)
    fltr_count.append(np.sum(abs(hmi_box.data)))
    date_array.append(hmi_box.date.datetime)

    #print('-------',f'{fol_nm}/{fnm}.jpg')
    plt.savefig(f'{box_pth}/{fnm}.jpg',dpi=300)
    plt.close()
    #ig = plt.figure(figsize=(5, 5))
    #suit_box=hmi_map.submap(coords)
    #suit_box.save(f'/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Aligned_HMI/{fnm}.fits',overwrite=True)
plt.figure()
plt.plot(date_array,fltr_count)
plt.xlabel('Time')
plt.ylabel('HMI count')
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig('hmi_box.png',dpi=300)
plt.show()