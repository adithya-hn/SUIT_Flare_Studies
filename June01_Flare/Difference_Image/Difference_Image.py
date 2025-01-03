import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
from sunkit_image.coalignment import calculate_match_template_shift,apply_shifts
from datetime import timedelta
import timeit
import pathlib
from colormap import filterColor
from astropy.coordinates import SkyCoord
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
#import ImagesToMovie_pkg
import matplotlib.image as mpimg
from PIL import Image
import pandas as pd
from scipy.ndimage import shift
from scipy.ndimage import gaussian_filter

Filter='NB08'
search_fold=f'/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Processed_Data/Aligned_images/{Filter}/' #Custom Folder

d=6
print(f'Searching for {Filter} images in {search_fold} folder')
fdir =search_fold 
files = glob.glob(fdir + '*3'+Filter+'.fits')
files=sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
files=files
print('Total files:',len(files))

ref_img=sunpy.map.Map(files[0])
#base_fold='/home1/Data/Adithya/POC_Works/Jitter/'
fl_date = datetime.datetime.fromisoformat(str(ref_img.date))
fol_nm=str(fl_date.day).zfill(2)+'_'+str(fl_date.month).zfill(2)+'_'+str(fl_date.year).zfill(2)

jpg_fold=fol_nm+'/'+'FD_imgs'
algn_dir=fol_nm+'/'+'Aligned_Fits'

#pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold+f'/{Filter}').mkdir(parents=True, exist_ok=True)

sigma=1
box_count=[]
full_img_count=[]
er_box_count=[]
date_array=[]
for i in range(5,len(files)):
    suitMap=sunpy.map.Map(files[4])
    suitMap_=sunpy.map.Map(files[i])
    img_head=suitMap_.fits_header
    
    suit_data=gaussian_filter(1000*suitMap.data/suitMap.meta.get('MEAS_EXP'),sigma=sigma)
    suit_data_=gaussian_filter(1000*suitMap_.data/suitMap_.meta.get('MEAS_EXP'),sigma=sigma)
    coords = SkyCoord(Tx=(-650, -300) * u.arcsec,Ty=(-140, -420) * u.arcsec,frame=suitMap_.coordinate_frame,)
    
    

    #alignedMap=apply_shifts(suitMap,yshift=100*u.pixel,xshift=100*u.pixel,clip=False)
    alned_data=suit_data
    alned_data_=suit_data_
    diff_img=alned_data_-alned_data #n -n-5
    alignedMap_=sunpy.map.Map(diff_img,img_head)
    alignedMap=alignedMap_.submap(coords)
    fl_nm=jpg_fold+f'/{Filter}'+'/'+os.path.basename(files[i+d])[:-4]+'jpg'
    fig=plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, projection=alignedMap)
    #initial_point = SkyCoord(-230*u.arcsec, -510*u.arcsec, frame=alignedMap.coordinate_frame)
    #point_plot = ax.scatter(initial_point.Tx.value, initial_point.Ty.value, color='white', s=90,marker='+')
    
    
    #'''
    box_coords = SkyCoord(Tx=(-580, -390) * u.arcsec, Ty=(-160, -320) * u.arcsec, frame=alignedMap.coordinate_frame)
    alignedMap.draw_quadrangle(box_coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
    er_coords = SkyCoord(Tx=(-400, -310) * u.arcsec, Ty=(-325, -415) * u.arcsec, frame=alignedMap.coordinate_frame)
    alignedMap.draw_quadrangle(er_coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')

    suit_box=alignedMap.submap(box_coords)    
    er_box=alignedMap.submap(er_coords)
    er_Main_img=suitMap_.submap(er_coords)

    box_count.append(np.sum(suit_box.data))
    full_img_count.append(np.sum(alignedMap.data))
    l,h=np.shape(er_box.data)
    tot_er_npix=l*h
    #print(np.std(er_Main_img.data))
    er_box_count.append(np.sum(abs(er_box.data)))
    date_array.append(suitMap.date)
    #'''

    #alignedMap.draw_limb(axes=ax)

    #alignedMap.draw_grid(axes=ax)
    alignedMap.plot(cmap=filterColor[Filter], vmin=0,vmax=300) #nb08 =0.5,if not 2
    ax.plot_coord(SkyCoord(-510*u.arcsec, -215*u.arcsec, frame=alignedMap.coordinate_frame),color='white',marker='X',markersize=12)
    plot_str=str(alignedMap.date)+' - '+str(suitMap.date)
    ax.text(50,50, plot_str, color='white', fontsize=10)
    #ax.legend(loc="lower center")

    plt.xticks([])
    plt.yticks([])
    fig.canvas.draw()
    plt.draw()
    plt.colorbar()
    #plt.title('NB03 diff: '+str(suitMap_.date))
    #plt.axis('off')
    #plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    
    plt.tight_layout()
    plt.savefig(fl_nm)

    plt.close()    
    print(i,' / ',len(files))

np.savetxt(f'Diff_img_data_{Filter}.dat',np.c_[date_array,box_count,full_img_count,er_box_count],header='Date  Box_count   Full_image_count     Error',fmt='%s')

