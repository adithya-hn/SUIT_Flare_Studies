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
from astropy.coordinates import SkyCoord
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
#import ImagesToMovie_pkg
import matplotlib.image as mpimg
from PIL import Image
import pandas as pd
from scipy.ndimage import shift
from scipy.ndimage import gaussian_filter
from colormap import filterColor
import map_to_image

Filter='NB08'
search_fold=f'/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/july10_11ut_to_18ut/Processed/Aligned_images/{Filter}/' #Custom Folder

d=2
print(f'Searching for {Filter} images in {search_fold} folder')
fdir =search_fold 
files = glob.glob(fdir + '*3'+Filter+'.fits')
files=sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
print('Total files:',len(files))

ref_img=sunpy.map.Map(files[0])
coords_ = SkyCoord(Tx=(-225, -25) * u.arcsec, Ty=(-300, -100)  * u.arcsec,frame=ref_img.coordinate_frame)
ref_img=ref_img.submap(coords_)
#ref_img.peek()
fl_date = datetime.datetime.fromisoformat(str(ref_img.date))
fol_nm=str(fl_date.day).zfill(2)+'_'+str(fl_date.month).zfill(2)+'_'+str(fl_date.year).zfill(2)

jpg_fold=fol_nm+'/'+'FD_imgs'
algn_dir=fol_nm+'/'+'Aligned_Fits'

pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold+f'/{Filter}').mkdir(parents=True, exist_ok=True)


sigma=1
for i in range(len(files)-d):
    suitMap=sunpy.map.Map(files[i])
    suitMap_=sunpy.map.Map(files[i+d])
    
    pathlib.Path(f'Flare_imgs/{Filter}').mkdir(parents=True, exist_ok=True)
    fl_pth=f'Flare_imgs/{Filter}/{(os.path.basename(files[i+d]))[:-4]}jpg'

    coords = SkyCoord(Tx=(-225, -25) * u.arcsec, Ty=(-300, -100) * u.arcsec,frame=suitMap.coordinate_frame,)
    diff_img=gaussian_filter((suitMap_.data*1000/suitMap_.meta.get('CMD_EXPT')),sigma=2)-gaussian_filter((suitMap.data*1000/suitMap.meta.get('CMD_EXPT')),sigma=2)
    map_to_image.sunpy_map_to_jpg(suitMap_,fl_pth)
    diff_map=sunpy.map.Map(diff_img,suitMap.fits_header)
    Diff_Map=diff_map.submap(coords)
    fnm_=os.path.basename(files[i+d])
    fnm='Diff_'+fnm_
    fl_nm=jpg_fold+f'/{Filter}'+'/'+os.path.basename(files[i+d])[:-4]+'jpg'
    fig=plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, projection=Diff_Map)
    Diff_Map.plot(cmap='hmimag',vmin=-200,vmax=200) #nb08 =0.5
    ax.plot_coord(SkyCoord(-125*u.arcsec, -180*u.arcsec, frame=Diff_Map.coordinate_frame),color='yellow',marker='X',markersize=12)
    #Diff_Map.draw_limb(axes=ax)
    #alignedMap.draw_grid(axes=ax)
    plot_str=str(Diff_Map.date)+' - '+str(suitMap.date)
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
    print(i,' / ',len(files) )



