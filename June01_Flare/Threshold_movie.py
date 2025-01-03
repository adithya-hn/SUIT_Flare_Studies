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

search_fold='/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Processed_Data/Aligned_images/NB03/' #Custom Folder

Filter='NB03'
d=2
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

jpg_fold=fol_nm+'/'+'Threshold_imgs'
algn_dir=fol_nm+'/'+'Aligned_Fits'

#pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold+f'/{Filter}').mkdir(parents=True, exist_ok=True)

sigma=1
Thresh_val=16
for i in range(len(files)):
    suitMap=sunpy.map.Map(files[i])
    img_head=suitMap.fits_header
    suit_data=suitMap.data#gaussian_filter(suitMap.data,sigma=sigma)
    alned_data=suit_data/suitMap.meta.get('MEAS_EXP')
    Thresh_alned_data=np.where(alned_data>Thresh_val,alned_data,0)
    alignedMap=sunpy.map.Map(Thresh_alned_data,img_head)
    fl_nm=jpg_fold+f'/{Filter}'+'/'+os.path.basename(files[i])[:-4]+'jpg'
    fig=plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection=alignedMap)
    alignedMap.plot(cmap=filterColor['NB03'], vmin=Thresh_val-1, vmax=np.percentile(alignedMap.data,99.9))
    #alignedMap.draw_limb(axes=ax)
    #alignedMap.draw_grid(axes=ax)
    plot_str=str(alignedMap.date)+' - '+str(suitMap.date)
    ax.text(50,50, plot_str, color='white', fontsize=10)
    plt.draw()
    plt.colorbar()
    #plt.title('NB03 diff: '+str(suitMap_.date))
    #plt.axis('off')
    #plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    #plt.tight_layout()
    plt.savefig(fl_nm)
    plt.close()    
    print(i,' / ',len(files))



