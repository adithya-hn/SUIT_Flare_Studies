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
from sunpy.time import parse_time
from astropy.time import Time

#Threshold values:

nb3T=11000
nb4T=11500
nb8T=3900
nb6T=86000
nb7T=295000
nb3Mx=14000
nb4Mx=15000
nb8Mx=4300
nb6Mx=95000
nb7Mx=310000
Filters=['NB08']
fltr2='NB03'
for fltr in Filters:
    plot_data=[]
    tot_count=[]
    dates=[]

    search_fold=f'/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/{fltr}/' #Custom Folder
    base_fold=f'/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/{fltr2}/'
    
    print(f'Searching for {fltr} images in {search_fold} folder')
    fdir =search_fold 
    files = glob.glob(fdir + '*3'+fltr+'.fits')
    files=sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    files=files
    b_files=glob.glob(base_fold + '*3'+fltr2+'.fits')
    b_files=sorted(b_files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

    print('Total files:',len(files))
    base_time_array=[]
    for f in range(len(b_files)):
        #print(datetime.datetime.strptime(os.path.basename(b_files[f]).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
        base_time_array.append(datetime.datetime.strptime(os.path.basename(b_files[f]).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    base_time_array=Time(parse_time(base_time_array)) #np.array(base_time_array)

    fol_nm=os.getcwd() #str(fl_date.day).zfill(2)+'_'+str(fl_date.month).zfill(2)+'_'+str(fl_date.year).zfill(2)
    print(fol_nm)

    jpg_fold=fol_nm+'/'+'Threshold_imgs'

    #pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)
    pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
    pathlib.Path(jpg_fold+f'/{fltr}').mkdir(parents=True, exist_ok=True)

    sigma=1
    if fltr=='NB03':
        Thresh_val=nb3T
        Tmax=nb3Mx
    
    if fltr=='NB04':
        Thresh_val=nb4T
        Tmax=nb4Mx

    if fltr=='NB08':
        Thresh_val=nb8T
        Tmax=nb8Mx

    if fltr=='NB06':
        Thresh_val=nb6T
        Tmax=nb6Mx
    
    if fltr=='NB07':
        Thresh_val=nb7T
        Tmax=nb7Mx
    
    for i in range(len(files)):
        suitMap=sunpy.map.Map(files[i])
        base_time=Time(parse_time(suitMap.date))
        idx=np.argmin(np.abs(base_time_array - base_time))

        #print('--->',np.argmin(np.abs(base_time_array - base_time)))#,base_time_array[np.argmin(np.abs(base_time_array - base_time))],base_time)
        BaseMap=sunpy.map.Map(b_files[idx])
        base_data=BaseMap.data*1000/int(BaseMap.meta.get('MEAS_EXP'))
        Base_img=sunpy.map.Map(base_data,BaseMap.fits_header)

        img_head=suitMap.fits_header
        suit_data=suitMap.data#gaussian_filter(suitMap.data,sigma=sigma)
        alned_data=suit_data*1000/int(suitMap.meta.get('MEAS_EXP'))
        th_lvs=[4100,4300]
        th_lvs2=[12000,13000]
        
        Thresh_alned_data=np.where(alned_data>Thresh_val,alned_data,0)
        alignedMap=sunpy.map.Map(Thresh_alned_data,img_head)

        fl_nm=jpg_fold+f'/{fltr}'+'/'+os.path.basename(files[i])[:-4]+'jpg'
        fig=plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection=alignedMap)

        Base_img.plot(cmap='gray',autoalign=True, vmin=2000, vmax=14000)

        alignedMap.draw_contours(axes=ax, levels=th_lvs,zorder=1,colors=['blue','yellow'],alpha=0.9)
        Base_img.draw_contours(axes=ax,levels=th_lvs2,colors=['red', 'green'],alpha=0.7,zorder=2)
        plot_data.append(np.count_nonzero(Thresh_alned_data))
        tot_count.append(np.sum(Thresh_alned_data))
        dates.append(suitMap.date)
        #alignedMap.draw_limb(axes=ax)
        #alignedMap.draw_grid(axes=ax)
        plot_str='Ca II h: '+str(base_time_array[idx])
        ax.text(50,50, plot_str, color='white', fontsize=10)
        plt.draw()
        plt.colorbar()
        #plt.title('NB03 diff: '+str(suitMap_.date))
        #plt.axis('off')
        #plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        #plt.tight_layout()
        plt.savefig(fl_nm)
        plt.close()    
        #print(i,' / ',len(files))
    plot_data=np.array(plot_data)
    tot_count=np.array(tot_count)
    np.savetxt(f'{fltr}_threshold_count.csv',np.c_[dates,plot_data,tot_count],delimiter=',',fmt='%s')



