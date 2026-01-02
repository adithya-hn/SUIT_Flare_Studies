
'''
Created on 4 Nov 2025
Author: adithya-hn
Purpose: Modular version of overplot contours of SUIT Mg II k  images to AIA


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
import seaborn as sns
from skimage import filters, measure
from skimage.morphology import disk, closing, opening,remove_small_objects

import sys
sys.path.append('/home/adithya/Adithya_repos/pil_detection/utils')
sys.path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()
# Import local libraries
from region_detection import pos_neg_detection
from pil_detection import detection
from video_loading import video_loader

import warnings
import logging

warnings.simplefilter('ignore')
log_ = logging.getLogger('sunpy')
log_.setLevel('WARNING')
logging.getLogger('reproject').setLevel(logging.WARNING)

#------------------------------------------------------------------------------#


def get_suit_scale_rebined_map(ref_fd_1600,ref_mg_rot):
    scale=ref_fd_1600.scale[0].value/ref_mg_rot.scale[0].value
    fd_new_dem=[ref_fd_1600.data.shape[1]*scale,ref_fd_1600.data.shape[0]*scale]*u.pixel
    ref_aia_resmp=ref_fd_1600.resample(fd_new_dem)
    blo = SkyCoord(ref_mg_rot.bottom_left_coord.Tx, ref_mg_rot.bottom_left_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    tro = SkyCoord(ref_mg_rot.top_right_coord.Tx, ref_mg_rot.top_right_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    return ref_aia_resmp.submap(blo,top_right=tro)  

def rebin_suit_map(aia_map,ref_mg_rot):
    scale=ref_mg_rot.scale[0].value/aia_map.scale[0].value
    fd_new_dem=[aia_map.data.shape[1]*scale,aia_map.data.shape[0]*scale]*u.pixel
    ref_aia_resmp=aia_map.resample(fd_new_dem)
    blo = SkyCoord(ref_mg_rot.bottom_left_coord.Tx, ref_mg_rot.bottom_left_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    tro = SkyCoord(ref_mg_rot.top_right_coord.Tx, ref_mg_rot.top_right_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    return ref_aia_resmp.submap(blo,top_right=tro) 

def draw_suit_contours_on_sdo(suitMap,aia_map,base_map,thresh_sig,fltr,aia_dt,base_channel,drot_pil):

    pathlib.Path(jpg_fold+f'/{fltr}_{base_channel}').mkdir(parents=True, exist_ok=True)

    diff_img=(suitMap.data*1000/int(suitMap.meta.get('CMD_EXPT')))-base_map.data
    hist_data=np.array(diff_img.flatten(),dtype='int')
    median_val=np.median(hist_data)
    mad_sig=np.median(np.abs(hist_data-np.median(hist_data)))/0.6745
    v_Thresh=median_val+thresh_sig*mad_sig
    img=np.where(diff_img>v_Thresh,diff_img,0)
    binary_image = diff_img > v_Thresh# True where pixel value > threshold
    labels = measure.label(binary_image,connectivity=2)
    if labels.max()>1:
        cleaned = remove_small_objects(labels, min_size=16)
    else:
        cleaned=labels
    mask=cleaned>0
    th_diff_img=mask.astype(int)*diff_img
    norm_mg_Map=Map(th_diff_img,suitMap.fits_header)
    mask = np.isnan(drot_pil.data) | (drot_pil.data < 0.1)
    masked_pil = np.ma.array(drot_pil.data, mask=mask)
    masked_pil.data[~masked_pil.mask] = 1.0
    mask_pil_map=Map(masked_pil,drot_pil.meta)

    fl_nm=jpg_fold+f'/{fltr}_{base_channel}'+'/'+os.path.basename(suitMap.meta['F_NAME'])[:-4]+'jpg'
    
    fig=plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection=aia_map)
    sns_cl = sns.color_palette("coolwarm",as_cmap=True)
    im=ax.imshow(aia_map.data,cmap=sns_cl,alpha=1,vmin=-2000,vmax=2000)
    #aia_map.plot(axes=ax,cmap='gray',title=str(aia_dt))
    norm_mg_Map.draw_contours(axes=ax, levels=v_Thresh,lws=0.5,colors=['green'],alpha=0.5)
    mask_pil_map.plot(axes=ax)
    ax.set_title(str(aia_dt))
    # ax.set_xlim(100,520)
    # ax.set_ylim(180,520)
    plt.savefig(fl_nm,dpi=300)
    plt.close()  
    return  




if __name__=='__main__':

#-----------------Intial paths and params--------------------------

    suit_aligned_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/aligned_crop/'
    aia_imgs='/media/adithya/Adi_disk4/SUIT_flare_work/case4_jul10/data/aia/cut_outs/171_cutouts/'
    hmi_imgs   = '/media/adithya/Adi_disk4/SUIT_flare_work/case4_jul10/data/HMI/HMI_cutouts/'
    #-----
    
    fol_nm=os.getcwd() #Custom folder to save contour images

    # # Create class object by initializing path to samples and the HARP number of the sample input
    dt = pos_neg_detection()
    pil_dt = detection(hmi_imgs)
    #pil_dt.check_outpath(output_path)
    STRENGTH_FILTER = 100
    BUFFER_SIZE = 4
    PRESERVED_FLUX_RATIO = 0.95
    MIN_MPIL_SIZE = 14
    pil_orig, label_orig = pil_dt.PIL_detect(pos_gauss = STRENGTH_FILTER, neg_gauss= -STRENGTH_FILTER, size_kernel = BUFFER_SIZE)
    
    jpg_fold=fol_nm+'/'+'Contour_imgs'
    fltr='NB04'
    thresh_sig=5
    base_channel='HMI'
    save_aligned_fits='yes'
    save_fits='no'
    save_pngs='yes'     #aligned pngs
    draw_contours='yes'
    
    nb4_fls = glob.glob(suit_aligned_files + '*3'+f'{fltr}.fits')
    #aia_fls = glob.glob(aia_imgs + '*.fits')
    hmi_fls = glob.glob(hmi_imgs + '*.fits')
    

    print('No of SUIT files:',len(nb4_fls))
    #print('No of AIA images: ',len(aia_fls))
    print('No of HMI images: ',len(hmi_fls))

    
    #rebin aia image
    
    suit_sq=Map(nb4_fls,sequence=True)
    ref_suit_map=Map(suit_sq[0])

    data_stack = np.stack([(suit_sq[i].data*1000/suit_sq[i].meta.get('CMD_EXPT'))for i in range(5)])
    base_img=np.median(data_stack, axis=0)
    base_map=Map(base_img,suit_sq[0].meta)
    

    aia_dt=[]
    hmi_dt=[]

    # for j in range(len(aia_fls)):
    #     aia_map=Map(aia_fls[j])
    #     aia_dt.append(aia_map.date.datetime)
    # aia_dt_array=Time(parse_time(aia_dt))

    for j in range(len(hmi_fls)):
        hmi_map=Map(hmi_fls[j])
        hmi_dt.append(hmi_map.date.datetime)
    hmi_dt_array=Time(parse_time(hmi_dt))

    for i in tqdm (range(len(nb4_fls))):
        suit_map=suit_sq[i]
        base_time=Time(parse_time(suit_map.date))
        # idx1=np.argmin(np.abs(aia_dt_array - base_time))
        idx=np.argmin(np.abs(hmi_dt_array - base_time))
        hmi_map_=Map(hmi_fls[idx])
        data_num = pil_dt.check_header(hmi_map_)
        masked_pil = dt.mask_pil(label_orig[data_num])
        strength_label = pil_dt.filter_by_strength(threshold = PRESERVED_FLUX_RATIO)
        masked_filter_RoPIs = dt.mask_pil(strength_label[data_num])
        thin_df, thin_binary = pil_dt.thin_strength_label(strength_label)
        masked_thinned_MPIL = dt.mask_pil(thin_binary[data_num])
        pil_msk=Map(masked_thinned_MPIL,hmi_map_.meta)
        
        hmi_map_dt=hmi_map_.date
        if "CROTA2" in hmi_map_.meta:
            map_rot_angl=int(hmi_map_.meta.get('CROTA2'))
            if map_rot_angl>5:
                hmi_map=hmi_map_.rotate(angle=-map_rot_angl*u.deg)
                pil_msk=pil_msk.rotate(angle=-map_rot_angl*u.deg)
            #hmi_dt=hmi_map.date
            
        else:
            hmi_map=hmi_map_
        if i==0:
            #ref_aia_map=aia_map#Map(aia_fls[idx1])
            ref_hmi_map=hmi_map#Map(hmi_fls[idx])
        # fig = plt.figure(figsize=(8,6))
        # hmi_map.plot()
        # pil_msk.plot()
        # plt.imshow(masked_thinned_MPIL, 'spring', interpolation='none', alpha=1) # 'spring' represents the bright pink color
        # plt.xlabel('Carrington Longitude [deg]',fontsize = 12)
        # plt.ylabel('Latitude [deg]',fontsize = 12)
        # plt.show()
        #----
        

        if base_channel=='HMI':
            hmi_rebinned=get_suit_scale_rebined_map(hmi_map,ref_suit_map)
            pil_rebin=get_suit_scale_rebined_map(pil_msk,ref_suit_map)
            with propagate_with_solar_surface():
                hmi_map_drot=(hmi_rebinned.reproject_to(ref_hmi_map.wcs,dask_method='none'))
                drot_pil=(pil_rebin.reproject_to(ref_hmi_map.wcs,dask_method='none'))
            draw_suit_contours_on_sdo(suit_map,hmi_map_drot,base_map,thresh_sig,fltr,hmi_map_dt,base_channel,drot_pil)

        # if base_channel=='171':
        #     aia_rebinned=get_suit_scale_rebined_map(aia_map,ref_suit_map)
        #     with propagate_with_solar_surface():
        #         aia_map_drot=(aia_rebinned.reproject_to(ref_aia_map.wcs,dask_method='none'))
        #     draw_suit_contours_on_sdo(suit_map,aia_map_drot,base_map,thresh_sig,fltr,aia_dt,base_channel)
    

