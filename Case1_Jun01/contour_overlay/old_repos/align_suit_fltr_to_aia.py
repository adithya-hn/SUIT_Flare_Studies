
'''
Created on 25 oct 2025
Author: adithya-hn
Purpose: Modular version of overplot contours of SUIT Mg II k and Ca II H on AIA 1600 images and align SUIT images to AIA


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
def get_sunspot_threshold(map):
    valid = map.data[(map.data > 100)]
    valid_int = valid.astype(int)  
    mode_val = stats.mode(valid_int, keepdims=True).mode[0]/map.meta.get('CMD_EXPT')*1000  #normalised mode value    
    th_lvs=np.array([mode_val*0.7,mode_val*1.5,mode_val*3] )
    return mode_val,th_lvs

def get_mg_threshold(ref_mg_rot):
    valid = ref_mg_rot.data[(ref_mg_rot.data > 100)]
    valid_int = valid.astype(int)  
    mode_val = stats.mode(valid_int, keepdims=True).mode[0]/ref_mg_rot.meta.get('CMD_EXPT')*1000  #normalised mode value    
    th_lvs=np.array([mode_val*0.7,mode_val*1.5,mode_val*3] )
    return mode_val,th_lvs

def get_ca_threshold(samp_ca_map_rot):
    valid_ca = samp_ca_map_rot.data[(samp_ca_map_rot.data > 100)]
    valid_ca_int = valid_ca.astype(int)  
    ca_mode_val = stats.mode(valid_ca_int, keepdims=True).mode[0]/samp_ca_map_rot.meta.get('CMD_EXPT')*1000  #normalised mode value    
    th_lvs2=np.array([ca_mode_val*0.7,ca_mode_val*1.03,ca_mode_val*1.15])
    return ca_mode_val,th_lvs2

def get_rebinned_crop_map(ref_fd_1600,ref_mg_rot):
    scale=ref_fd_1600.scale[0].value/ref_mg_rot.scale[0].value
    fd_new_dem=[ref_fd_1600.data.shape[1]*scale,ref_fd_1600.data.shape[0]*scale]*u.pixel
    ref_aia_resmp=ref_fd_1600.resample(fd_new_dem)
    blo = SkyCoord(ref_mg_rot.bottom_left_coord.Tx, ref_mg_rot.bottom_left_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    tro = SkyCoord(ref_mg_rot.top_right_coord.Tx, ref_mg_rot.top_right_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    return ref_aia_resmp.submap(blo,top_right=tro)  

def draw_suit_contours_on_sdo(Map_sq,jpg_fold,fltr,ref_baseMap,th_lvs):

    pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
    pathlib.Path(jpg_fold+f'/{fltr}').mkdir(parents=True, exist_ok=True)

    for i in tqdm (range(len(Map_sq))):
        suitMap=Map_sq[i]

        norm_data=suitMap.data*1000/int(suitMap.meta.get('CMD_EXPT'))
        norm_mg_Map=Map(gaussian_filter(norm_data, sigma=1),suitMap.fits_header)

        fl_nm=jpg_fold+f'/{fltr}'+'/'+os.path.basename(suitMap.meta['F_NAME'])[:-4]+'jpg'
        fig=plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection=ref_baseMap)
        ref_baseMap.plot(axes=ax,cmap='gray',title=str(ref_baseMap.date))
        norm_mg_Map.draw_contours(axes=ax, levels=th_lvs,lws=0.5,colors=['red','magenta','green'])
        #plt.colorbar()
        plt.savefig(fl_nm,dpi=300)
        plt.close()  

    return  print('--------Done----------')

def co_align_maps(suit_raw_files,jpg_fold,ref_fd_img_pth,fltr,nb1_fl,tx1,tx2,ty1,ty2,save_pngs=None,save_fits='yes',draw_contours='yes'):
    
    save_aligned_pth=str(pathlib.Path(suit_raw_files).parents[0])+'/1600_aligned/'  #Custom Folder to save aligned SUIT images

    #---------------aia alignemnt--------
    ref_map=Map(nb1_fl[0]) # or specific map
    suit_pos = get_horizons_coord(-21, ref_map.date)
    ref_map.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
    ref_rot=ref_map.rotate(angle=ref_map.meta["CROTA2"] * u.deg,missing=0)
    ref_fd_img=Map(ref_fd_img_pth)
    if ref_fd_img.meta['CROTA2']>90:
        ref_fd_img=ref_fd_img.rotate(angle=180 * u.deg,missing=0)

    #------get contour thresholds-------
    nb1_mode_val,th_nb1=get_sunspot_threshold(ref_rot)# 
    print(f'Mode value {fltr}: ', nb1_mode_val)
    print(f'{fltr} threshold levels:', th_nb1.astype(int)) 
                                                                                                                                                                                                                                                                                        

    #--AIA/HMI template--
    ref_baseMap=get_rebinned_crop_map(ref_fd_img,ref_rot)
    rect_template=SkyCoord(Tx=(tx1,tx2 )* u.arcsec, Ty=(ty1,ty2 )* u.arcsec, frame=ref_baseMap.coordinate_frame)
    ref_template=ref_baseMap.submap(rect_template)

    #----

    Map_sq=[]
    nb1_map_set=[]

    #--Adding aia base image to sequence for coalignment--

    for i in range(len(nb1_fl)):
        if i==0:
            nb1_map_set.append(ref_baseMap)
        nb1_maps_in=Map(nb1_fl[i])
        nb1_maps_in.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
        nb1_map_set.append(nb1_maps_in.rotate(angle=nb1_maps_in.meta["CROTA2"]*u.deg,missing=0))
    nb1_map_sq=MapSequence(nb1_map_set)

    print('Aligning images...')
    print('Reference image wavelength:',nb1_map_sq[0].meta["WAVELNTH"])

    nb1_aln_maps=mc_coalign(nb1_map_sq,template=ref_template,clip=False)

    print('Done co-aligning')

    #--correcting header--
    for j in range(len(nb1_aln_maps)):
        alnMap=Map(nb1_aln_maps[j].data,nb1_aln_maps[j].meta)
        alnMap.meta['CRPIX1']=nb1_aln_maps[0].meta['CRPIX1']
        alnMap.meta['CRPIX2']=nb1_aln_maps[0].meta['CRPIX2']
        
        if alnMap.meta["WAVELNTH"]==1600:
            continue

        if alnMap.meta["WAVELNTH"]==6173.0:
            continue
        else:
            Map_sq.append(alnMap)
        
        if save_fits=='yes':
            pathlib.Path(save_aligned_pth).mkdir(parents=True, exist_ok=True)
            try:
                fname=str(nb1_aln_maps[j].meta.get('F_NAME'))[:-5]+'.png'
                alnMap_=alnMap.rotate(angle=-Map(nb1_fl[j]).meta['CROTA2']*u.deg,missing=0)
                alnMap_.save(save_aligned_pth+nb1_aln_maps[j].meta.get('F_NAME'),overwrite=True)
            except:
                print('No Mg II file name found, could not save file')
        if save_pngs =='yes':
            try:
                alnMap.plot()
                plt.savefig(fname)
                plt.close()
            except:
                pass
 
    if draw_contours == 'yes':
        draw_suit_contours_on_sdo(Map_sq,jpg_fold,fltr,ref_baseMap,th_nb1)


#----------------------------------------------------------------------------------------



if __name__=='__main__':

#-----------------Intial paths and params--------------------------

    Filters=['6173']
    suit_filters=['NB01','NB02','NB05','NB06','NB07','NB09']
    suit_raw_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/raw/'
    ref_fd_img_pth='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/aia/hmi.ic_45s.20240601_070945_TAI.2.continuum.fits'
    tx1,ty1=-520,-400
    tx2,ty2=-320,-230

    save_aligned_fits='yes'
    save_fits='no'
    save_pngs='yes'     #aligned pngs
    draw_contours='yes'
    fol_nm=os.getcwd() #Custom folder to save contour images
    
    jpg_fold=fol_nm+'/'+'Contour_imgs'
    nb1_fl = glob.glob(suit_raw_files + '*3'+'NB01.fits')
    nb2_fl = glob.glob(suit_raw_files + '*3'+'NB02.fits')
    nb5_fl = glob.glob(suit_raw_files + '*3'+'NB05.fits')
    nb6_fl = glob.glob(suit_raw_files + '*3'+'NB06.fits')
    nb7_fl = glob.glob(suit_raw_files + '*3'+'NB07.fits')

    print('Total SUIT NB01 files:',len(nb1_fl))
    print('Total SUIT NB02 files:',len(nb2_fl))
    print('Total SUIT NB05 files:',len(nb5_fl))
    print('Total SUIT NB06 files:',len(nb6_fl))
    print('Total SUIT NB07 files:',len(nb7_fl))
 
    print('---------------')


    for fltr in suit_filters:
        print(f'Starting with {fltr}')
        fltr_fl = glob.glob(suit_raw_files + '*3'+f'{fltr}.fits')
        fltr_fl=sorted(fltr_fl, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
        if len(fltr_fl)==0:
            print(f'No {fltr} filter files')
            continue
        co_align_maps(suit_raw_files,jpg_fold,ref_fd_img_pth,fltr,fltr_fl,tx1,tx2,ty1,ty2,save_pngs=None,save_fits='yes',draw_contours='yes')

#-------------------------------------------


    
