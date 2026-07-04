
'''
Created on 4 Nov 2025
Author: adithya-hn
Purpose: Modular version of overplot contours of SUIT Mg II k  images to AIA modified to draw HMIbox light cuvrves for list  of boxes


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
from skimage import filters, measure
from skimage.measure import label, regionprops
from skimage.morphology import disk, closing, opening,remove_small_objects
import seaborn as sns
import pickle
from skimage.morphology import disk, closing,opening,dilation
from matplotlib.widgets import RectangleSelector
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

import warnings
import logging
warnings.simplefilter('ignore')
log_ = logging.getLogger('sunpy')
log_.setLevel('WARNING')
logging.getLogger('reproject').setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.ERROR)

#------------------------------------------------------------------------------#
def select_roi_with_mouse(sunpy_map, cmap=None, norm=None):
    """
    DESCRIPTION: To select RoI template.
    INPUT: Sunpy map.
    RETURNS: Sunpy submap. 
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=sunpy_map)
    ax.set_title("Select ROI (click and drag) then close the window")
    sunpy_map.plot(axes=ax)
    coords = []
    def onselect(eclick, erelease):
        coords.append((eclick.xdata, eclick.ydata, erelease.xdata, erelease.ydata))

    toggle_selector = RectangleSelector(ax, onselect, useblit=True,
                      button=[1], minspanx=5, minspany=5, spancoords='pixels',
                      interactive=True)
    plt.show()


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

def draw_suit_contours_on_sdo(suitMap,aia_map,base_map,thresh_sig,fltr,aia_dt,base_channel):
    
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

    fl_nm=jpg_fold+f'/{fltr}_{base_channel}'+'/'+os.path.basename(suitMap.meta['F_NAME'])[:-4]+'jpg'
    fig=plt.figure(figsize=(16,16))
    plt.rcParams["font.size"]=50
    plt.rcParams["axes.labelsize"]=50
    plt.rcParams["xtick.labelsize"]=50
    plt.rcParams["ytick.labelsize"]=50
    plt.rcParams["legend.fontsize"]=50
    plt.rcParams["figure.titlesize"]=50
    plt.rcParams["axes.titlesize"]=50

    
    points=np.loadtxt('boxes.csv',delimiter=',')
    pts=np.array(points,dtype=float)
    sns_cl = sns.color_palette("coolwarm",as_cmap=True)
    sns_cl2=sns.color_palette('colorblind')
    sns_cl3=sns.color_palette('bright')
    ax = fig.add_subplot(111, projection=aia_map)
    
    reg_lc=[]
    dt_lc=[]
    width  = 1 * u.arcsec
    height = 1 * u.arcsec
    for i in range(len(pts)):
        point1 = SkyCoord(pts[i][0]*u.arcsec, pts[i][1]*u.arcsec, frame=aia_map.coordinate_frame)
        point2 = SkyCoord(pts[i][2]*u.arcsec, pts[i][3]*u.arcsec, frame=aia_map.coordinate_frame)
        coords = SkyCoord(Tx=(point1.Tx.value, point2.Tx.value) * u.arcsec,Ty=(point1.Ty.value, point2.Ty.value) * u.arcsec,frame=aia_map.coordinate_frame,)
        aia_map.draw_quadrangle(coords,axes=ax,edgecolor=sns_cl3[i],linestyle="-",linewidth=2,label=f'Box {i}')
        ax.legend(fontsize=15)
        suit_box=aia_map.submap(coords)
        reg_lc.append(np.sum(abs(suit_box.data)))
    #dt_lc.append(reg_lc)


    im=ax.imshow(aia_map.data,cmap=sns_cl,alpha=1,vmin=-2000,vmax=2000)
    
    if base_channel=='HMI':

        mask = np.abs(aia_map.data)>500 
        kernel=np.array([[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]]) #disk(3) #
        #print(kernel)
        mask_=dilation(mask,kernel)
        mask_=dilation(mask_,kernel)
        mask_=dilation(mask_,kernel)
        mask_=dilation(mask_,kernel)

        inv=np.invert(mask_)
        data = aia_map.data.astype(float)
        data[inv] = np.nan                 # assign NaN where mask is True
        masked_map = Map(data, aia_map.meta)
        masked_map.draw_contours(axes=ax,levels=0,colors='k',linewidths=1.5,alpha=0.6)
        
    norm_mg_Map.draw_contours(axes=ax, levels=v_Thresh,linewidths=2,colors="#106D10",alpha=1)
    ax.set_xlim(400,900)
    ax.set_ylim(200,600)
    ax.text(0.5, 0.95,str(aia_dt)[:-4],transform=ax.transAxes,ha="center", va="center",fontsize=42)
    #suitMap.draw_contours(axes=ax, levels=31000,lws=0.5,colors=['red'],alpha=0.5)
    #plt.colorbar(im, ax=ax)

    plt.savefig(fl_nm,dpi=200)
    plt.close()  

    return  reg_lc




if __name__=='__main__':
#-----------------Intial paths and params--------------------------
    suit_aligned_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop/'
    aia_imgs='/media/adithya/Adi_disk4/SUIT_flare_work/case5_Jul10/data/aia/cut_outs/171_cutouts/'
    hmi_imgs='/media/adithya/Adi_disk4/SUIT_flare_work/case5_jul10/data/HMI/HMI_cutouts/'
    fol_nm=os.getcwd() #Custom folder to save contour images
    jpg_fold=fol_nm+'/'+'Contour_imgs'
    fltr='NB04'
    thresh_sig=4
    base_channel='HMI'


    save_aligned_fits='yes'
    save_fits='no'
    save_pngs='yes'     #aligned pngs
    draw_contours='yes'
    
    nb4_fls = glob.glob(suit_aligned_files + '*3'+f'{fltr}.fits')
    aia_fls = glob.glob(aia_imgs + '*.fits')
    hmi_fls = glob.glob(hmi_imgs + '*.fits')
    

    print('No of SUIT files:',len(nb4_fls))
    print('No of AIA images: ',len(aia_fls))
    print('No of HMI images: ',len(hmi_fls))
    
    suit_sq=Map(nb4_fls,sequence=True)
    ref_suit_map=Map(suit_sq[0])
    
    

    data_stack = np.stack([(suit_sq[i].data*1000/suit_sq[i].meta.get('CMD_EXPT'))for i in range(5)])
    base_img=np.median(data_stack, axis=0)
    base_map=Map(base_img,suit_sq[0].meta)
    

    aia_dt=[]
    hmi_dt=[]
    fig_pickle=[]
    lc_mtx=[]
    dt_array=[]

    for j in range(len(aia_fls)):
        aia_map=Map(aia_fls[j])
        aia_dt.append(aia_map.date.datetime)
    aia_dt_array=Time(parse_time(aia_dt))

    for j in range(len(hmi_fls)):
        hmi_map=Map(hmi_fls[j])
        hmi_dt.append(hmi_map.date.datetime)
    hmi_dt_array=Time(parse_time(hmi_dt))
   

    for i in tqdm (range(len(nb4_fls))): #
        suit_map=suit_sq[i]
        base_time=Time(parse_time(suit_map.date))
        idx1=np.argmin(np.abs(aia_dt_array - base_time))
        aia_map=Map(aia_fls[idx1])
        aia_dt=aia_map.date
        idx=np.argmin(np.abs(hmi_dt_array - base_time))
        hmi_map_=Map(hmi_fls[idx])
        
        if "CROTA2" in hmi_map_.meta:        # case-insensitive in FITS
            #print("CROTA2 exists:", hmi_map_.meta["CROTA2"])
            map_rot_angl=int(hmi_map_.meta.get('CROTA2'))
            if map_rot_angl>5:
                hmi_map=hmi_map_.rotate(angle=-map_rot_angl*u.deg)
                print('Rotation applied')
            else:
                hmi_map=hmi_map_
        else:
            #print("CROTA2 not present")
            hmi_map=hmi_map_
        if i==0:
            ref_aia_map=aia_map#Map(aia_fls[idx1])
            ref_hmi_map=hmi_map#Map(hmi_fls[idx])
        #ref_hmi_map.peek()
        aia_rebinned=get_suit_scale_rebined_map(aia_map,ref_suit_map)
        hmi_rebinned=get_suit_scale_rebined_map(hmi_map,ref_suit_map)
        dt_array.append(hmi_rebinned.date.datetime)
        #print(base_time,hmi_map.date)
        with propagate_with_solar_surface():
            #aia_map_drot=(aia_rebinned.reproject_to(ref_aia_map.wcs,dask_method='none')) #,parallel=True,dask_method='memmap'
            hmi_map_drot=(hmi_rebinned.reproject_to(ref_hmi_map.wcs,dask_method='none')) 

        lc_mtx.append(draw_suit_contours_on_sdo(suit_map,hmi_map_drot,base_map,thresh_sig,fltr,aia_dt,base_channel))
        
    
    out = np.column_stack([dt_array, lc_mtx])
    np.savetxt('hmi_lc.csv',out,delimiter=',',header='time,reg1,reg2,reg3,reg4,reg5,reg6,reg7',fmt='%s',comments='')
