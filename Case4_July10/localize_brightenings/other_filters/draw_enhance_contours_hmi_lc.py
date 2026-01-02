
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
import matplotlib.patches as patches


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


def draw_suit_contours_on_sdo(suitMap,base_map,thresh_sig,fltr):
    
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
    # fig=plt.figure(figsize=(16,16))
    # plt.rcParams["font.size"]=50
    # plt.rcParams["axes.labelsize"]=50
    # plt.rcParams["xtick.labelsize"]=50
    # plt.rcParams["ytick.labelsize"]=50
    # plt.rcParams["legend.fontsize"]=50
    # plt.rcParams["figure.titlesize"]=50
    # plt.rcParams["axes.titlesize"]=50
    fig=plt.figure()
    points=np.loadtxt(f'{fltr}_boxes.csv',delimiter=',')
    pts=np.array(points,dtype=float)
 
    sns_cl3=sns.color_palette('bright')
    ax = fig.add_subplot(111, projection=base_map)
    
    bl = SkyCoord(pts[0]*u.arcsec, pts[1]*u.arcsec, frame=suitMap.coordinate_frame)
    tr = SkyCoord(pts[2]*u.arcsec, pts[3]*u.arcsec, frame=suitMap.coordinate_frame)

    # Convert to pixel coordinates
    # Convert to pixel coords
    bl_pix = suitMap.world_to_pixel(bl)
    tr_pix = suitMap.world_to_pixel(tr)
    blp=(bl_pix[0].value,bl_pix[1].value)*u.pix
    trp=(tr_pix[0].value,tr_pix[1].value)*u.pix

    width  = tr_pix[0].value - bl_pix[0].value
    height = tr_pix[1].value - bl_pix[1].value

    
    suit_box = suitMap.submap(bottom_left=blp, top_right=trp)
    rect = patches.Rectangle((bl_pix[0].value,bl_pix[1].value), width, height,
                         fill=False, linewidth=1,
                         transform=ax.get_transform('pixel'))
 
    #coords = SkyCoord(Tx=(point1.Tx.value, point2.Tx.value) * u.arcsec,Ty=(point1.Ty.value, point2.Ty.value) * u.arcsec,frame=suitMap.coordinate_frame,)
    #suitMap.draw_quadrangle(coords,axes=ax,edgecolor=sns_cl3[0],linestyle="-",linewidth=1,label=f'Box {i}')
    ax.legend(fontsize=15)
    ax.add_patch(rect)
    #suit_box=suitMap.submap(coords)
    reg_lc=(np.sum(abs(suit_box.data)))
    suitMap.plot(axes=ax)
    norm_mg_Map.draw_contours(axes=ax, levels=v_Thresh,linewidths=.8,colors="#106D10",alpha=1)
    #ax.text(0.5, 0.95,str(aia_dt)[:-4],transform=ax.transAxes,ha="center", va="center",fontsize=42)
    plt.savefig(fl_nm,dpi=300)
    plt.close()  

    return  reg_lc




if __name__=='__main__':
#-----------------Intial paths and params--------------------------
    suit_aligned_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/aligned_crop/'
    
    fol_nm=os.getcwd() #Custom folder to save contour images
    jpg_fold=fol_nm+'/'+'Contour_imgs'
    fltr='NB02'
    thresh_sig=3
    base_channel='HMI'


    save_aligned_fits='yes'
    save_fits='no'
    save_pngs='yes'     #aligned pngs
    draw_contours='yes'
    
    nb4_fls = glob.glob(suit_aligned_files + '*3'+f'{fltr}.fits')
    print('No of SUIT files:',len(nb4_fls))
 
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


    for i in tqdm (range(len(nb4_fls))): #
        suit_map=suit_sq[i]
        base_time=Time(parse_time(suit_map.date))
        lc_mtx.append(draw_suit_contours_on_sdo(suit_map,base_map,thresh_sig,fltr))
        
    
    out = np.column_stack([dt_array, lc_mtx])
    np.savetxt(f'{fltr}_lc.csv',out,delimiter=',',header='time,reg1,reg2,reg3,reg4,reg5,reg6',fmt='%s',comments='')
