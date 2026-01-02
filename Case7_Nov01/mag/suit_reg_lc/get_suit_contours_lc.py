
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
from matplotlib.patches import Rectangle
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

def draw_suit_contours_on_sdo(suitMap):
    
    pathlib.Path(jpg_fold+f'/{fltr}_{base_channel}').mkdir(parents=True, exist_ok=True)


    SuitMap=Map(suitMap.data*1000/suitMap.meta.get('CMD_EXPT'),suitMap.meta)
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
    ax = fig.add_subplot(111, projection=SuitMap)
    reg_lc=[]
    dt_lc=[]
    width  = 1 * u.arcsec
    height = 1 * u.arcsec
    SuitMap.plot(axes=ax)
    for i in range(len(pts)):
        point1 = SkyCoord(pts[i][0]*u.arcsec, pts[i][1]*u.arcsec, frame=SuitMap.coordinate_frame)
        point2 = SkyCoord(pts[i][2]*u.arcsec, pts[i][3]*u.arcsec, frame=SuitMap.coordinate_frame)
        x_bl, y_bl = SuitMap.world_to_pixel(point1)
        x_tr, y_tr = SuitMap.world_to_pixel(point2)
        x0 = min(x_bl, x_tr)
        y0 = min(y_bl, y_tr)
        width  = abs(x_tr - x_bl)
        height = abs(y_tr - y_bl)

        rect = Rectangle(
            (x0.value, y0.value),
            width.value,
            height.value,
            edgecolor='red',
            facecolor='none',
            linewidth=2)


        #coords = SkyCoord(Tx=(point1.Tx.value, point2.Tx.value) * u.arcsec,Ty=(point1.Ty.value, point2.Ty.value) * u.arcsec,frame=SuitMap.coordinate_frame,)
        #SuitMap.draw_quadrangle(coords,axes=ax,edgecolor=sns_cl3[i],linestyle="-",linewidth=2,label=f'Box {i}')
        ax.add_patch(rect)
        ax.legend(fontsize=15)
        suit_box=SuitMap.submap(bottom_left=(x0, y0) * u.pix,top_right=(x0 + width, y0 + height) * u.pix)
        reg_lc.append(np.sum(abs(suit_box.data)))
    dt_lc.append(reg_lc)
    plt.savefig(fl_nm,dpi=200)
    plt.close()  

    return  reg_lc




if __name__=='__main__':
#-----------------Intial paths and params--------------------------
    suit_aligned_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/aligned_crop/'
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
    print('No of SUIT files:',len(nb4_fls))
    suit_sq=Map(nb4_fls,sequence=True)
    ref_suit_map=Map(suit_sq[0])

    aia_dt=[]
    hmi_dt=[]
    fig_pickle=[]
    lc_mtx=[]
    dt_array=[]


    for i in tqdm (range(len(nb4_fls))): #
        suit_map=suit_sq[i]
        dt_array.append(suit_map.date.datetime)
        lc_mtx.append(draw_suit_contours_on_sdo(suit_map))
        
    
    out = np.column_stack([dt_array, lc_mtx])
    np.savetxt('suit_lc.csv',out,delimiter=',',header='time,reg1,reg2,reg3,reg4,reg5,reg6,reg7',fmt='%s',comments='')
