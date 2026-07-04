from multiprocessing import Pool, cpu_count
from functools import partial
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
from tqdm import tqdm


import warnings
import logging
warnings.simplefilter('ignore')
log_ = logging.getLogger('sunpy')
log_.setLevel('WARNING')
logging.getLogger('reproject').setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.ERROR)

def process_one_frame(i, suit_sq, hmi_fls, hmi_dt_array, ref_suit_map, ref_hmi_wcs):
    # --- SUIT map ---
    suit_map = suit_sq[i]
    base_time = Time(parse_time(suit_map.date))

    # --- closest HMI ---
    idx = np.argmin(np.abs(hmi_dt_array - base_time))
    hmi_map_ = Map(hmi_fls[idx])

    # --- CROTA2 handling ---
    if "CROTA2" in hmi_map_.meta:
        map_rot_angl = int(hmi_map_.meta.get("CROTA2"))
        if map_rot_angl > 5:
            hmi_map = hmi_map_.rotate(angle=-180 * u.deg)
        else:
            hmi_map = hmi_map_
    else:
        hmi_map = hmi_map_
    

    # --- rebin to SUIT scale ---
    hmi_rebinned = get_suit_scale_rebined_map(hmi_map, ref_suit_map)

    # --- differential rotation ---
    with propagate_with_solar_surface():
        hmi_map_drot = hmi_rebinned.reproject_to(
            ref_hmi_wcs,
            dask_method="none"
        )

    # --- plotting / saving ---
    hmi_date=hmi_rebinned.date.datetime
    draw_suit_contours_on_sdo(suit_map, hmi_map_drot,'NB04','HMI',hmi_date)

    # return minimal data only
    return hmi_rebinned.date.datetime


def get_suit_scale_rebined_map(ref_fd_1600,ref_mg_rot):
    scale=ref_fd_1600.scale[0].value/ref_mg_rot.scale[0].value
    fd_new_dem=[ref_fd_1600.data.shape[1]*scale,ref_fd_1600.data.shape[0]*scale]*u.pixel
    ref_aia_resmp=ref_fd_1600.resample(fd_new_dem)
    blo = SkyCoord(ref_mg_rot.bottom_left_coord.Tx, ref_mg_rot.bottom_left_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    tro = SkyCoord(ref_mg_rot.top_right_coord.Tx, ref_mg_rot.top_right_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    return ref_aia_resmp.submap(blo,top_right=tro)  


def draw_suit_contours_on_sdo(suitMap,aia_map,fltr,base_channel,hmi_date):
    pathlib.Path(jpg_fold+f'/{fltr}_{base_channel}').mkdir(parents=True, exist_ok=True)
    points=np.loadtxt('boxes.csv',delimiter=',')
    pts=np.array(points,dtype=float)
    sns_cl3=sns.color_palette('bright')
    
    reg_lc=[]
    for i in range(len(pts)):
        fl_nm=jpg_fold+f'/{fltr}_{base_channel}'+f'/box_{i}_'+os.path.basename(suitMap.meta['F_NAME'])[:-4]+'jpg'
        fig=plt.figure(figsize=(20,16))
        
        ax = fig.add_subplot(111, projection=aia_map)
        ax.set_title(str(hmi_date))
        point1 = SkyCoord(pts[i][0]*u.arcsec, pts[i][1]*u.arcsec, frame=aia_map.coordinate_frame)
        point2 = SkyCoord(pts[i][2]*u.arcsec, pts[i][3]*u.arcsec, frame=aia_map.coordinate_frame)
        coords = SkyCoord(Tx=(point1.Tx.value, point2.Tx.value) * u.arcsec,Ty=(point1.Ty.value, point2.Ty.value) * u.arcsec,frame=aia_map.coordinate_frame,)
        #aia_map.draw_quadrangle(coords,axes=ax,edgecolor=sns_cl3[i],linestyle="-",linewidth=2,label=f'Box {i}')
        #ax.legend(fontsize=15)
        suit_box=aia_map.submap(coords)
        im=suit_box.plot(axes=ax)
        plt.colorbar(im)
        plt.savefig(fl_nm,dpi=200)
        plt.close()  
    return  





if __name__=='__main__':
#-----------------Intial paths and params--------------------------
    suit_aligned_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/aligned_crop/'

    hmi_imgs='/media/adithya/Adi_disk4/SUIT_flare_work/case6_oct09/data/HMI/HMI_cutouts/'
    fol_nm=os.getcwd() #Custom folder to save contour images
    jpg_fold=fol_nm+'/'+'hmi_cutout_imgs'
    fltr='NB04'
    thresh_sig=5
    base_channel='HMI'

    nb4_fls = glob.glob(suit_aligned_files + '*3'+f'{fltr}.fits')
    hmi_fls = glob.glob(hmi_imgs + '*.fits') 
    print('No of SUIT files:',len(nb4_fls))
    print('No of HMI images: ',len(hmi_fls))
    
    

    hmi_dt=[]
    lc_mtx=[]
    dt_array=[]

    for j in range(len(hmi_fls)):
        hmi_map=Map(hmi_fls[j])
        hmi_dt.append(hmi_map.date.datetime)
    hmi_dt_array=Time(parse_time(hmi_dt))
    suit_sq=Map(nb4_fls,sequence=True)
    ref_suit_map= Map(suit_sq[0])
    base_time=Time(parse_time(ref_suit_map.date))
    idx=np.argmin(np.abs(hmi_dt_array - base_time))
    ref_hmi_map=Map(hmi_fls[idx])
   
    if "CROTA2" in ref_hmi_map.meta:
        map_rot_angl = int(ref_hmi_map.meta.get("CROTA2"))
        if map_rot_angl > 5:
            
            ref_hmi_map = ref_hmi_map.rotate(angle=-180 * u.deg)
        else:
            ref_hmi_map = ref_hmi_map
    else:
        ref_hmi_map = ref_hmi_map
   
    nproc = min(cpu_count() - 1, 8)  # safe limit
    print('no of processors: ',nproc)
    worker = partial(
        process_one_frame,
        suit_sq=suit_sq,
        hmi_fls=hmi_fls,
        hmi_dt_array=hmi_dt_array,
        ref_suit_map=ref_suit_map,
        ref_hmi_wcs=ref_hmi_map.wcs)

    with Pool(processes=nproc) as pool:
        dt_array = list(
            tqdm(
                pool.imap(worker, range(len(nb4_fls))),
                total=len(nb4_fls)
            )
        )
