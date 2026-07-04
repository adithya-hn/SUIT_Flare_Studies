
'''
Created on 28 Dec 2025
Author: adithya-hn
Purpose: co align SUIT Mg II h  image to image align with AIA1600


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
from matplotlib.widgets import RectangleSelector
from sunpy.coordinates import RotatedSunFrame
import datetime
from astropy.time import TimeDelta
from astropy.io import fits

import warnings
import logging

warnings.simplefilter('ignore')
log_ = logging.getLogger('sunpy')
log_.setLevel('WARNING')
logging.getLogger('reproject').setLevel(logging.WARNING)

#------------------------------------------------------------------------------#

def tracked_sunspot(first_map,tx,ty,current_map):
    coord = SkyCoord(tx * u.arcsec, ty * u.arcsec    , frame=first_map.coordinate_frame)
    t1 = Time(first_map.date)   # observation time of first map
    t2 = Time(current_map.date)   # observation time of second map

    durations = (t2 - t1).to('min')
    diffrot_point = SkyCoord(RotatedSunFrame(base=coord, duration=durations))
    diff_rotated=diffrot_point.transform_to(current_map.coordinate_frame)
    blo = SkyCoord((diff_rotated.Tx.value-50)*u.arcsec, (diff_rotated.Ty.value-50)*u.arcsec, frame=current_map.coordinate_frame)
    tro = SkyCoord((diff_rotated.Tx.value+50)*u.arcsec, (diff_rotated.Ty.value+50)*u.arcsec, frame=current_map.coordinate_frame)
    return current_map.submap(blo,top_right=tro)  



def get_rebinned_crop_map(ref_fd_1600,ref_mg_rot):
    scale=ref_fd_1600.scale[0].value/ref_mg_rot.scale[0].value
    fd_new_dem=[ref_fd_1600.data.shape[1]*scale,ref_fd_1600.data.shape[0]*scale]*u.pixel
    ref_aia_resmp=ref_fd_1600.resample(fd_new_dem)
    blo = SkyCoord(ref_mg_rot.bottom_left_coord.Tx, ref_mg_rot.bottom_left_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    tro = SkyCoord(ref_mg_rot.top_right_coord.Tx, ref_mg_rot.top_right_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    return ref_aia_resmp.submap(blo,top_right=tro)  

def get_mg_threshold(ref_mg_rot):
    valid = ref_mg_rot.data[(ref_mg_rot.data > 100)]
    valid_int = valid.astype(int)  
    mode_val = stats.mode(valid_int, keepdims=True).mode[0]/ref_mg_rot.meta.get('CMD_EXPT')*1000  #normalised mode value    
    th_lvs=np.array([mode_val*0.7,mode_val*1.5,mode_val*3] )
    return mode_val,th_lvs

def draw_contours_and_save(BaseMap,suit_aligned_map,suit_thresh_levs,save_path):
    pathlib.Path(save_path).mkdir(parents=True, exist_ok=True) #make sure dir exists
    fl_nm=save_path+'/aia_1600_'+str(suit_aligned_map.meta["f_name"])[:-5]+'.jpg'
    fig=plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection=BaseMap)
    BaseMap.plot(axes=ax,cmap='gray',title=str(BaseMap.date))
    norm_data=suit_aligned_map.data*1000/int(suit_aligned_map.meta.get('CMD_EXPT'))
    norm_mg_Map=Map(gaussian_filter(norm_data, sigma=1),suit_aligned_map.fits_header)
    norm_mg_Map.draw_contours(axes=ax, levels=suit_thresh_levs,lws=0.5,colors=['red','pink','green'])
    plot_str='AIA 1600: '+str(BaseMap.date) +'\n'+ 'Mg II h: '+str(suit_aligned_map.date) 
    ax.text(50,50, plot_str, color='white', fontsize=10)
    plt.draw()
    #plt.colorbar()
    plt.savefig(fl_nm,dpi=300)
    plt.close()  

#--------------------------------------------------------------------------------------------------------

suit_raw_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/raw/'
aia_roi_files='/media/adithya/Adi_disk4/SUIT_flare_work/case6_oct09/data/aia/cut_outs/1600_cutouts/'

fltr='NB04'

fltr_fl = glob.glob(suit_raw_files + '*3'+f'{fltr}.fits')
fltr_fl=sorted(fltr_fl, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

aia1600s=glob.glob(aia_roi_files + '*.fits')#
aia1600s.sort()
#times_str = np.genfromtxt('aia_times.csv', delimiter=',', dtype=str)

dt_list = []
for t in aia1600s:
    header=fits.getheader(t)
    dt_list.append(datetime.datetime.fromisoformat(header['DATE-OBS']))
times=Time(parse_time(dt_list))
datetime_array = times#.to_datetime()
#print(datetime_array)

print('No. of SUIT images: ',len(fltr_fl))
print('No. of AIA 1600 images: ',len(aia1600s))

suit_map=Map(fltr_fl[0])
aia_map_=Map(aia1600s[0])

for i in tqdm(range(len(fltr_fl))):
    suit_map=Map(fltr_fl[i])
    suit_map_rot=suit_map.rotate(angle=suit_map.meta["CROTA2"] * u.deg,missing=0)
    base_time=Time(parse_time(suit_map_rot.date))
    idx=np.argmin(np.abs(datetime_array - base_time))
    dt = np.array([(t - base_time).to_value('s')  for t in datetime_array])
    # keep only times ≤ base_time
    dt[dt > 0] = -np.inf
    idx = np.argmax(dt)   # largest negative (or zero) offset
    #print(idx,i)
    aia_map=Map(aia1600s[idx])
    aia_rebin_map=get_rebinned_crop_map(aia_map,suit_map_rot)
    sq=MapSequence([aia_rebin_map,suit_map_rot])

    template= tracked_sunspot(aia_map_,10,0,aia_rebin_map)

    mg_aln_maps=mc_coalign(sq,layer_index=0,template=template,clip=False)#,func=np.asinh)

    thresh_lvs=get_mg_threshold(suit_map_rot)
    save_path='../data/aln_1600_conts'
    fits_path='../data/aligned_fits/'
    png_path='../data/aligned_png/'
    draw_contours_and_save(aia_rebin_map,mg_aln_maps[1],thresh_lvs[1],save_path)

    alnMap=Map(mg_aln_maps[1].data,mg_aln_maps[1].meta)
    alnMap.meta['CRPIX1']=mg_aln_maps[0].meta['CRPIX1']
    alnMap.meta['CRPIX2']=mg_aln_maps[0].meta['CRPIX2']
    pathlib.Path(fits_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(png_path).mkdir(parents=True, exist_ok=True)
    try:

        fname=png_path+str(mg_aln_maps[1].meta.get('F_NAME'))[:-5]+'.png'
        alnMap_=alnMap.rotate(angle=-alnMap.meta['p_angle']*u.deg,missing=0)
        alnMap_.save(fits_path+mg_aln_maps[1].meta.get('F_NAME'),overwrite=True)
    except:
        print('could not save file')
    
    try:
        alnMap.plot()
        plt.savefig(fname)
        plt.close()
    except:
        pass
   
    


# co align SUIT AIA
# draw contour check align

