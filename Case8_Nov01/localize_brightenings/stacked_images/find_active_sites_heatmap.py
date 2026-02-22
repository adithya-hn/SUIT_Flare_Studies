

'''
Created on Dec 1st 2025
Author=adithya-hn
Purpose: get the coordonates of the 
'''


import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from sunpy.map import MapSequence
from sunpy.net import Fido
import glob
import datetime
from datetime import timedelta
import timeit
import pathlib
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
from PIL import Image
import pandas as pd
from sunpy.time import parse_time
from astropy.time import Time
from tqdm import tqdm
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from matplotlib.widgets import RectangleSelector
from skimage.morphology import disk, closing, opening,remove_small_objects
from skimage import filters, measure
import matplotlib.patches as patches
import seaborn as sns
from skimage.morphology import disk, closing,opening,dilation
import numpy.ma as ma
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
#set_pub_style()


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

    if not coords:
        raise RuntimeError("ROI selection cancelled or failed.")

    x1, y1, x2, y2 = coords[0]
    
    bottom_left = (min(x1, x2), min(y1, y2)) * u.pix
    top_right = (max(x1, x2), max(y1, y2)) * u.pix

    submap = sunpy_map.submap(bottom_left=bottom_left, top_right=top_right)
    return submap

def rebin_suit_map(aia_map,ref_mg_rot):
    scale=ref_mg_rot.scale[0].value/aia_map.scale[0].value
    fd_new_dem=[aia_map.data.shape[1]*scale,aia_map.data.shape[0]*scale]*u.pixel
    ref_aia_resmp=aia_map.resample(fd_new_dem)
    blo = SkyCoord(ref_mg_rot.bottom_left_coord.Tx, ref_mg_rot.bottom_left_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    tro = SkyCoord(ref_mg_rot.top_right_coord.Tx, ref_mg_rot.top_right_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    return ref_aia_resmp.submap(blo,top_right=tro) 

def draw_suit_contours_on_sdo(suitMap,base_map,thresh_sig,heat,clr):
    
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
    heat=heat+mask.astype(int)
    return heat

#---------------------------------------------------------------------------------------------


peak_time=datetime.datetime(2024, 11,1 ,14,31, 0)
start_time=datetime.datetime(2024, 11,1 ,14,18, 0)
suit_aligned_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/data/aligned_crop/'
pk_img="/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/data/aligned_crop/SUT_T24_1592_000636_Lev1.0_2024-11-01T14.28.46.392_0973NB04.fits"
aia_fl='/media/adithya/Adi_disk4/SUIT_flare_work/case8_nov01/data/aia/cut_outs/171_cutouts/aia.lev1_euv_12s.2024-11-01T113110Z.171.image_lev1.fits'
fltr='NB04'
thresh_sig=5

#---------------------------------------------------------------------------------------------

fol_nm=os.getcwd() #Custom folder to save contour images
jpg_fold=fol_nm+'/'+'Contour_imgs'
suit_fls = glob.glob(suit_aligned_files + '*3'+f'{fltr}.fits')

nb4_fls=[]
for f in suit_fls:
    timestamp = datetime.datetime.strptime(os.path.basename(f).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f")
    if (timestamp <= start_time) :#& (timestamp>=start_time):
        nb4_fls.append(f)

print('No of SUIT files:',len(nb4_fls))

suit_sq=Map(nb4_fls,sequence=True)
data_stack = np.stack([(suit_sq[i].data*1000/suit_sq[i].meta.get('CMD_EXPT'))for i in range(5)])
base_img=np.median(data_stack, axis=0)
base_map=Map(base_img,suit_sq[0].meta)

heat = np.zeros(base_map.data.shape)

tx_array=[]
ty_array=[]
pos=[]
base_drot=base_map.rotate(angle=(int(base_map.meta.get('p_angle')))*u.deg)

peak_map_=Map(pk_img)
peak_map=Map(peak_map_.data*1000/peak_map_.meta.get('CMD_EXPT'),peak_map_.meta)
print(np.max(peak_map.data))

aia_map_=Map(aia_fl)
aia_map =rebin_suit_map(aia_map_,base_drot)
fig=plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection=base_map)

peak_map.draw_contours(axes=ax,levels=[np.max(peak_map.data)*.6],colors='k')
sns_cl = sns.color_palette("coolwarm",as_cmap=True)
#aia_map.plot(axes=ax)
base_map.plot(axes=ax)

for i in range(len(nb4_fls)): #nb4_fls
    suit_map=suit_sq[i]
    #pos.append(draw_suit_contours_on_sdo(suit_map,base_map,thresh_sig,'red'))
    heat=draw_suit_contours_on_sdo(suit_map,base_map,thresh_sig,heat,'red')
heat_masked = ma.masked_where(heat == 0, heat)
im=ax.imshow(heat_masked, origin='lower', cmap='rainbow', alpha=0.6,label='stacked contours')

plt.colorbar(im,label='No. of stacked contours')
plt.savefig(f'{fltr}_stacked_contours.png',dpi=300)#aia_171_
plt.show()

