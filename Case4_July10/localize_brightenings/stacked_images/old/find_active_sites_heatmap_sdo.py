

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
import matplotlib.colors as colors
from astropy.visualization import ImageNormalize, AsinhStretch
import sys
sys.path.append('/home/adithya/Adithya_repos/pil_detection/utils')
sys.path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()
# Import local libraries
from region_detection import pos_neg_detection
from pil_detection import detection
from video_loading import video_loader

sns_cl = sns.color_palette("coolwarm",as_cmap=True)

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
    heat=heat+diff_img*mask.astype(int)
    return heat

#---------------------------------------------------------------------------------------------


peak_time=datetime.datetime(2024, 7,10 ,6,0, 0)
start_time=datetime.datetime(2024, 7,10 ,5,44, 0)
suit_aligned_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/aligned_crop_fits/'
pk_img="/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/aligned_crop_fits/SUT_T24_0956_000465_Lev1.0_2024-07-10T05.56.20.250_0973NB04.fits"
base_image='SDO' #suit
if base_image=='SDO':
    aia_fl='/media/adithya/Adi_disk4/SUIT_flare_work/case4_jul10/data/HMI/HMI_cutouts/hmi.m_45s.20240710_040000_TAI.2.magnetogram.fits' #'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/aia/aia.lev1_euv_12s.2024-07-10T040010Z.171.image_lev1 (2).fits' #
    hmi_imgs   = '/media/adithya/Adi_disk4/SUIT_flare_work/case4_jul10/data/HMI/HMI_cutouts/test_set/'
fltr='NB04'
thresh_sig=5


#---------------------------------------------------------------------------------------------



fol_nm=os.getcwd() #Custom folder to save contour images
jpg_fold=fol_nm+'/'+'Contour_imgs'
suit_fls = glob.glob(suit_aligned_files + '*3'+f'{fltr}.fits')
dt = pos_neg_detection()
pil_dt = detection(hmi_imgs)
#pil_dt.check_outpath(output_path)
STRENGTH_FILTER = 100
BUFFER_SIZE = 4
PRESERVED_FLUX_RATIO = 0.95
MIN_MPIL_SIZE = 14
pil_orig, label_orig = pil_dt.PIL_detect(pos_gauss = STRENGTH_FILTER, neg_gauss= -STRENGTH_FILTER, size_kernel = BUFFER_SIZE)
hmi_map_=Map(aia_fl)
data_num = pil_dt.check_header(hmi_map_)
masked_pil = dt.mask_pil(label_orig[data_num])
strength_label = pil_dt.filter_by_strength(threshold = PRESERVED_FLUX_RATIO)
masked_filter_RoPIs = dt.mask_pil(strength_label[data_num])
thin_df, thin_binary = pil_dt.thin_strength_label(strength_label)
masked_thinned_MPIL = dt.mask_pil(thin_binary[data_num])
pil_msk=Map(masked_thinned_MPIL,hmi_map_.meta)
nb4_fls=[]

hmi_map_dt=hmi_map_.date
if "CROTA2" in hmi_map_.meta:
    map_rot_angl=int(hmi_map_.meta.get('CROTA2'))
    if map_rot_angl>5:
        hmi_map=hmi_map_.rotate(angle=-map_rot_angl*u.deg)
        pil_msk=pil_msk.rotate(angle=-map_rot_angl*u.deg)
    #hmi_dt=hmi_map.date
    
else:
    aia_map_=hmi_map_


for f in suit_fls:
    timestamp = datetime.datetime.strptime(os.path.basename(f).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f")
    if (timestamp <= peak_time) :#& (timestamp>=start_time):
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
fig=plt.figure(figsize=(16,8))
ax = fig.add_subplot(111, projection=base_map)
norm = ImageNormalize( base_map.data/1e3,stretch=AsinhStretch(a=0.2))
#sdo_norm=
if base_image=='SDO':
    mp=aia_map.reproject_to(peak_map.wcs)
    #im2=mp.plot(axes=ax)
    im2=aia_map.plot(axes=ax)
    #im2=ax.imshow(aia_map.data/1e3,cmap='sdoaia171',alpha=1,vmin=-2,vmax=2,transform=ax.get_transform(aia_map.wcs))
    #im2=ax.imshow(aia_map.data/1e3,cmap='sdoaia171',alpha=1,transform=ax.get_transform(aia_map.wcs))
    # im2=base_map.plot(axes=ax)
    #im2=ax.imshow(base_map.data/1e3,cmap='suit_nb04',norm=norm)
else:
    base_map.plot(axes=ax,cmap='suit_nb04')
peak_map.draw_contours(axes=ax,levels=[np.max(peak_map.data)*.6],colors='magenta',linewidth=.4)
sns_cl = sns.color_palette("bwr",as_cmap=True)
mask = np.isnan(pil_msk.data) | (pil_msk.data < 0.08)
masked_pil = np.ma.array(pil_msk.data, mask=mask)
masked_pil.data[~masked_pil.mask] = 1.0
mask_pil_map=Map(masked_pil,pil_msk.meta)

#aia_map.plot(axes=ax)
#   print(ax.get_xlim,ax.get_ylim)

from matplotlib.colors import ListedColormap

pil_cmap = ListedColormap(['yellow'])


for i in range(len(nb4_fls)): #nb4_fls
    suit_map=suit_sq[i]
    #pos.append(draw_suit_contours_on_sdo(suit_map,base_map,thresh_sig,'red'))
    heat=draw_suit_contours_on_sdo(suit_map,base_map,thresh_sig,heat,'red')
heat_masked = ma.masked_where(heat == 0, heat)
im=ax.imshow(heat_masked/1000, origin='lower', cmap='rainbow', alpha=0.4)
#ax.imshow(masked_thinned_MPIL,'spring',transform=ax.get_transform(hmi_map_.wcs)) #cmap=pil_cmap,Wistia'
mask_pil_map.plot(axes=ax,cmap=pil_cmap)
#print(ax.get_xlim)

plt.colorbar(im,label=r'Excess Intensity ($\times 10^3$)')
#plt.colorbar(im2,label=r'Mg II h Intensity ($\times 10^3$ DN)')
plt.colorbar(im2,label=r'LOS Magnetic field ($\times 10^3$ Gauss)')
#plt.colorbar(im2,label=r'Intensity ')
ax.set_xlabel("Solar X (arcsec)")
ax.set_ylabel("Solar Y (arcsec)")
ax.set_xlim(0,600)
ax.set_ylim(0,600) #hmi
# ax.set_xlim(20,550)
# ax.set_ylim(150,500)
plt.savefig(f'HMI_{fltr}_stacked_contours.png',dpi=400)
plt.show()

