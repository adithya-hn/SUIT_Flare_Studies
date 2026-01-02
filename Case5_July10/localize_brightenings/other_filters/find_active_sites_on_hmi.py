

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
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()



suit_aligned_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop/'
hmi_map_raw=Map('/media/adithya/Adi_disk4/SUIT_flare_work/case5_jul10/data/HMI/HMI_cutouts/hmi.m_45s.20240710_133730_TAI.2.magnetogram.fits')
fol_nm=os.getcwd() #Custom folder to save contour images
jpg_fold=fol_nm+'/'+'Contour_imgs'
fltr='NB02'
nb4_fls = glob.glob(suit_aligned_files + '*3'+f'{fltr}.fits')
print('No of SUIT files:',len(nb4_fls))



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

def get_suit_scale_rebined_map(ref_fd_1600,ref_mg_rot):
    scale=ref_fd_1600.scale[0].value/ref_mg_rot.scale[0].value
    fd_new_dem=[ref_fd_1600.data.shape[1]*scale,ref_fd_1600.data.shape[0]*scale]*u.pixel
    ref_aia_resmp=ref_fd_1600.resample(fd_new_dem)
    blo = SkyCoord(ref_mg_rot.bottom_left_coord.Tx, ref_mg_rot.bottom_left_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    tro = SkyCoord(ref_mg_rot.top_right_coord.Tx, ref_mg_rot.top_right_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    return ref_aia_resmp.submap(blo,top_right=tro)  


def draw_suit_contours_on_sdo(suitMap,base_map,thresh_sig,clr):
    
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
    if clr=='green':
        norm_mg_Map.draw_contours(axes=ax, levels=v_Thresh,linewidths=2,colors="#106D10",alpha=0.1)
    if clr=='red':
        norm_mg_Map.draw_contours(axes=ax, levels=v_Thresh,linewidths=2,colors="#ED2913",alpha=0.5)

    return fig

#---------------------------------------------------------------------------------------------

suit_sq=Map(nb4_fls,sequence=True)
data_stack = np.stack([(suit_sq[i].data*1000/suit_sq[i].meta.get('CMD_EXPT'))for i in range(5)])
base_img=np.median(data_stack, axis=0)
base_map=Map(base_img,suit_sq[0].meta)
thresh_sig=3
if hmi_map_raw.meta['CROTA2'] >5:
    hmi_map_raw=hmi_map_raw.rotate(angle=180*u.deg)


hmi_map=get_suit_scale_rebined_map(hmi_map_raw,base_map)
#for i in tqdm (range(len(nb4_fls))):
tx_array=[]
ty_array=[]
pos=[]
sel_frams=[2,2,31,49,51,68]

fig=plt.figure(figsize=(16,16))
ax = fig.add_subplot(111, projection=hmi_map)
#base_map.plot(axes=ax)
sns_cl = sns.color_palette("coolwarm",as_cmap=True)
sns_cl3=sns.color_palette('bright')
im=ax.imshow(hmi_map.data,cmap=sns_cl,alpha=1,vmin=-2000,vmax=2000)
points=np.loadtxt('boxes.csv',delimiter=',')
pts=np.array(points,dtype=float)
plt.rcParams["font.size"]=50
plt.rcParams["axes.labelsize"]=50
plt.rcParams["xtick.labelsize"]=50
plt.rcParams["ytick.labelsize"]=50
plt.rcParams["legend.fontsize"]=20
plt.rcParams["figure.titlesize"]=50
plt.rcParams["axes.titlesize"]=50

for i in range(pts.shape[0]):
    point1 = SkyCoord(pts[i][0]*u.arcsec, pts[i][1]*u.arcsec, frame=hmi_map.coordinate_frame)
    point2 = SkyCoord(pts[i][2]*u.arcsec, pts[i][3]*u.arcsec, frame=hmi_map.coordinate_frame)
    coords = SkyCoord(Tx=(point1.Tx.value, point2.Tx.value) * u.arcsec,Ty=(point1.Ty.value, point2.Ty.value) * u.arcsec,frame=hmi_map.coordinate_frame,)
    hmi_map.draw_quadrangle(coords,axes=ax,edgecolor=sns_cl3[i],linestyle="-",linewidth=2,label=f'Box {i}')



mask = np.abs(hmi_map.data)>500 
kernel=np.array([[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]]) #disk(3) #
#print(kernel)
mask_=dilation(mask,kernel)
mask_=dilation(mask_,kernel)
mask_=dilation(mask_,kernel)
mask_=dilation(mask_,kernel)

inv=np.invert(mask_)
data = hmi_map.data.astype(float)
data[inv] = np.nan                 # assign NaN where mask is True
masked_map = Map(data, hmi_map.meta)
masked_map.draw_contours(axes=ax,levels=0*u.G,colors='k',linewidths=1.5,alpha=0.6)


for i in range(len(suit_sq)):
    suit_map=suit_sq[i]
    pos.append(draw_suit_contours_on_sdo(suit_map,base_map,thresh_sig,'green'))

#draw_suit_contours_on_sdo(pre_flare_map,base_map,thresh_sig,'red')

ax.set_xlim(0,520)
ax.set_ylim(80,420)
# ax.set_xlim(350,900)
# ax.set_ylim(300,900)
plt.legend()
plt.savefig('HMI_plots.pdf',dpi=300)
plt.show()
