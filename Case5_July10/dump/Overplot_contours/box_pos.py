

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

suit_aligned_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop/'
fol_nm=os.getcwd() #Custom folder to save contour images
jpg_fold=fol_nm+'/'+'Contour_imgs'
fltr='NB04'
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



def draw_suit_contours_on_sdo(suitMap,base_map,thresh_sig,fltrl):
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
    fig=plt.figure(figsize=(16,16))
    ax = fig.add_subplot(111, projection=suitMap)
    suitMap.plot(axes=ax)
    coords = []
    def onselect(eclick, erelease):
        coords.append((eclick.xdata, eclick.ydata, erelease.xdata, erelease.ydata))

    toggle_selector = RectangleSelector(ax, onselect, useblit=True,
                      button=[1], minspanx=5, minspany=5, spancoords='pixels',
                      interactive=True)
    norm_mg_Map.draw_contours(axes=ax, levels=v_Thresh,linewidths=2,colors="#106D10",alpha=1)
    plt.show()

    if not coords:
        raise RuntimeError("ROI selection cancelled or failed.")

    x1, y1, x2, y2 = coords[0]
    
    bottom_left = (min(x1, x2), min(y1, y2)) * u.pix
    top_right = (max(x1, x2), max(y1, y2)) * u.pix

    bl = suitMap.pixel_to_world(bottom_left[0], bottom_left[1] )
    tr = suitMap.pixel_to_world(top_right[0] , top_right[1] )
    print([float(bl.Tx.value) ,float(bl.Ty.value),float(tr.Tx.value) ,float(tr.Ty.value)])
    print('---------')
    return [float(bl.Tx.value) ,float(bl.Ty.value),float(tr.Tx.value) ,float(tr.Ty.value)]

#---------------------------------------------------------------------------------------------

suit_sq=Map(nb4_fls,sequence=True)
data_stack = np.stack([(suit_sq[i].data*1000/suit_sq[i].meta.get('CMD_EXPT'))for i in range(5)])
base_img=np.median(data_stack, axis=0)
base_map=Map(base_img,suit_sq[0].meta)
thresh_sig=5
    
#for i in tqdm (range(len(nb4_fls))):
tx_array=[]
ty_array=[]
pos=[]
sel_frams=[2,2,31,49,51,68]
for i in sel_frams:
    suit_map=suit_sq[i]
    pos.append(draw_suit_contours_on_sdo(suit_map,base_map,thresh_sig,fltr))
print(pos)

