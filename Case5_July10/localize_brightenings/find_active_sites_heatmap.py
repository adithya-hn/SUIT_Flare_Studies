

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
    # contours = measure.find_contours(mask, 0.5)
    
    # # add contour pixels to heatmap
    # for c in contours:
    #     c = np.round(c).astype(int)
    #     heat[c[:,0], c[:,1]] += 1   # increment count
    # #th_diff_img=mask.astype(int)*diff_img
    # #norm_mg_Map=Map(th_diff_img,suitMap.fits_header)
    # # if clr=='green':
    # #     norm_mg_Map.draw_contours(axes=ax, levels=v_Thresh,linewidths=2,colors="#106D10",alpha=0.1)
    # # if clr=='red':
    # #     norm_mg_Map.draw_contours(axes=ax, levels=v_Thresh,linewidths=2,colors="#ED2913",alpha=0.5)

    return heat

#---------------------------------------------------------------------------------------------

suit_aligned_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop/'
peak_time=datetime.datetime(2024, 7,10 ,15,25, 0)
start_time=datetime.datetime(2024, 7,10 ,15,25, 0)
pk_img="/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.35.37.129_0983NB04.fits"
    
fol_nm=os.getcwd() #Custom folder to save contour images
jpg_fold=fol_nm+'/'+'Contour_imgs'
fltr='NB04'
suit_fls = glob.glob(suit_aligned_files + '*3'+f'{fltr}.fits')
#files =sorted(suit_fls, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

nb4_fls=[]
for f in suit_fls:
    # # extract timestamp from your filename
    # ts_str = name.split("_")[1] + name.split("_")[2].split(".")[0]
    timestamp = datetime.datetime.strptime(os.path.basename(f).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f")

    if (timestamp <= peak_time) :#& (timestamp>=start_time):
        nb4_fls.append(f)

print('No of SUIT files:',len(nb4_fls))



suit_sq=Map(nb4_fls,sequence=True)
data_stack = np.stack([(suit_sq[i].data*1000/suit_sq[i].meta.get('CMD_EXPT'))for i in range(5)])
base_img=np.median(data_stack, axis=0)
base_map=Map(base_img,suit_sq[0].meta)
thresh_sig=5
heat = np.zeros(base_map.data.shape)

#for i in tqdm (range(len(nb4_fls))):
tx_array=[]
ty_array=[]
pos=[]
base_drot=base_map#.rotate(angle=(int(base_map.meta.get('p_angle')))*u.deg)

peak_map_=Map(pk_img)
peak_map=Map(peak_map_.data*1000/peak_map_.meta.get('CMD_EXPT'),peak_map_.meta)
print(np.max(peak_map.data))



fig=plt.figure(figsize=(16,16))
ax = fig.add_subplot(111, projection=base_drot)
base_drot.plot(axes=ax)
peak_map.draw_contours(axes=ax,levels=[np.max(peak_map.data)*.6])
sns_cl = sns.color_palette("coolwarm",as_cmap=True)

#np.max(peak_map.data)*.4,np.max(peak_map.data)*.5,
for i in range(len(nb4_fls)): #nb4_fls
    suit_map=suit_sq[i]
    #pos.append(draw_suit_contours_on_sdo(suit_map,base_map,thresh_sig,'red'))
    heat=draw_suit_contours_on_sdo(suit_map,base_map,thresh_sig,heat,'red')
heat_masked = ma.masked_where(heat == 0, heat)
ax.imshow(heat_masked, origin='lower', cmap='rainbow', alpha=0.6)
plt.show()


#peak_map.peek()
#5331



# coords = []
# rects = [] 
# def onselect(eclick, erelease):

#     x1, y1 = eclick.xdata, eclick.ydata
#     x2, y2 = erelease.xdata, erelease.ydata

#     xmin, xmax = sorted([x1, x2])
#     ymin, ymax = sorted([y1, y2])

#     width  = xmax - xmin
#     height = ymax - ymin

#     # Store box coordinates
#     bottom_left = (min(x1, x2), min(y1, y2)) * u.pix
#     top_right = (max(x1, x2), max(y1, y2)) * u.pix
#     bl = base_drot.pixel_to_world(bottom_left[0], bottom_left[1] )
#     tr = base_drot.pixel_to_world(top_right[0] , top_right[1] )

#     coords.append([float(bl.Tx.value) ,float(bl.Ty.value),float(tr.Tx.value) ,float(tr.Ty.value)])

#     # Draw a permanent rectangle
#     rect = patches.Rectangle(
#         (xmin, ymin), width, height,
#         fill=False,
#         color='red',
#         linewidth=1.5)
#     ax.add_patch(rect)
#     rects.append(rect)
#     plt.draw()


# # Create the selector
# selector = RectangleSelector(
#     ax, onselect,
#     useblit=True,
#     button=[1],
#     minspanx=5, minspany=5,
#     spancoords='pixels',
#     interactive=True
# )
# # Keyboard shortcuts (optional):
# def on_key(event):
#     if event.key == 'd':   # delete last box
#         if rects:
#             r = rects.pop()
#             r.remove()
#             coords.pop()
#             plt.draw()

#     if event.key == 'c':   # clear all boxes
#         for r in rects:
#             r.remove()
#         rects.clear()
#         coords.clear()
#         plt.draw()


# fig.canvas.mpl_connect('key_press_event', on_key)

# plt.show()
# fig.savefig(f"{fltr}_boxed_figure.png", dpi=300, bbox_inches='tight')
# if not coords:
#     raise RuntimeError("ROI selection cancelled or failed.")

# print('-----------')
# for box in coords:
#     print(box)

# np.savetxt(f"{fltr}_boxes.csv", coords, fmt="%.6f",delimiter=',', header="xmin xmax ymin ymax")

