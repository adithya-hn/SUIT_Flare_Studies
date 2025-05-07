import os
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import glob
import datetime
from matplotlib.path import Path
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
#import cv2
import imageio
from skimage import filters, measure
from skimage.measure import label, regionprops
from skimage.morphology import disk, closing
from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord
from sunpy.time import parse_time
from astropy.time import Time


flare_files = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/raw/'

rct_Tx1,rct_Ty1=-581,-500
rct_Tx2,rct_Ty2=-225,-174

Tx1_qs1,Ty1_qs1=-340,-520
Tx2_qs1,Ty2_qs1=-290,-470

Tx1_qs2,Ty1_qs2,Tx2_qs2,Ty2_qs2=-275,-240,-225,-225

Tx1_qs3,Ty1_qs3,Tx2_qs3,Ty2_qs3=-555,-490,-525,-465

# Load the Filter 1 image and compute its contours
#ref_mg_map = sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/raw/SUT_T24_0785_000396_Lev1.0_2024-06-01T08.46.29.783_0983NB03.fits')
ref_mg_map=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/raw/SUT_T24_0785_000396_Lev1.0_2024-06-01T07.14.12.134_0973NB03.fits')
ref_pos = get_horizons_coord(-21, ref_mg_map.date)
ref_mg_map.meta.update(get_observer_meta(ref_pos, rsun=ref_pos.rsun))

ref_mg_map_data = ref_mg_map.data * 1000 / ref_mg_map.meta.get('CMD_EXPT')
qs_coords1 =SkyCoord(Tx=(Tx1_qs1,Tx2_qs1) * u.arcsec, Ty=(Ty1_qs1,Ty2_qs1) * u.arcsec, frame=ref_mg_map.coordinate_frame)
qs_coords2 =SkyCoord(Tx=(Tx1_qs2,Tx2_qs2) * u.arcsec, Ty=(Ty1_qs2,Ty2_qs2) * u.arcsec, frame=ref_mg_map.coordinate_frame)
qs_coords3 =SkyCoord(Tx=(Tx1_qs3,Tx2_qs3) * u.arcsec, Ty=(Ty1_qs3,Ty2_qs3) * u.arcsec, frame=ref_mg_map.coordinate_frame)
#qs_coords2=SkyCoord(Tx=(Tx_qs1,Tx_qs2) * u.arcsec, Ty=(Ty_qs1,Ty_qs2) * u.arcsec, frame=ref_mg_map.coordinate_frame)
fig=plt.figure()
ax=fig.add_subplot(111,projection=ref_mg_map)
ref_mg_map.plot(vmin=1000,vmax=25000)
rect_coords=SkyCoord(Tx=(rct_Tx1,rct_Tx2) * u.arcsec, Ty=(rct_Ty1,rct_Ty2) * u.arcsec, frame=ref_mg_map.coordinate_frame)
ref_mg_map.draw_quadrangle(rect_coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='selection region')
ref_mg_map.draw_quadrangle(qs_coords1,axes=ax,edgecolor='yellow',linestyle="-",linewidth=2,label='QS region 1')
ref_mg_map.draw_quadrangle(qs_coords2,axes=ax,edgecolor='yellow',linestyle="-",linewidth=2,label='QS region 2')
ref_mg_map.draw_quadrangle(qs_coords3,axes=ax,edgecolor='yellow',linestyle="-",linewidth=2,label='QS region 3')
plt.show()
selcted_map=ref_mg_map.submap(rect_coords)
plt.hist((selcted_map.data).flatten(),bins=1000)
plt.show()

#'''
norm_sel_rct=selcted_map.data*1000/selcted_map.meta.get('CMD_EXPT')
median_=np.median(norm_sel_rct)
std_=np.std(norm_sel_rct)
quiet_sun_mask = norm_sel_rct < (median_)
print('Region thresh',median_,std_)
quiet_sun_level = np.median(norm_sel_rct[quiet_sun_mask])
qs_sigma=np.std(norm_sel_rct[quiet_sun_mask])
print('QS_level',quiet_sun_level,qs_sigma)
fig,axs=plt.subplots(1,2)
axs[0].imshow(norm_sel_rct,origin='lower')
axs[1].imshow(quiet_sun_mask,origin='lower')
plt.show()
#'''

test_map=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/raw/SUT_T24_0785_000396_Lev1.0_2024-06-01T07.14.12.134_0973NB03.fits')

norm_sel_rct=selcted_map.data*1000/selcted_map.meta.get('CMD_EXPT')

qs_map1 = ref_mg_map.submap(qs_coords1)
qs_map2 = ref_mg_map.submap(qs_coords2)
qs_map3 = ref_mg_map.submap(qs_coords3)
qs_data1 = qs_map1.data * 1000 / qs_map1.meta.get('CMD_EXPT')
qs_data2 = qs_map2.data * 1000 / qs_map2.meta.get('CMD_EXPT')
qs_data3 = qs_map3.data * 1000 / qs_map3.meta.get('CMD_EXPT')

print(np.median(qs_data1), np.mean(qs_data1), np.std(qs_data1))
print(np.median(qs_data2), np.mean(qs_data2), np.std(qs_data1))
print(np.median(qs_data3), np.mean(qs_data1), np.std(qs_data1))
qs_thresh=(np.median(qs_data1)+np.median(qs_data2)+np.median(qs_data3))/3
print('QS Threshold: ',qs_thresh )

quiet_sun_mask = norm_sel_rct < (qs_thresh*1.5)

quiet_sun_level = np.median(norm_sel_rct[quiet_sun_mask])
qs_sigma=np.std(norm_sel_rct[quiet_sun_mask])
print('QS_level',quiet_sun_level,qs_sigma)

fig,axs=plt.subplots(1,2)
axs[0].imshow(norm_sel_rct,origin='lower')
axs[1].imshow(quiet_sun_mask,origin='lower')
plt.show()

Thresh_val=qs_thresh*4 #np.median(qs_data1) * 3
print('Threshold: ', Thresh_val)
ny, nx = ref_mg_map_data.data.shape
# Create binary mask
binary_image = ref_mg_map_data > Thresh_val# True where pixel value > threshold
selem = disk(3)
binary_image=closing(binary_image, selem)
label_img = label(binary_image)
regions = sorted(regionprops(label_img), key=lambda r: r.area, reverse=True)
print('Number of regions:', len(regions))

# Process all Filter 2 images and overlay Filter 1 contours
filter2_files = glob.glob(flare_files + '*3NB03.fits')
filter1_files = glob.glob(flare_files + '*3NB08.fits')
filter2_files = sorted(filter2_files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
filter1_files = sorted(filter1_files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

nb3_time_array=[]
for f in range(len(filter2_files)):
    nb3_time_array.append(datetime.datetime.strptime(os.path.basename(filter2_files[f]).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
nb3_time_array=Time(parse_time(nb3_time_array))
# Create output directories
output_dir = 'Contour_Overlay_Results'
os.makedirs(output_dir, exist_ok=True)

# Initialize lists to store results
results = []

for filter1_file in filter1_files:
    # Load the Filter 2 image
    filter1_map=sunpy.map.Map(filter1_file)
    
    ca_time=Time(parse_time(filter1_map.date))
    idx=np.argmin(np.abs(nb3_time_array - ca_time))
    #print(filter1_file,filter2_files[idx])
    filter2_map = sunpy.map.Map(filter2_files[idx])
    filter2_pos = get_horizons_coord(-21, filter2_map.date)
    filter2_map.meta.update(get_observer_meta(filter2_pos, rsun=filter2_pos.rsun))
    selected_reg=filter2_map.submap(rect_coords)
    #ca_region=filter1_map.submap(rect_coords)
    rct_data = selected_reg.data * 1000 / filter2_map.meta.get('CMD_EXPT')
    filter2_data = filter2_map.data * 1000 / filter2_map.meta.get('CMD_EXPT')
    #filter1_data = ca_region.data * 1000 / filter1_map.meta.get('CMD_EXPT')
    Norm_map=sunpy.map.Map(rct_data,selected_reg.fits_header)
    #Norm_ca_map=sunpy.map.Map(ca_region.data,ca_region.fits_header)
    norm_map=sunpy.map.Map(filter2_data,filter2_map.fits_header)
    qs_data1 = norm_map.submap(qs_coords1).data
    qs_data2 = norm_map.submap(qs_coords2).data
    qs_data3 = norm_map.submap(qs_coords3).data
    #print(np.median(qs_data1), np.mean(qs_data1), np.std(qs_data1))
    #print(np.median(qs_data2), np.mean(qs_data2), np.std(qs_data1))
    #print(np.median(qs_data3), np.mean(qs_data1), np.std(qs_data1))
    qs_thresh=(np.median(qs_data1))#+np.median(qs_data2)+np.median(qs_data3))/3
    print('QS Threshold: ',qs_thresh )

    plg_msk_b = (Norm_map.data > qs_thresh*1.5) & (Norm_map.data < qs_thresh*4) # True where pixel value > threshold
    kernel = disk(3)
    plg_msk=closing(plg_msk_b, kernel)
    plg_cont = measure.find_contours(plg_msk, level=0.5)
    print('Total plages contour: ',len(plg_cont))
    #plg_hpc=[]
    



    flare_msk_b=Norm_map.data > qs_thresh*4
    flare_msk=closing(flare_msk_b,kernel)
    flare_cont = measure.find_contours(flare_msk, level=0.5)
    print(' Total flare contours: ', len(flare_cont))

    fig=plt.figure(figsize=(10, 10))
    ax=fig.add_subplot(111, projection=Norm_map)
    Norm_map.plot(axes=ax,autoalign=True,vmin=1000,vmax=16000)
    th_lv=[qs_thresh*1.5,(qs_thresh)*4]
    for pc in range(len(plg_cont)):
        plg_hpc=ref_mg_map.pixel_to_world(plg_cont[pc][:, 1]*u.pixel, plg_cont[pc][:, 0]*u.pixel)

    Norm_map.draw_contours(axes=ax, levels=th_lv,zorder=1,colors=['yellow','red'])
    Norm_map.draw_quadrangle(qs_coords1,axes=ax,edgecolor='green',linestyle="-",linewidth=2,label='QS region 1')
    Norm_map.draw_quadrangle(qs_coords2,axes=ax,edgecolor='green',linestyle="-",linewidth=2,label='QS region 2')
    Norm_map.draw_quadrangle(qs_coords3,axes=ax,edgecolor='green',linestyle="-",linewidth=2,label='QS region 3')

    output_filename = os.path.join(output_dir, os.path.basename(filter2_files[idx])[:-5] + '_overlay.jpg')
    #ax.set_colorbars()
    plt.title(f"Mg II k {filter2_map.date}", fontsize=16)
    plt.savefig(output_filename)
    plt.close()

    # Append results to the list
    Thresh_alned_data=np.where(Norm_map.data>qs_thresh*4,Norm_map.data,0)
    flare_area=np.count_nonzero(Thresh_alned_data)
    flare_count=np.sum(Thresh_alned_data)

    plage_th=np.where(((Norm_map.data>qs_thresh*2)&(Norm_map.data<qs_thresh*4)),Norm_map.data,0)
    plage_area=np.count_nonzero(plage_th)
    plage_count=np.sum(plage_th)


    results.append({'filter2_file': filter2_map.date,
                     'flare_tot_count': flare_count,
                     'flare_tot_area':flare_area,
                     'plage_tot_count':plage_count,
                     'plage_tot_area':plage_area,
                     'qs_threshold':qs_thresh,
                     'qs_smpl1':np.median(qs_data1),
                     'qs_smpl2':np.median(qs_data2),
                     'qs_smpl2':np.median(qs_data3)})

    #rint(f"Processed { os.path.basename(filter2_file)}: Counts under contours = {counts_under_contours:.2f}")

# Save results to a CSV file
import pandas as pd
results_df = pd.DataFrame(results)
results_df.to_csv(('nb03_contours.csv'), index=False)