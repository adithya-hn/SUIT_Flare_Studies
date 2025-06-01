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
import cv2
import imageio
from skimage import filters, measure
from skimage.measure import label, regionprops
from skimage.morphology import disk, closing
from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord



filter2_fold = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/cropped/crop_fits/NB03/'

Tx1_qs1,Ty1_qs1=-330,-480
Tx2_qs1,Ty2_qs1=-280,-450


# Load the Filter 1 image and compute its contours
ref_mg_map = sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/cropped/crop_fits/NB03/SUT_T24_0785_000396_Lev1.0_2024-06-01T08.46.29.783_0983NB03.fits')
mgk_map=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/cropped/crop_fits/NB03/SUT_T24_0785_000396_Lev1.0_2024-06-01T07.10.22.791_0973NB03.fits')
ref_mg_map_data = ref_mg_map.data * 1000 / ref_mg_map.meta.get('CMD_EXPT')

ref_pos = get_horizons_coord(-21, ref_mg_map.date)
ref_mg_map.meta.update(get_observer_meta(ref_pos, rsun=ref_pos.rsun))
mgk_map.meta.update(get_observer_meta(ref_pos, rsun=ref_pos.rsun))
qs_coords = SkyCoord(Tx=(Tx1_qs1,Tx2_qs1) * u.arcsec, Ty=(Ty1_qs1,Ty2_qs1) * u.arcsec, frame=mgk_map.coordinate_frame)

qs_coords1 = SkyCoord(Tx=(Tx1_qs1,Tx2_qs1) * u.arcsec, Ty=(Ty1_qs1,Ty2_qs1) * u.arcsec, frame=mgk_map.coordinate_frame)


qs_map = mgk_map.submap(qs_coords)
qs_data = qs_map.data * 1000 / qs_map.meta.get('CMD_EXPT')
print(np.median(qs_data), np.mean(qs_data), np.std(qs_data))
Thresh_val=np.median(qs_data) * 3
print('Threshold: ', Thresh_val)
ny, nx = ref_mg_map_data.data.shape
# Create binary mask
binary_image = ref_mg_map_data > Thresh_val# True where pixel value > threshold
selem = disk(3)
binary_image=closing(binary_image, selem)
label_img = label(binary_image)
regions = sorted(regionprops(label_img), key=lambda r: r.area, reverse=True)
print('Number of regions:', len(regions))
top_labels = [regions[0].label, regions[1].label, regions[2].label]#, regions[3].label] #, regions[4].label]

mask1 = label_img == top_labels[0]
mask2 = label_img == top_labels[1]
mask3 = label_img == top_labels[2]
#mask4 = label_img == top_labels[3]

msk=mask1+mask2#+mask3#+mask4
#plt.imshow(binary_image, cmap='gray', origin='lower')
plt.imshow(msk, cmap='gray', origin='lower')
#plt.plot(contours[1][:,1], contours[1][:,0], color='red', linewidth=1,label='Flare contour')
plt.colorbar()
plt.show()

contours = measure.find_contours(msk, level=0.5)


# Process all Filter 2 images and overlay Filter 1 contours
filter2_files = glob.glob(filter2_fold + '*3NB03.fits')
filter2_files = sorted(filter2_files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

# Create output directories
output_dir = 'Contour_Overlay_Results'
os.makedirs(output_dir, exist_ok=True)

fig=plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection=mgk_map)
#ca_map.plot(axes=ax, cmap='gray', vmin=0, vmax=5000)
im=plt.imshow(mgk_map.data, origin='lower', cmap='gray', vmin=0, vmax=30000)
mgk_map.draw_quadrangle(qs_coords1, axes=ax, edgecolor="blue", linestyle="-", linewidth=2, label='Mg-II h Flare contour')
hpc_coord=[]
for contour in contours:
    ax.plot(contour[:, 1], contour[:, 0], linewidth=1.5, color='red')
    hpc_coord.append(ref_mg_map.pixel_to_world(contour[:, 1]*u.pixel, contour[:, 0]*u.pixel))
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im, cax=cax) # Similar to fig.colorbar(im, cax = cax)
plt.savefig('Mg_II_k_map.jpg')
plt.show()


# Initialize lists to store results
results = []

for filter2_file in filter2_files:
    # Load the Filter 2 image
    filter2_map = sunpy.map.Map(filter2_file)
    filter2_map.meta.update(get_observer_meta(ref_pos, rsun=ref_pos.rsun))
    qs_coords2= SkyCoord(Tx=(Tx1_qs1,Tx2_qs1) * u.arcsec, Ty=(Ty1_qs1,Ty2_qs1) * u.arcsec, frame=filter2_map.coordinate_frame)
    qs_submap = filter2_map.submap(qs_coords2)
    ca_qs_data = qs_submap.data * 1000 / qs_submap.meta.get('CMD_EXPT')
    filter2_data = filter2_map.data * 1000 / filter2_map.meta.get('CMD_EXPT')
    qs_counts=np.mean(ca_qs_data)
    qs_area=(ca_qs_data).shape[0]*(ca_qs_data).shape[1]

    # Calculate the counts in Filter 2 under the Filter 1 contours
    counts_under_contours = np.sum(np.where(msk==1, filter2_data, 0))
    #counts_under_contours = np.sum(filter2_data[mask])

    # Overlay the Filter 1 contours on the Filter 2 image

    fig= plt.figure(figsize=(10, 10))
    ax=fig.add_subplot(111,projection=filter2_map)
    filter2_map.plot(vmin=0, vmax=30000)

    for i in range(len(hpc_coord)):
        hpc=  hpc_coord[i].transform_to(filter2_map.coordinate_frame)
        ax.plot_coord(hpc,color='red', linewidth=1,label='Flare contour')
    


    # Add text to display the counts
    ax.text(50, 50, f"Counts under contours: {counts_under_contours:.2f}", color='white', fontsize=12)
    filter2_map.draw_quadrangle(qs_coords2,axes=ax,edgecolor="blue",linestyle="-",linewidth=1,label='QS region',alpha=0.5)
    plt.colorbar()

    # Save the plot
    output_filename = os.path.join(output_dir, os.path.basename(filter2_file)[:-5] + '_overlay.jpg')
    #ax.set_colorbars()
    plt.title(f"Mg II k {filter2_map.date}", fontsize=16)
    plt.savefig(output_filename)
    plt.close()

    results.append({
        'filter2_file': filter2_map.date,
        'total_counts_under_contours': counts_under_contours,
        'qs_mean_counts':qs_counts,
        'contour_area':np.count_nonzero(msk),
        'QS_area':qs_area

    })


    print(f"Processed { os.path.basename(filter2_file)}: Counts under contours = {counts_under_contours:.2f}")

# Save results to a CSV file
import pandas as pd
results_df = pd.DataFrame(results)
results_df.to_csv(('nb03_contours.csv'), index=False)