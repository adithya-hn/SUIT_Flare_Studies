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
import matplotlib.patches as patches

filter2_fold = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/processed/Aligned_images/NB08/'

Tx1_qs1,Ty1_qs1=-10,-10
Tx2_qs1,Ty2_qs1=80,-100

# Load the Filter 1 image and compute its contours
ref_mg_map = sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/processed/Aligned_images/NB03/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.35.28.135_0983NB03.fits')
ref_mg_map_data = ref_mg_map.data * 1000 / ref_mg_map.meta.get('CMD_EXPT')
qs_coords = SkyCoord(Tx=(Tx1_qs1,Tx2_qs1) * u.arcsec, Ty=(Ty1_qs1,Ty2_qs1) * u.arcsec,  frame=ref_mg_map.coordinate_frame)


#-------------
ca_map = sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/processed/Aligned_images/NB08/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.34.38.918_0983NB08.fits') #not sure about this
ca_map_data = ca_map.data * 1000 / ca_map.meta.get('CMD_EXPT')
qs_coords1 =  SkyCoord(Tx=(Tx1_qs1,Tx2_qs1) * u.arcsec, Ty=(Ty1_qs1,Ty2_qs1) * u.arcsec,  frame=ca_map.coordinate_frame)

#-------------

qs_map = ref_mg_map.submap(qs_coords)
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
top_labels = [regions[0].label]#, regions[1].label, regions[2].label, regions[3].label] #, regions[4].label]

mask1 = label_img == top_labels[0]
#mask2 = label_img == top_labels[1]
#mask3 = label_img == top_labels[2]
#mask4 = label_img == top_labels[3]

msk=mask1#+mask2#+mask3#+mask4
#plt.imshow(binary_image, cmap='gray', origin='lower')
plt.imshow(msk, cmap='gray', origin='lower')
#plt.plot(contours[1][:,1], contours[1][:,0], color='red', linewidth=1,label='Flare contour')
plt.colorbar()
plt.show()

contours = measure.find_contours(msk, level=0.5)

print('Total Contours:', len(contours))
#---ca_map
fig=plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection=ca_map)
#ca_map.plot(axes=ax, cmap='gray', vmin=0, vmax=5000)
im=plt.imshow(ca_map_data, origin='lower', cmap='gray', vmin=0, vmax=5000)
ca_map.draw_quadrangle(qs_coords1, axes=ax, edgecolor="blue", linestyle="-", linewidth=2, label='Mg-II k Flare contour')
for contour in contours:
    ax.plot(contour[:, 1], contour[:, 0], linewidth=1.5, color='red')
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im, cax=cax) # Similar to fig.colorbar(im, cax = cax)
plt.savefig('Ca_II_H_map.jpg')
plt.show()

# Process all Filter 2 images and overlay Filter 1 contours
filter2_files = glob.glob(filter2_fold + '*3NB08.fits')
filter2_files = sorted(filter2_files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

# Create output directories
output_dir = 'Contour_Overlay_Results'
os.makedirs(output_dir, exist_ok=True)

# Initialize lists to store results
results = []

for filter2_file in filter2_files:
    # Load the Filter 2 image
    filter2_map = sunpy.map.Map(filter2_file)
    qs_submap = filter2_map.submap(qs_coords)
    mg_qs_data = qs_submap.data * 1000 / qs_submap.meta.get('CMD_EXPT')
    filter2_data = filter2_map.data * 1000 / filter2_map.meta.get('CMD_EXPT')
    qs_counts=np.mean(mg_qs_data)
    
    # Calculate the counts in Filter 2 under the Filter 1 contours
    counts_under_contours = np.sum(np.where(msk==1, filter2_data, 0))
    #counts_under_contours = np.sum(filter2_data[mask])
    pix1,pix2=filter2_map.world_to_pixel(qs_coords)
    print(pix1,pix2)
    #rpx=rct_pix.x.value
    #rpy=rct_pix.y.value

    #x_min = min(pix1.x.value, pix2.x.value)
    #y_min = min(pix1.y.value, pix2.y.value)
    #width = abs(pix1.x.value - pix2.x.value)
    #height = abs(pix1.y.value - pix2.y.value)

    # Overlay the Filter 1 contours on the Filter 2 image
    fig, ax = plt.subplots(figsize=(10, 10))
    

    vr=ax.imshow(filter2_data,origin='lower' ,cmap='gray', vmin=0, vmax=5000)  # Adjust vmin/vmax as needed
    fig.colorbar(vr,ax=ax , orientation='vertical')
    #rect = patches.Rectangle((x_min, y_min), width, height,  linewidth=2, edgecolor='blue', facecolor='none')
    #ax.add_patch(rect)

    #for contour in contours:
    
    plt.plot(contours[0][:, 1], contours[0][:, 0], linewidth=1.5, color='blue')
    #plt.plot(contours[1][:, 1], contours[1][:, 0], linewidth=1.5, color='blue')
    #plt.plot(contours[2][:, 1], contours[2][:, 0], linewidth=1.5, color='blue')
    #plt.plot(contours[3][:, 1], contours[3][:, 0], linewidth=1.5, color='blue')

    # Add text to display the counts
    ax.text(50, 50, f"Counts under contours: {counts_under_contours:.2f}", color='white', fontsize=12)

    # Save the plot
    output_filename = os.path.join(output_dir, os.path.basename(filter2_file)[:-5] + '_overlay.jpg')
    #ax.set_colorbars()
    plt.title(f"ca II h {filter2_map.date}", fontsize=16)
    plt.savefig(output_filename)
    plt.close()

    # Append results to the list
    results.append({
        'filter2_file': filter2_map.date,
        'total_counts_under_contours': counts_under_contours,
        'qs_mean_counts':qs_counts,
        'contour_area':np.count_nonzero(msk)

    })

    print(f"Processed { os.path.basename(filter2_file)}: Counts under contours = {counts_under_contours:.2f}")

# Save results to a CSV file
import pandas as pd
results_df = pd.DataFrame(results)
results_df.to_csv(('nb08_contours.csv'), index=False)