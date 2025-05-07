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


# Threshold values for Filter 1 (NB03) and Filter 2 (NB04)
nb3T = 12000  # Threshold for Filter 1 (NB03)
nb4T = 10000  # Threshold for Filter 2 (NB04)

flare_files = '/Analysis/Projects_Data/Flare_Data/Nov01_Flare_Data1/Processed/Aligned_images/NB03/'


# Load the Filter 1 image and compute its contours
ref_mg_map = sunpy.map.Map('/Analysis/Projects_Data/Flare_Data/Nov01_Flare_Data1/Processed/Aligned_images/NB03/SUT_T24_1592_000636_Lev1.0_2024-11-01T02.16.26.605_0983NB03.fits')
ref_mg_map_data = ref_mg_map.data * 1000 / ref_mg_map.meta.get('CMD_EXPT')
qs_coords =SkyCoord(Tx=(-600, -530) * u.arcsec, Ty=(400,300) * u.arcsec, frame=ref_mg_map.coordinate_frame)

qs_map = ref_mg_map.submap(qs_coords)
qs_data = qs_map.data * 1000 / qs_map.meta.get('CMD_EXPT')
print(np.median(qs_data), np.mean(qs_data), np.std(qs_data))
Thresh_val=np.median(qs_data) * 3.5
print('Threshold: ', Thresh_val)
ny, nx = ref_mg_map_data.data.shape
# Create binary mask
binary_image = ref_mg_map_data > Thresh_val# True where pixel value > threshold
selem = disk(3)
binary_image=closing(binary_image, selem)
label_img = label(binary_image)
regions = sorted(regionprops(label_img), key=lambda r: r.area, reverse=True)
print('Number of regions:', len(regions))
top_labels = [regions[0].label]#, regions[1].label, regions[2].label]#, regions[3].label, regions[4].label]

mask1 = label_img == top_labels[0]
#mask2 = label_img == top_labels[1]
#mask3 = label_img == top_labels[2]
#mask4 = label_img == top_labels[3]

msk=mask1#+mask2+mask3#+mask4
#plt.imshow(binary_image, cmap='gray', origin='lower')
plt.imshow(msk, cmap='gray', origin='lower')
#plt.plot(contours[1][:,1], contours[1][:,0], color='red', linewidth=1,label='Flare contour')
plt.colorbar()
plt.show()

contours = measure.find_contours(msk, level=0.5)
print('Number of contours:', len(contours))


# Process all Filter 2 images and overlay Filter 1 contours
filter2_files = glob.glob(flare_files + '*3NB03.fits')
filter2_files = sorted(filter2_files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

# Create output directories
output_dir = 'Contour_Overlay_Results'
os.makedirs(output_dir, exist_ok=True)

# Initialize lists to store results
results = []

for filter2_file in filter2_files:
    # Load the Filter 2 image
    filter2_map = sunpy.map.Map(filter2_file)
    filter2_data = filter2_map.data * 1000 / filter2_map.meta.get('CMD_EXPT')
    scale_img= (filter2_data/255).astype(np.uint8)

    # Calculate the counts in Filter 2 under the Filter 1 contours
    counts_under_contours = np.sum(np.where(msk==1, filter2_data, 0))
    #counts_under_contours = np.sum(filter2_data[mask])

    # Overlay the Filter 1 contours on the Filter 2 image
    fig, ax = plt.subplots(figsize=(10, 10))
    

    vr=ax.imshow(filter2_data,origin='lower' ,cmap='gray', vmin=0, vmax=16000)  # Adjust vmin/vmax as needed
    fig.colorbar(vr,ax=ax , orientation='vertical')

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
    plt.title(f"Mg II k {filter2_map.date}", fontsize=16)
    plt.savefig(output_filename)
    plt.close()

    # Append results to the list
    results.append({
        'filter2_file': filter2_map.date,
        'counts_under_contours': counts_under_contours
    })

    print(f"Processed { os.path.basename(filter2_file)}: Counts under contours = {counts_under_contours:.2f}")

# Save results to a CSV file
import pandas as pd
results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(output_dir, 'results.csv'), index=False)