import os
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import glob
import datetime
from matplotlib.path import Path
from astropy.wcs import WCS
from skimage import measure
from astropy.coordinates import SkyCoord
import astropy.units as u


# Threshold values for Filter 1 (NB03) and Filter 2 (NB04)
nb3T = 12000  # Threshold for Filter 1 (NB03)
nb4T = 10000  # Threshold for Filter 2 (NB04)

# Define the specific timestamp for which contours will be drawn
specific_timestamp = "2024-07-10T15.35.28.135"  # Replace with your desired timestamp

# Paths to Filter 1 and Filter 2 data
filter1_fold = '/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_/Aligned_images/NB03/'
filter2_fold = '/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_/Aligned_images/NB03/'

# Find the specific file for Filter 1
filter1_files = glob.glob(filter1_fold + '*3NB03.fits')
filter1_files = sorted(filter1_files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

specific_file_filter1 = None
for file in filter1_files:
    file_timestamp = os.path.basename(file).split('_')[5]
    #print(file_timestamp)
    if file_timestamp == specific_timestamp:
        specific_file_filter1 = file
        break

if specific_file_filter1 is None:
    print(f"No file found for the specific timestamp: {specific_timestamp}")
    exit()

# Load the Filter 1 image and compute its contours
filter1_map = sunpy.map.Map(specific_file_filter1)

filter1_data = filter1_map.data * 1000 / filter1_map.meta.get('CMD_EXPT')
qs_coord=SkyCoord(Tx=(-10, 90) * u.arcsec, Ty=(-20, -120) * u.arcsec, frame=filter1_map.coordinate_frame)
qs_map=filter1_map.submap(qs_coord)
qs_data = qs_map.data * 1000 / qs_map.meta.get('CMD_EXPT')
print(np.median(qs_data),np.mean(qs_data),np.std(qs_data))
Thresh_val=np.median(qs_data)*3
print('Threshold: ',Thresh_val)

filter1_data[:,0:100]=0
filter1_data[:,300:500]=0
#print(filter1_data)
msk = np.where(filter1_data > Thresh_val, 1, 0)
print(np.sum(msk))

#breakpoint()

# Create binary mask
binary_image = filter1_data > Thresh_val  # True where pixel value > threshold
contours = measure.find_contours(binary_image, level=0.5)  # 0.5 for binary images
# Select the largest (outermost) contour by length
largest_contour = max(contours, key=len)
print(f"Found {len(contours)} contours")


# Process all Filter 2 images and overlay Filter 1 contours
filter2_files = glob.glob(filter2_fold + '*3NB03.fits')
filter2_files = sorted(filter2_files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

# Create output directories
output_dir = 'NB03_cont_imgs'
os.makedirs(output_dir, exist_ok=True)

# Initialize lists to store results
results = []

for filter2_file in filter2_files:
    # Load the Filter 2 image
    filter2_map = sunpy.map.Map(filter2_file)
    filter2_data = filter2_map.data * 1000 / filter2_map.meta.get('CMD_EXPT')
    QS_coord=SkyCoord(Tx=(-10, 90) * u.arcsec, Ty=(-20, -120) * u.arcsec, frame=filter2_map.coordinate_frame)
    QS_map=filter2_map.submap(QS_coord)
    QS_data = QS_map.data * 1000 / QS_map.meta.get('CMD_EXPT')

    
    # Calculate the counts in Filter 2 under the Filter 1 contours
    counts_under_contours = np.mean(np.where(msk==1, filter2_data, 0))
    #counts_under_contours = np.sum(filter2_data[mask])

    # Overlay the Filter 1 contours on the Filter 2 image
    fig, ax = plt.subplots(figsize=(10, 10))
    vr=ax.imshow(filter2_data,origin='lower' ,cmap='gray', vmin=0, vmax=14000)  # Adjust vmin/vmax as needed
    fig.colorbar(vr,ax=ax , orientation='vertical')

    #for contour in contours:
    
    plt.plot(largest_contour[:, 1], largest_contour[:, 0], linewidth=1.5, color='blue')

    # Add text to display the counts
    ax.text(50, 50, f"Counts under contours: {counts_under_contours:.2f}", color='white', fontsize=12)

    # Save the plot
    output_filename = os.path.join(output_dir, os.path.basename(filter2_file)[:-5] + '_overlay.jpg')
    #ax.set_colorbars()
    
    plt.savefig(output_filename)
    plt.close()

    # Append results to the list
    results.append({
        'filter2_file': filter2_map.date,
        'mean_counts_under_contours': counts_under_contours,
        'QS_mean_counts': np.mean(QS_data)
    })

    print(f"Processed { os.path.basename(filter2_file)}: Counts under contours = {counts_under_contours:.2f}")

# Save results to a CSV file
import pandas as pd
results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(output_dir, 'results.csv'), index=False)