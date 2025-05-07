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


flare_files = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/raw/Jul10_1/'

Tx_er1=-400
Ty_er1=-400
Tx_er2=-310
Ty_er2=-475
# Load the Filter 1 image and compute its contours
ref_mg_map = sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/raw/Jul10_1/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.35.28.135_0983NB03.fits')
ref_mg_map_data = ref_mg_map.data * 1000 / ref_mg_map.meta.get('CMD_EXPT')
qs_coords =SkyCoord(Tx=(Tx_er1,Tx_er2) * u.arcsec, Ty=(Ty_er1,Ty_er2) * u.arcsec, frame=ref_mg_map.coordinate_frame)

test_map=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/raw/Jul10_1/SUT_T24_0956_000465_Lev1.0_2024-07-10T13.31.01.128_0983NB03.fits')

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
top_labels = [regions[0].label, regions[1].label, regions[2].label] #, regions[3].label, regions[4].label]

mask1 = label_img == top_labels[0]
mask2 = label_img == top_labels[1]
mask3 = label_img == top_labels[2]
#mask4 = label_img == top_labels[3]

msk=mask1+mask2+mask3#+mask4
#plt.imshow(binary_image, cmap='gray', origin='lower')
#plt.imshow(msk, cmap='gray', origin='lower')

##plt.plot(contours[1][:,1], contours[1][:,0], color='red', linewidth=1,label='Flare contour')

#plt.colorbar()
#plt.show()

contours = measure.find_contours(msk, level=0.5)
print('Number of contours:', len(contours))
hpc_coord=[]
for i in range(len(contours)):
    hpc_coord.append(ref_mg_map.pixel_to_world(contours[i][:, 1]*u.pixel, contours[i][:, 0]*u.pixel))
    #hpc_coord[i]= hpc_coord[i].transform_to(ref_mg_map.coordinate_frame) # not much chage
#'''
fig=plt.figure(figsize=(10, 10))
ax=fig.add_subplot(111, projection=test_map)
test_map.plot()
#plt.imshow(ref_mg_map_data,origin='lower' ,cmap='gray', vmin=0, vmax=16000)  # Adjust vmin/vmax as needed
#plt.plot(contours[0][:, 1], contours[0][:, 0], color='blue', linewidth=5,label='Flare contour') 
hpc=  hpc_coord[1].transform_to(test_map.coordinate_frame)
ax.plot_coord(hpc, color='red', linewidth=1,label='Flare contour')
plt.show()
#'''
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
    
    ny, nx = filter2_map.data.shape
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
    points = np.vstack((X.ravel(), Y.ravel())).T
    masks=[]
    m_msk=np.zeros((ny,nx))
    fig=plt.figure(figsize=(10, 10))
    ax=fig.add_subplot(111, projection=filter2_map)
    filter2_map.plot(axes=ax)
    for i in range(3):
        # Convert the contour coordinates to the Filter 2 image's coordinate frame
        print('Contour:', i)    
        hpc=  hpc_coord[i].transform_to(filter2_map.coordinate_frame)
        #plt.plot(contours[0][:, 1], contours[0][:, 0], color='blue', linewidth=5,label='Flare contour') 
        ax.plot_coord(hpc,color='red', linewidth=1,label='Flare contour')

        x_pix,y_pix = filter2_map.world_to_pixel(hpc_coord[i])
        # Create a polygon path from the contour in pixel space
        contour_path = Path(np.vstack([x_pix.value, y_pix.value]).T)
        # Create mask
        mask = contour_path.contains_points(points).reshape((ny, nx))
        masks.append(mask)
        m_msk+=mask

    plt.show()

   
    # Calculate the counts in Filter 2 under the Filter 1 contours
    counts_under_contours = np.sum(np.where(m_msk==1, filter2_data, 0))
    #counts_under_contours = np.sum(filter2_data[mask])

    # Overlay the Filter 1 contours on the Filter 2 image
    fig, ax = plt.subplots(figsize=(10, 10))
    

    vr=ax.imshow(filter2_data,origin='lower' ,cmap='gray', vmin=0, vmax=16000)  # Adjust vmin/vmax as needed
    fig.colorbar(vr,ax=ax , orientation='vertical')

    #for contour in contours:
    
    #plt.plot(contours[0][:, 1], contours[0][:, 0], linewidth=1.5, color='blue')
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
results_df.to_csv(('nb03_contours.csv'), index=False)