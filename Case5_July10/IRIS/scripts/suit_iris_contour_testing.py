import numpy as np
import matplotlib.pyplot as plt
import sunpy.map
from skimage import measure
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib.path import Path
from astropy.io import fits
from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord

# Load two maps
iris_img = fits.open("../data/raw/flare_data/iris_l2_20240710_152851_3620108477_SJI_2796_t000.fits")
suit_map = sunpy.map.Map("/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/july10_11ut_to_18ut/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.35.28.135_0983NB03.fits")
suit_pos = get_horizons_coord(-21, suit_map.date)
suit_map.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
        
sji_imgs=iris_img[0].data
iris_map=sunpy.map.Map(sji_imgs[21], iris_img[0].header)

# Step 1: Create contours from Map A's data
threshold = 25000  # adjust as needed
contours_pix = measure.find_contours(suit_map.data, threshold)
largest_contour = max(contours_pix, key=len)
print(f"Found {len(contours_pix)} contours")
print(f"Largest contour length: {largest_contour}")

hpc_coords2 = suit_map.pixel_to_world(largest_contour[:, 1]*u.pixel, largest_contour[:, 0]*u.pixel)
#hpc_coords2 = suit_map.pixel_to_world(360*u.pixel, 350*u.pixel)
#print(hpc_coords2)
hpc_coords2= hpc_coords2.transform_to(iris_map.coordinate_frame) # not much chage
#print(hpc_coords2)
#iris_WCS= iris_map.wcs
#iris_pix= iris_WCS.world_to_pixel(hpc_coords2)
#print(iris_pix)

# Set up the figure with Map B
fig = plt.figure()
ax = plt.subplot(projection=iris_map)
iris_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
suit_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent,autoalign=True,alpha=0.5)
#plt.imshow(iris_map.data,axes=ax, origin='lower',cmap='gray')
'''
#suit_map.draw_contours(axes=ax, levels=[25000,40000],zorder=2,colors=['pink','green'],alpha=0.7)
# Loop through each contour
for contour in contours_pix:
    # contour[:, 1] = x pixels, contour[:, 0] = y pixels
    hpc_coords1 = suit_map.pixel_to_world(contour[:, 1]*u.pixel, contour[:, 0]*u.pixel)

    # Transform to Map B’s coordinate frame if different
    if suit_map.coordinate_frame != iris_map.coordinate_frame:
        hpc_coords = hpc_coords1.transform_to(iris_map.coordinate_frame)
        #print('--')

    # Plot on Map B
    ax.plot_coord(hpc_coords, color='red', linewidth=1)'''
ax.plot_coord(hpc_coords2, color='red', linewidth=1)
#ax.plot(largest_contour[:, 1]*u.pixel, largest_contour[:, 0]*u.pixel, color='red', linewidth=1)
plt.title("Contours from Map A overlaid on Map B")
plt.show()
