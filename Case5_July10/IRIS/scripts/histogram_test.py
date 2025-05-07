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
from matplotlib.patches import PathPatch

from skimage.filters import threshold_otsu

ref_mg_map = sunpy.map.Map('/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/NB03/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.35.28.135_0983NB03.fits')
image = ref_mg_map.data * 1000 / ref_mg_map.meta.get('CMD_EXPT')

thresh = threshold_otsu(image)
qs_coords = SkyCoord(Tx=(-10, 90) * u.arcsec, Ty=(-55, -110) * u.arcsec, frame=ref_mg_map.coordinate_frame)
qs_map = ref_mg_map.submap(qs_coords)
qs_data = qs_map.data * 1000 / qs_map.meta.get('CMD_EXPT')
print(np.median(qs_data), np.mean(qs_data), np.std(qs_data))
Thresh_val=np.median(qs_data) * 3
print('Threshold: ', Thresh_val, thresh)

binary_image = image > Thresh_val
plt.figure(figsize=(8, 5))
plt.imshow(binary_image, cmap='gray', origin='lower')
plt.title("Binary Image using otsu Thresholding")
plt.colorbar(label='Pixel Intensity')
plt.savefig('binary_image_3md.png', dpi=300)
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(image.ravel(), bins=256, range=(np.min(image), np.max(image)), color='gray')
plt.title("Histogram of Image")
plt.xlabel("Pixel Intensity")
plt.ylabel("Frequency")

plt.grid(True)


plt.axvline(thresh, color='red', linestyle='--', label=f"Otsu Threshold = {thresh:.2f}")
plt.axvline(Thresh_val, color='blue', linestyle='--', label=f"3xmedian Threshold = {Thresh_val:.2f}")
plt.legend()
plt.savefig('histogram_image_3median.png', dpi=300)

plt.show()