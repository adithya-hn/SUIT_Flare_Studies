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

ref_mg_map = sunpy.map.Map('/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Processed_Data/Aligned_images/NB03/SUT_T24_0785_000396_Lev1.0_2024-06-01T08.46.29.783_0983NB03.fits')
image = ref_mg_map.data * 1000 / ref_mg_map.meta.get('CMD_EXPT')

thresh = threshold_otsu(image)
qs_coords = SkyCoord(Tx=(-400, -310) * u.arcsec, Ty=(-400, -475) * u.arcsec, frame=ref_mg_map.coordinate_frame)
qs_map = ref_mg_map.submap(qs_coords)
qs_data = qs_map.data * 1000 / qs_map.meta.get('CMD_EXPT')
print(np.median(qs_data), np.mean(qs_data), np.std(qs_data))
Thresh_val=np.median(qs_data) * 4
print('Threshold: ', Thresh_val, thresh)

binary_image = image > thresh
plt.figure(figsize=(8, 5))
plt.imshow(binary_image, cmap='gray', origin='lower')
plt.title("Binary Image using otsu Thresholding")
plt.colorbar(label='Pixel Intensity')
plt.savefig('binary_image_otsu.png', dpi=300)
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(image.ravel(), bins=256, range=(np.min(image), np.max(image)), color='gray')
plt.title("Histogram of Image")
plt.xlabel("Pixel Intensity")
plt.ylabel("Frequency")

plt.grid(True)


plt.axvline(thresh, color='red', linestyle='--', label=f"Otsu Threshold = {thresh:.2f}")
plt.axvline(Thresh_val, color='blue', linestyle='--', label=f"4xmedian Threshold = {Thresh_val:.2f}")
plt.legend()
plt.savefig('histogram_image_4median.png', dpi=300)

plt.show()