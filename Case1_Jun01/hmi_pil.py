from sunpy.map import Map
from scipy.ndimage import gaussian_filter
import numpy as np
import matplotlib.pyplot as plt
from skimage.morphology import skeletonize

hmi_map = Map("/media/adithya/Adi_disk4/SUIT_flare_work/June01_Flare/AIA_Data/HMI/HMI_cutouts/hmi.m_45s.20240601_070730_TAI.2.magnetogram.fits")  # full-disk or patch
data = hmi_map.data


smoothed_data = gaussian_filter(data, sigma=1)


# Get sign of magnetic field
sign_map = np.sign(smoothed_data)

# Calculate where the sign changes between neighboring pixels
pil_mask = np.zeros_like(data, dtype=bool)
pil_mask[:-1, :] |= (sign_map[:-1, :] * sign_map[1:, :] < 0)  # vertical edges
pil_mask[:, :-1] |= (sign_map[:, :-1] * sign_map[:, 1:] < 0)  # horizontal edges

pil_skeleton = skeletonize(pil_mask)

plt.figure(figsize=(8, 6))
hmi_map.plot(cmap='gray')  # adjust limits as needed
plt.contour(pil_skeleton, colors='red', linewidths=0.8)
plt.title("HMI Magnetogram with PIL Overlay")
#plt.colorbar()
plt.show()
