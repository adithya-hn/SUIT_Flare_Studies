

import os
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import glob
import datetime
from matplotlib.path import Path
from astropy.wcs import WCS
from skimage import measure
from matplotlib.widgets import Slider

'''
#Ca II h set_1
suit_map1=sunpy.map.Map('/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/NB08/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.21.40.774_0973NB08.fits')
suit_map2=sunpy.map.Map('/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/NB08/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.34.38.918_0983NB08.fits')
'''

#Ca II h set_2
suit_map1=sunpy.map.Map('/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/NB08/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.08.45.268_0983NB08.fits')
suit_map2=sunpy.map.Map('/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/NB08/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.15.12.550_0973NB08.fits')

#local peak  /Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/NB08/SUT_T24_0956_000465_Lev1.0_2024-07-10T14.12.10.637_0983NB08.fits
#smooth drop /Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/NB08/SUT_T24_0956_000465_Lev1.0_2024-07-10T14.43.12.206_0983NB08.fits
#continuum   /Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/NB08/SUT_T24_0956_000465_Lev1.0_2024-07-10T14.50.30.225_0983NB08.fits
#Edge        /Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/NB08/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.05.06.261_0983NB08.fits
#dip         /Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/NB08/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.08.45.268_0983NB08.fits
#small peek  /Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/NB08/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.15.12.550_0973NB08.fits


suit_data1=suit_map1.data
suit_data2=suit_map2.data

diffed_data=suit_data2-suit_data1
diff_map=sunpy.map.Map(diffed_data,suit_map1.meta)
#diff_map.peek(cmap='coolwarm',colorbar=True)

# Create figure and axis
fig, ax = plt.subplots(figsize=(6, 6))
plt.subplots_adjust(bottom=0.25)  # Space for slider

# Initial threshold
initial_threshold = np.max(diff_map.data) / 2

# Display image
masked_data = np.where(diff_map.data > initial_threshold, diff_map.data, np.nan)
im = ax.imshow(masked_data, cmap='coolwarm', origin='lower', vmin=np.min(diff_map.data), vmax=np.max(diff_map.data))
ax.set_title("Difference image")

# Create slider axis and slider
ax_slider = plt.axes([0.25, 0.1, 0.65, 0.03])  # [left, bottom, width, height]
threshold_slider = Slider(ax_slider, "Threshold", np.min(diff_map.data), np.max(diff_map.data), valinit=initial_threshold)

# Update function for slider
def update(val):
    threshold = threshold_slider.val
    masked_data = np.where(diff_map.data > threshold, diff_map.data, np.nan)
    im.set_data(masked_data)
    fig.canvas.draw_idle()

# Connect slider to update function
threshold_slider.on_changed(update)

# Show plot
plt.colorbar(im)
plt.show()

