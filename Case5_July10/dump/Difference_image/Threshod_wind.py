import sunpy.map
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider

# Load a sample SunPy map (Replace with your own FITS file)
sample_data =sunpy.map.Map('/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/NB08/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.21.40.774_0973NB08.fits')
smap = sunpy.map.Map(sample_data)

# Create figure and axis
fig, ax = plt.subplots(figsize=(6, 6))
plt.subplots_adjust(bottom=0.25)  # Space for slider

# Initial threshold
initial_threshold = np.max(smap.data) / 2

# Display image
masked_data = np.where(smap.data > initial_threshold, smap.data, np.nan)
im = ax.imshow(masked_data, cmap='gray', origin='lower', vmin=np.min(smap.data), vmax=np.max(smap.data))
ax.set_title("Adjust Threshold with Slider")

# Create slider axis and slider
ax_slider = plt.axes([0.25, 0.1, 0.65, 0.03])  # [left, bottom, width, height]
threshold_slider = Slider(ax_slider, "Threshold", np.min(smap.data), np.max(smap.data), valinit=initial_threshold)

# Update function for slider
def update(val):
    threshold = threshold_slider.val
    masked_data = np.where(smap.data > threshold, smap.data, np.nan)
    im.set_data(masked_data)
    fig.canvas.draw_idle()

# Connect slider to update function
threshold_slider.on_changed(update)

# Show plot
plt.show()
