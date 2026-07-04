import sunpy.map
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button
from tkinter import Tk, filedialog

# Function to open file dialog and return selected FITS file
def select_fits():
    root = Tk()
    root.withdraw()  # Hide Tkinter root window
    file_path = filedialog.askopenfilename(title="Select a FITS file", filetypes=[("FITS files", "*.fits")])
    return file_path

# Load default SunPy sample images (can be replaced later)
image1 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
image2 = sunpy.map.Map(sunpy.data.sample.AIA_193_IMAGE)

# Create figure and axis
fig, ax = plt.subplots(figsize=(6, 6))
plt.subplots_adjust(bottom=0.35)  # Space for buttons and slider

# Initial blending factor
alpha_init = 0.5
blended_data = (1 - alpha_init) * image1.data + alpha_init * image2.data
im = ax.imshow(blended_data, cmap='gray', origin='lower', vmin=np.min(blended_data), vmax=np.max(blended_data))
ax.set_title("Blend Between Two Images")

# Create slider axis and slider
ax_slider = plt.axes([0.25, 0.15, 0.65, 0.03])
blend_slider = Slider(ax_slider, "Blend", 0, 1, valinit=alpha_init)

# Create button axes
ax_button1 = plt.axes([0.25, 0.05, 0.2, 0.075])
ax_button2 = plt.axes([0.55, 0.05, 0.2, 0.075])

# Create buttons
button1 = Button(ax_button1, "Load Image 1")
button2 = Button(ax_button2, "Load Image 2")

# Function to update blending
def update(val):
    alpha = blend_slider.val
    blended_data = (1 - alpha) * image1.data + alpha * image2.data
    im.set_data(blended_data)
    fig.canvas.draw_idle()

# Function to load new image1
def load_image1(event):
    global image1
    file_path = select_fits()
    if file_path:
        image1 = sunpy.map.Map(file_path)
        update(blend_slider.val)  # Update blend with new image

# Function to load new image2
def load_image2(event):
    global image2
    file_path = select_fits()
    if file_path:
        image2 = sunpy.map.Map(file_path)
        update(blend_slider.val)  # Update blend with new image

# Connect slider and buttons to functions
blend_slider.on_changed(update)
button1.on_clicked(load_image1)
button2.on_clicked(load_image2)

# Show plot
plt.show()

