import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import glob
import matplotlib.animation as animation
import matplotlib.image as mpimg
import numpy as np
from datetime import datetime

path_dir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/OVER_PLOT_CONTOURS/Contour_imgs/1600/'
# Replace 'my_image.png' with the actual path to your image file

imgs=sorted(glob.glob(path_dir+'*.jpg'))



# --- Set up plot ---
#fig, (ax_lc, ax_img) = plt.subplots(1, 2, figsize=(15, 6))
fig, (ax_lc, ax_img) = plt.subplots(1, figsize=(15, 6), gridspec_kw={'width_ratios': [1.7, 1.3]})

ax2=ax_lc.twinx()
#ax_lc.set_aspect(1.5) 
# Plot light curve

ax_lc.set_ylabel("Counts")
#ax_lc.legend()
ax_lc.set_title("Light Curve")

# Load first image
img = mpimg.imread(imgs[0])
img_disp = ax_img.imshow(img)
ax_img.axis('off')


# --- Update function ---
def update(i):


    img = mpimg.imread(imgs[i])
    img_disp.set_data(img)
   

    return img_disp

# --- Animate ---
ani = animation.FuncAnimation(fig, update, frames=len(imgs), interval=300, blit=False)
#plt.gca().xaxis.set_major_formatter(time_formatter)
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
plt.tight_layout()
ani.save("lightcurve_with_sliding_image.mp4", writer='ffmpeg',dpi=300, fps=5)
plt.show()