import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
import glob
from matplotlib import rcParams
#rcParams['font.family'] = 'Times New Roman'  
plt.rcParams.update({"text.usetex": True,})
#-----------------------------------------------------------
# USER INPUT
#-----------------------------------------------------------
image_folder = "nb4_imgs/"             # path to your image directory
pattern = "*.png"
n_rows, n_cols = 3, 5

# default crop (left, right, top, bottom)
default_crop = (200,220,8,125)  

# special crop rules
first_column_crop = (200,220,8,125)  #(110,220,30,125)    # example: more crop on left for column-1
last_row_crop    = (200,220,8,125)  #(200,220,30,70)     # example: more crop on top for last row

fig_size = (10,6)
save_name = "c4_nb04_transients.pdf"
#-----------------------------------------------------------

def crop_image(im, crop):
    l,r,t,b = crop
    w,h = im.size
    return im.crop((l, t, w-r, h-b))

files = sorted(glob.glob(image_folder + pattern))[:n_rows*n_cols]
images = []

for idx, file in enumerate(files):
    img = Image.open(file)

    row = idx // n_cols
    col = idx % n_cols

    crop = default_crop

    # apply special crop if in last row
    if row == n_rows-1:
        crop = last_row_crop

    # apply special crop if in first column
    if col == 0:
        crop = first_column_crop

    img = crop_image(img, crop)
    images.append(img)


def remove_axis(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_visible(True)          # keep axes but hide frame
    for spine in ax.spines.values():
        spine.set_visible(False)

# plotting
fig, axes = plt.subplots(n_rows, n_cols, figsize=fig_size)
plt.subplots_adjust(left=0, right=.5, top=0.49, bottom=0,wspace=0, hspace=0)   # reduce gaps

for i, ax in enumerate(axes.flatten()):
    if i < len(images):
        ax.imshow(images[i])
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_frame_on(False)      # <<< IMPORTANT
        ax.set_axis_off()           # <<< IMPORTANT

        ax.text(0.2, 0.94, f"{i+1})",
                transform=ax.transAxes,
                fontsize=6,
                va='top',
                ha='left',
                family= " serif",
                font='Times New Roman',
                color='k')

    else:
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_frame_on(False)
        ax.set_axis_off()

plt.savefig(save_name, dpi=400, bbox_inches='tight')