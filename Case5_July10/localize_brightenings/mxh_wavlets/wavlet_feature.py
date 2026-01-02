import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import cm
import glob
import os
import scipy.ndimage as ndi

# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------
fltr='NB04'
fits_dir = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop_/"
files = sorted(glob.glob(fits_dir + f"*{fltr}.fits"))

binary_file = f"{fltr}_transient_features_binary.fits"
binary = fits.getdata(binary_file)        # (N,Y,X)

# Optional morphological cleanup of mask
# (remove single-pixel noise)
binary = ndi.binary_opening(binary, structure=np.ones((3,3,3)))

# Make output folder
out_dir = f"{fltr}_transient_contour_pngs"
os.makedirs(out_dir, exist_ok=True)

ntime = len(files)

# ------------------------------------------------------------
# PROCESS EACH FRAME
# ------------------------------------------------------------
for i in range(ntime):
    # load original image
    img = fits.getdata(files[i])

    # current binary mask
    mask = binary[i].astype(bool)

    # --------------------------------------------------------
    # CONTOUR PLOTTING
    # --------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6,6))

    # plot base image
    im=ax.imshow(img, origin='lower', cmap='gray',vmin=6000,vmax=25000)
    cbar = fig.colorbar(im, ax=ax)

    # if any transient pixels exist
    if mask.any():
        ax.contour(mask, levels=[0.5], colors='red', linewidths=1.2)

    # title = frame index or filename
    ax.set_title(f"Frame {i:03d}")
    ax.set_xticks([])
    ax.set_yticks([])

    # save
    outname = os.path.join(out_dir, f"frame_{i:03d}.png")
    plt.savefig(outname, dpi=150, bbox_inches='tight')
    plt.close()

print("DONE: PNGs saved in", out_dir)
