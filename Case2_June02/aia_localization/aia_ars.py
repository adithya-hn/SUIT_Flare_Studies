"""
Find bright active regions in full-disc AIA 131 FITS images, draw boxes around them,
and save results.

Dependencies:
  pip install sunpy astropy scikit-image scipy matplotlib pandas
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.ndimage import gaussian_filter
from skimage.measure import label, regionprops
import pandas as pd

import sunpy.map
import astropy.units as u

# --- User parameters ---
#fits_folder = Path("/media/adithya/Adi_disk4/SUIT_flare_work/case1_Jun01/data/aia/aia_fd")   # change this
fits_folder = Path("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case2_June02/data/aia") 
out_folder  = Path("./aia_boxes_output")
out_folder.mkdir(exist_ok=True, parents=True)
sigma_smooth = 3.0          # gaussian smoothing (pixels)
threshold_mode = "percentile"# "percentile" or "nsigma"
percentile_val = 97          # used if percentile mode
nsigma = 1.0                 # used if nsigma mode (mean + nsigma*std)
min_area = 10000             # minimum area (pixels) to keep a region
pad_pixels = 60              # pad each bbox by this many pixels
max_regions_to_plot = 20     # safety cap

# --- helper functions ---
def compute_threshold(data):
    masked = np.nan_to_num(data, nan=0.0)
    if threshold_mode == "percentile":
        return np.percentile(masked, percentile_val)
    else:
        m = masked.mean()
        s = masked.std()
        return m + nsigma * s

def detect_regions(data, sigma=sigma_smooth, thr=None):
    sm = gaussian_filter(data, sigma=sigma)
    if thr is None:
        thr = compute_threshold(sm)
    mask = sm > thr
    lab = label(mask)
    props = regionprops(lab, intensity_image=sm)
    # filter by area & sort by max intensity
    props = [p for p in props if p.area >= min_area]
    props.sort(key=lambda p: p.max_intensity, reverse=True)
    return props, sm, thr

# --- main processing loop on all fits files in folder ---
rows = []
for fits in sorted(fits_folder.glob("*.131.image_lev1.fits")):
    try:
        m = sunpy.map.Map(fits)
    except Exception as e:
        print(f"Skipping {fits.name}: {e}")
        continue

    data = m.data.astype(float)  # 2D array
    props, smoothed, thr = detect_regions(data)

    print(f"{fits.name}: found {len(props)} regions (after area filter)")
    # limit regions to plot to a sane number
    props_to_use = props[:max_regions_to_plot]

    # Make a figure using pixel coordinates (imshow) so rectangles align simply
    fig, ax = plt.subplots(figsize=(8,8))
    ax.imshow(data, origin='lower', vmin=np.nanpercentile(data, 1), vmax=np.nanpercentile(data, 99.8))
    ax.set_title(f"{fits.name}  AIA 131\nthreshold={thr:.1f}, detected={len(props)}")
    ax.set_xlabel("X (pixels)")
    ax.set_ylabel("Y (pixels)")

    for i, p in enumerate(props_to_use, start=1):
        minr, minc, maxr, maxc = p.bbox  # note: regionprops bbox gives (min_row, min_col, max_row, max_col)
        # swap to (x,y) pixel ordering for plotting (col->x, row->y)
        x0 = minc - pad_pixels
        y0 = minr - pad_pixels
        width  = (maxc - minc) + 2*pad_pixels
        height = (maxr - minr) + 2*pad_pixels

        # Clip to image
        x0 = max(0, x0)
        y0 = max(0, y0)
        mid_ptx=x0+width/2
        mid_pty=y0+height/2
        if x0 + width > data.shape[1]:
            width = data.shape[1] - x0
        if y0 + height > data.shape[0]:
            height = data.shape[0] - y0

        # draw rectangle in pixel coordinates
        rect = Rectangle((x0, y0), width, height, linewidth=1.6, edgecolor='red', facecolor='none')
        ax.plot(mid_ptx,mid_pty,'r+')
        ax.add_patch(rect)
        ax.text(x0+2, y0+10, f"R{i}", color='yellow', fontsize=8, bbox=dict(facecolor='black', alpha=0.4, pad=1))

        # Save a row with pixel bbox and world coords for corners (use SunPy map conversion)
        # Convert four corner pixel coordinates to world (helioprojective) SkyCoord
        # Note: SunPy Map expects pixel coords as (x_pixels, y_pixels) where x=col, y=row
        top_left_world = m.pixel_to_world(x0 * u.pixel, (y0+height) * u.pixel)   # top-left (higher y)
        bottom_right_world = m.pixel_to_world((x0+width) * u.pixel, y0 * u.pixel)
        cen = m.pixel_to_world((x0+width/2) * u.pixel, (y0+height/2) * u.pixel)

        rows.append({
            "file": fits.name,
            "region_id": i,
            "x0_pix": int(x0),
            "y0_pix": int(y0),
            "width_pix": int(width),
            "height_pix": int(height),
            "top_left_hpc_Tx": int(top_left_world.Tx.to(u.arcsec).value),
            "top_left_hpc_Ty": int (top_left_world.Ty.to(u.arcsec).value),
            "bottom_right_hpc_Tx":int(bottom_right_world.Tx.to(u.arcsec).value),
            "bottom_right_hpc_Ty": int(bottom_right_world.Ty.to(u.arcsec).value),
            "Centre_hpc_Tx": int(cen.Tx.to(u.arcsec).value),
            "Centre_hpc_Ty": int(cen.Ty.to(u.arcsec).value),
            "max_intensity": float(p.max_intensity),
            "area_pix": int(p.area)
        })

    # Save figure
    figpath = out_folder / (fits.stem + "_boxes.png")
    plt.tight_layout()
    fig.savefig(figpath, dpi=200)
    plt.close(fig)
    print("Saved:", figpath)

# Save CSV listing all detections
df = pd.DataFrame(rows)
if not df.empty:
    csvpath = "aia131_regions.csv"
    df.to_csv(csvpath, index=False)
    #print("Saved region table to:", csvpath)
else:
    print("No regions saved (empty).")
