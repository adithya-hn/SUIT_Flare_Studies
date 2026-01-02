import numpy as np
import pywt
import glob
from astropy.io import fits
from tqdm import tqdm
import scipy.ndimage as ndi

# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------
fltr='NB04'
fits_dir = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop_/"
files = sorted(glob.glob(fits_dir + f"*{fltr}.fits"))

wavelet = "haar"           # good for impulsive transients
levels_keep = [1,2,3,4]   # fast timescales
sigma_thresh = 3.0        # significance threshold


# ------------------------------------------------------------
# LOAD DATA CUBE
# ------------------------------------------------------------

cube = []
for f in files:
    data = fits.getdata(f)
    cube.append(data)
cube = np.array(cube, dtype=float)         # (N, Y, X)
cube = cube[:62,:,:]
ntime, ny, nx = cube.shape

# normalise temporal background for stability
# cube = cube / np.nanmedian(cube, axis=0)

transient_cube = np.zeros_like(cube)

def detect_transient_times(ts, wavelet="db4", levels_keep=[1,2], sigma=3.0):
    """
    ts: 1-D intensity time series
    returns: boolean array (ntime,) marking transient times
    """
    coeffs = pywt.wavedec(ts, wavelet)
    ntime = len(ts)

    transient_flag = np.zeros(ntime, dtype=bool)

    for level in levels_keep:
        cD = coeffs[-level]   # cD1 is last element
        med = np.median(cD)
        mad = np.median(np.abs(cD - med))
        thresh = med + sigma * 1.4826 * mad

        idx = np.where(np.abs(cD) > thresh)[0]

        # map DWT indices back to time
        factor = int(np.ceil(ntime / len(cD)))
        for i in idx:
            t0 = i * factor
            t1 = min((i+1)*factor, ntime)
            transient_flag[t0:t1] = True

    return transient_flag

for y in range(ny):
    for x in range(nx): 
        ts = cube[:, y, x]
        if np.all(~np.isfinite(ts)):
            continue

        ts = ts / np.nanmedian(ts)
        flag = detect_transient_times(ts,wavelet=wavelet)
        transient_cube[flag, y, x] = 1

hdu_trans = fits.PrimaryHDU(transient_cube)
hdu_trans.writeto(f"{fltr}_db_transient_cube.fits", overwrite=True)

for t in range(ntime):
    transient_cube[t] = ndi.binary_opening(
        transient_cube[t],
        structure=np.ones((3,3))
    )

hdu_trans = fits.PrimaryHDU(transient_cube)
hdu_trans.writeto(f"{fltr}_db_transient_reg_cube.fits", overwrite=True)

# labels_all = []

# for t in range(ntime):
#     labels, nlab = ndi.label(transient_cube[t])
#     labels_all.append(labels)

# from skimage.measure import regionprops

# event_catalog = []

# for t in range(ntime):
#     labels = labels_all[t]
#     props = regionprops(labels, intensity_image=cube[t])

#     for p in props:
#         event_catalog.append({
#             "time_index": t,
#             "area_px": p.area,
#             "mean_intensity": p.mean_intensity,
#             "centroid_yx": p.centroid
#         })
