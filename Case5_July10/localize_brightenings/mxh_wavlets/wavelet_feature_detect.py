import numpy as np
import pywt
import glob
from astropy.io import fits
from tqdm import tqdm

# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------
fltr='NB04'
fits_dir = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop_/"
files = sorted(glob.glob(fits_dir + f"*{fltr}.fits"))

# choose wavelet and scales
wavelet_name = "gaus1"               # Mexican hat (good for transient detection)
n_scales = 12                       # number of scales
max_scale = 6                       # high frequency = small scales (< ~3–4)
min_scale = 1

scales = np.arange(min_scale, min_scale + n_scales)

# ------------------------------------------------------------
# LOAD DATA CUBE
# ------------------------------------------------------------
cube = []
for f in files:
    data = fits.getdata(f)
    cube.append(data)
cube = np.array(cube, dtype=float)         # (N, Y, X)
ntime, ny, nx = cube.shape

# normalise temporal background for stability
cube = cube / np.nanmedian(cube, axis=0)

# ------------------------------------------------------------
# WAVELET TRANSFORM PER PIXEL
# output: transient_cube (N,Y,X)
# ------------------------------------------------------------
transient_cube = np.zeros_like(cube)

# reshape to (Y*X,N) for vectorised loop
flat = cube.reshape(ntime, -1)      # (N, P)
flat = flat.T                       # (P, N)
n_pix = flat.shape[0]

# container for reconstructed transient component
flat_rec = np.zeros_like(flat)

for idx in tqdm(range(n_pix), desc="Wavelet"):
    ts = flat[idx]

    # skip dead pixels
    if np.all(~np.isfinite(ts)):
        continue

    # Continuous wavelet transform
    coeffs, freqs = pywt.cwt(ts, scales, wavelet_name)

    # coeffs shape: (SCALES, N)
    power = coeffs**2

    # choose high-frequency scales only (small scales)
    # e.g., scale <= max_scale
    hf_mask = scales <= max_scale
    hf_coeffs = coeffs[hf_mask]

    # inverse CWT approximation: average selected coeffs
    # (not exact inverse CWT; this isolates transient component)
    rec = np.mean(hf_coeffs, axis=0)

    flat_rec[idx] = rec

# reshape transient
transient_cube = flat_rec.T.reshape(ntime, ny, nx)

# ------------------------------------------------------------
# OPTIONAL: THRESHOLD & BINARY FEATURE MAPS
# ------------------------------------------------------------
# sigma threshold on transient signal
sigma = 3.0
binary_features = np.zeros_like(transient_cube, dtype=np.uint8)

for i in range(ntime):
    frame = transient_cube[i]
    med = np.nanmedian(frame)
    std = np.nanstd(frame)
    mask = frame > (med + sigma * std)
    binary_features[i][mask] = 1

# ------------------------------------------------------------
# SAVE OUTPUT
# ------------------------------------------------------------
# Example: save transient and binary map as FITS
hdu_trans = fits.PrimaryHDU(transient_cube)
hdu_trans.writeto(f"{fltr}_transient_cube.fits", overwrite=True)

hdu_bin = fits.PrimaryHDU(binary_features)
hdu_bin.writeto(f"{fltr}_transient_features_binary.fits", overwrite=True)

