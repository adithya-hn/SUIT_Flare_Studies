import numpy as np
import pywt
import glob
from astropy.io import fits
from tqdm import tqdm
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------
fltr='NB04'
fits_dir = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop/"
files = sorted(glob.glob(fits_dir + f"*{fltr}.fits"))

print('No of files ',len(files))

# choose wavelet and scales
# wavelet_name = "db1"               # Mexican hat (good for transient detection)
# n_scales = 12                       # number of scales
# max_scale = 6                       # high frequency = small scales (< ~3–4)
# min_scale = 1

# scales = np.arange(min_scale, min_scale + n_scales)

# ------------------------------------------------------------
# LOAD DATA CUBE
# ------------------------------------------------------------
cube = []
for f in files:
    data = fits.getdata(f)
    cube.append(data)
cube = np.array(cube, dtype=float) 

# hdu_trans = fits.PrimaryHDU(cube)
# hdu_trans.writeto(f"{fltr}_aligned_cube.fits", overwrite=True)
#cube =fits.getdata ('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/localize_brightenings/NB04_aligned_cube.fits')#np.array(cube, dtype=float)         # (N, Y, X)
ntime, ny, nx = cube.shape

# normalise temporal background for stability
#cube = cube / np.nanmedian(cube, axis=0)
y, x = 300, 300
ts = cube[:, y, x]/np.median(cube[:, y, x])
print(ts.shape)
wavelet = "db2"
max_level = pywt.dwt_max_level(len(ts), pywt.Wavelet(wavelet).dec_len)
print('Max lev ',max_level)

coeffs = pywt.wavedec(ts, wavelet, level=max_level)

def upsample_coeff(c, length):
    up = np.zeros(length)
    step = length // len(c)
    up[::step] = c
    return up

def expand_to_length(c, length):
    factor = int(np.ceil(length / len(c)))
    expanded = np.repeat(c, factor)
    return expanded[:length]

details = coeffs[1:]   # remove approximation
n_levels = len(details)
length = len(ts)

wavelet_map = np.zeros((n_levels, length))

for i, cD in enumerate(details):
    wavelet_map[i] = expand_to_length(cD, ntime)

fig, ax = plt.subplots(figsize=(8, 5))

im = ax.imshow(
    np.abs(wavelet_map),
    aspect="auto",
    origin="lower",
    cmap="inferno"
)

ax.set_xlabel("Time (frame index)")
ax.set_ylabel("DWT level (low → high freq)")
ax.set_title("Daubechies (db4) Wavelet Coefficient Map")

cbar = plt.colorbar(im, ax=ax)
cbar.set_label("|Coefficient|")

plt.tight_layout()
plt.savefig('Wavelet.png',dpi=300)
plt.show()

plt.plot(ts)
plt.savefig('timeseries.png',dpi=300)
plt.show()


"""
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

"""