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
wavelet_name = "mexh"               # Mexican hat (good for transient detection)
n_scales = 12                       # number of scales
max_scale = 6                       # high frequency = small scales (< ~3–4)
min_scale = 1

scales = np.arange(min_scale, min_scale + n_scales)
'''
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
hdu_trans = fits.PrimaryHDU(cube)
hdu_trans.writeto(f"{fltr}_aligned_cube.fits", overwrite=True)'''
cube = fits.getdata(f"{fltr}_aligned_cube.fits")
# ------------------------------------------------------------
# WAVELET TRANSFORM PER PIXEL
# output: transient_cube (N,Y,X)
# ------------------------------------------------------------
transient_cube = np.zeros_like(cube)

'''# reshape to (Y*X,N) for vectorised loop
flat = cube.reshape(ntime, -1)      # (N, P)
flat = flat.T                       # (P, N)
n_pix = flat.shape[0]

# container for reconstructed transient component
flat_rec = np.zeros_like(flat)

for idx in tqdm(range(n_pix), desc="Wavelet"):
    ts = flat[idx]
'''

# choose pixel (y, x)
y, x = 290, 260

ts = cube[:, y, x]    # intensity vs time
ts = ts #/ np.nanmedian(ts)   # normalisation (recommended)
scales = np.arange(1, 32)     # adjust based on cadence

coeffs, freqs = pywt.cwt(ts, scales, wavelet_name)

power = np.abs(coeffs)**2     # wavelet power


import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 5))

im = ax.imshow(
    power,
    aspect="auto",
    origin="lower",
    cmap="inferno",
    extent=[0, len(ts), scales[0], scales[-1]]
)

ax.set_xlabel("Time (frame index)")
ax.set_ylabel("Wavelet scale")
ax.set_title("Wavelet Power Map (Mg II h transient)")

cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Power")

plt.tight_layout()
plt.show()
plt.plot(ts)
plt.show()