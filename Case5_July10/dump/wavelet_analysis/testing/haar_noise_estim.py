import numpy as np
import pywt
import glob
from astropy.io import fits
from tqdm import tqdm
import scipy.ndimage as ndi
from astropy.time import Time
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------
fltr='NB04'
fits_dir = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/wavelet_analysis/derotated_common_crop/"
files = sorted(glob.glob(fits_dir + f"*{fltr}.fits"))

# wavelet = "haar"           # good for impulsive transients
# levels_keep = [1,2,3,4,5,6]   # fast timescales
# sigma_thresh = 1.0        # significance threshold


# ------------------------------------------------------------
# LOAD DATA CUBE
# ------------------------------------------------------------

cube = []
dt_array= []
for f in files:
    data = fits.getdata(f)
    header=fits.getheader(f)
    cube.append(data)
    date_obs=header['DATE-OBS']
    t = Time(date_obs, format='isot', scale='utc')
    dt_array.append(t.to_datetime())
cube = np.array(cube, dtype=float)         # (N, Y, X)
cube = cube[:,:,:]

dt_array = np.array(dt_array)
#print(dt_array)
time_sec = dt_array.astype("datetime64[s]").astype(np.int64)

dt = np.diff(time_sec)

gap_idx = np.where(dt > 90)[0]##
segments = []

start = 0
for g in gap_idx:
    end = g + 1
    segments.append((start, end))
    start = end

segments.append((start, len(time_sec)))
valid_indices = []

for s, e in segments:
    if e - s >= 3:               # need at least 3 points
        valid_indices.extend(range(s+1, e-1))

valid_indices = np.array(valid_indices)
# Haar = first difference
D = cube[1:] - cube[:-1]        # shape (nt-1, ny, nx)

# keep only valid differences
valid_D = D[valid_indices - 1]  # shift because diff reduces length by 1
noise_samples = valid_D.reshape(-1)
sigma = 1.4826 * np.median(np.abs(noise_samples - np.median(noise_samples)))
std= np.std(noise_samples)
mn=np.mean(noise_samples)
# x = noise_samples[np.isfinite(noise_samples)]
# bin_width = 10.0

# xmin = np.floor(x.min() / bin_width) * bin_width  ``
# xmax = np.ceil(x.max() / bin_width) * bin_width

# bins = np.arange(xmin, xmax + bin_width, bin_width)
# plt.figure(figsize=(6,4))
# plt.hist(x, bins=bins, histtype='step', linewidth=1.8)
# plt.xlabel("Haar first difference (counts)")
# plt.ylabel("Number of samples")
# plt.title("Noise distribution (bin width = 10 counts)")
# plt.tight_layout()
# plt.savefig('Noise_distrib.png',dpi=300)
# plt.show()
print('mean: ',mn,'std: ',std,'MAD sigma: ', sigma)
sigma_map = 1.4826 * np.median(
    np.abs(valid_D - np.median(valid_D, axis=0)),
    axis=0
)  # shape (ny, nx)
trigger = np.abs(D) > 4.0 * sigma_map
triggers=trigger.astype(int)
ntime, ny, nx = cube.shape
plt.imshow(sigma_map)
plt.colorbar()
plt.savefig('noisemap.png',dpi=300)
plt.show()

plt.plot(valid_indices)
plt.savefig('indicies.png',dpi=300)
plt.show()

hdu_trans = fits.PrimaryHDU(valid_D)
hdu_trans.writeto(f"{fltr}_full_first_diff_cube.fits", overwrite=True)

hdu_trans = fits.PrimaryHDU(triggers)
hdu_trans.writeto(f"{fltr}_full_trigger_cube.fits", overwrite=True)


