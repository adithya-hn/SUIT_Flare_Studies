import numpy as np
import pywt
import glob
from astropy.io import fits
from tqdm import tqdm
import scipy.ndimage as ndi
from astropy.time import Time
import matplotlib.pyplot as plt
from skimage import filters, measure
from skimage.measure import label, regionprops
from skimage.morphology import disk, closing, opening,remove_small_objects
import seaborn as sns
sns_cl3=sns.color_palette('bright')
# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------
fltr='NB04'
fits_dir = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop_/"
files = sorted(glob.glob(fits_dir + f"*{fltr}.fits"))

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
nt, ny, nx = cube.shape
dt_array = np.array(dt_array)
#print(dt_array)
time_sec = dt_array.astype("datetime64[s]").astype(np.int64)

dt = np.diff(time_sec)

gap_idx = np.where(dt > 111)[0]##
segments = []
print('data_gap indices',gap_idx)
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
scales = np.array([1, 2, 3, 4, 5, 6,7,8])
responses = {}

for s in scales:
    Ds = cube[s:] - cube[:-s]   # shape (nt-s, ny, nx)
    responses[s] = Ds

sigma_s = {}

for s in scales:
    Ds_valid = responses[s][valid_indices[:-s]]
    noise = Ds_valid.reshape(-1)
    sigma_s[s] = 1.4826 * np.median(np.abs(noise - np.median(noise)))

significance = {}

for s in scales:
    significance[s] = responses[s] / sigma_s[s]

k = 6.0 

triggers = {
    s: np.abs(significance[s]) > k
    for s in scales
}




scale_triggers = {}

for s in scales:
    full = np.zeros((nt, ny, nx), dtype=bool)
    trig = triggers[s]                  # shape (nt-s, ny, nx)

    offset = s # // 2                     # centered alignment
    full[offset:offset + trig.shape[0]] = trig

    scale_triggers[s] = full
    hdu_trans = fits.PrimaryHDU(scale_triggers[s].astype(int));hdu_trans.writeto(f"{fltr}_t{s}_cube.fits", overwrite=True)
    
dominant_map = np.zeros((nt, ny, nx), dtype=int)

for s in scales:
    dominant_map[scale_triggers[s]] = s

# # any scale triggers
# transient_mask = np.zeros_like(dominant_map, dtype=bool)

# temporal = transient_mask.copy()
# temporal[1:-1] &= (transient_mask[:-2] | transient_mask[2:])

# scale_persistent = stack >= 2
# final_transient = temporal | scale_persistent

# from scipy.ndimage import label



# for t_idx in range(78):
#     labels, _ = label(final_transient[t_idx])
#     sizes = np.bincount(labels.ravel())

#     min_area = 4  # pixels
#     valid = sizes >= min_area
#     valid[0] = False

#     final_transient[t_idx] = valid[labels]
#     baseline = np.median(cube[t_idx-2:t_idx], axis=0)
#     peak_excess = cube[t_idx] - baseline
#     energy = np.sum(peak_excess[event_mask])