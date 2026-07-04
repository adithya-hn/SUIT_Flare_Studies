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
from scipy.signal import convolve
sns_cl3=sns.color_palette('bright')
# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------
fltr='NB04'
fits_dir = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/alined_img_by_img/"
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
cube = cube[:67,:,:]
hdu_trans = fits.PrimaryHDU(cube)
hdu_trans.writeto(f"{fltr}_cube.fits", overwrite=True)
nt, ny, nx = cube.shape
dt_array = np.array(dt_array)

def haar_kernel(scale):
    return np.concatenate([np.ones(scale), -np.ones(scale)])

def sym_haar(scale):
    return np.r_[-0.5*np.ones(scale),
                 np.ones(scale),
                 -0.5*np.ones(scale)]

def haar_convolution(ts, scales):
    out = {}
    for s in scales:
        h = haar_kernel(s)
        conv = convolve(ts, h, mode='same')
        out[s] = conv
    return out

def sym_haar_convolution(ts, scales):
    out = {}
    for s in scales:
        h = sym_haar(s)
        conv = convolve(ts, h, mode='same')
        out[s] = conv
    return out

# scales = [1, 2, 4, 8, 16]   # ≈ 1–20 min for 85 s cadence
scales = [1, 2, 3, 4, 5,6,7,8]   # ≈ 1–20 min for 85 s cadence

def estimate_noise(ts):
    d1 = np.diff(ts)
    sigma = 1.4826 * np.median(np.abs(d1 - np.median(d1)))
    return sigma

def detect_triggers(ts, scales, nsigma=5):
    sigma0 = estimate_noise(ts)
    triggers = {}

    # convs = haar_convolution(ts, scales)
    convs = sym_haar_convolution(ts, scales)
    for s in scales:
        #sigma_s = np.sqrt(2*s) * sigma0
        sigma_s = np.sqrt(1.5*s) * sigma0
        triggers[s] = convs[s] > nsigma * sigma_s

    return triggers, convs

nt, ny, nx = cube.shape
trigger_cube = {s: np.zeros((nt, ny, nx), dtype=bool) for s in scales}

for y in range(ny):
    for x in range(nx):
        ts = cube[:, y, x]
        if np.any(np.isnan(ts)):
            continue

        triggers, _ = detect_triggers(ts, scales)
        for s in scales:
            trigger_cube[s][:, y, x] = triggers[s]

for s in scales:           
    hdu_trans = fits.PrimaryHDU(trigger_cube[s].astype(int))
    hdu_trans.writeto(f"{fltr}_{s}_sym_haar_cube.fits", overwrite=True)

