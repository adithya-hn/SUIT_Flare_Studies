import numpy as np
from astropy.io import fits
import glob
from datetime import datetime
import pywt
from skimage.measure import label, regionprops
import matplotlib.pyplot as plt
import pandas as pd



def load_fits_sequence(folder, pattern="*NB04.fits"):
    files = sorted(glob.glob(f"{folder}/{pattern}"))
    data = []
    times = []

    for f in files:
        with fits.open(f) as hdul:
            data.append(hdul[0].data.astype(float))
            hdr = hdul[0].header
            times.append(datetime.fromisoformat(hdr['DATE-OBS']))

    return np.array(data), np.array(times)

data_cube, time_array = load_fits_sequence("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop")
# shape = (nt, ny, nx)

from scipy.signal import detrend

def detrend_cube(data_cube):
    nt, ny, nx = data_cube.shape
    out = np.empty_like(data_cube)

    for y in range(ny):
        for x in range(nx):
            out[:, y, x] = detrend(data_cube[:, y, x], type='linear')

    return out

data_cube_dt = detrend_cube(data_cube)

def wavelet_power_cube(data_cube, scales, wavelet='gaus2'):
    """
    Computes wavelet power for each pixel along time axis
    """
    nt, ny, nx = data_cube.shape
    power_cube = np.zeros((len(scales), nt, ny, nx))

    for y in range(ny):
        for x in range(nx):
            coeffs, _ = pywt.cwt(
                data_cube[:, y, x],
                scales,
                wavelet
            )
            power_cube[:, :, y, x] = np.abs(coeffs) ** 2

    return power_cube

scales = np.arange(2, 30)   # adjust for cadence
power_cube = wavelet_power_cube(data_cube_dt, scales)

power_sum = np.sum(power_cube, axis=0)
# shape = (nt, ny, nx)
def mad_threshold(data, nsigma=5):
    med = np.nanmedian(data)
    mad = np.nanmedian(np.abs(data - med))
    return med + nsigma * 1.4826 * mad

threshold = mad_threshold(power_sum, nsigma=5)

detection_mask = power_sum > threshold
# shape = (nt, ny, nx)

def detect_regions(mask_2d, min_area=5):
    labels = label(mask_2d)
    regions = [r for r in regionprops(labels) if r.area >= min_area]
    return regions

def save_contour_images(data_cube, detection_mask, outdir="contours"):
    import os
    os.makedirs(outdir, exist_ok=True)

    nt = data_cube.shape[0]

    for t in range(nt):
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.imshow(data_cube[t], origin='lower', cmap='gray')

        regions = detect_regions(detection_mask[t])
        for r in regions:
            y0, x0 = r.centroid
            ax.plot(x0, y0, 'ro', markersize=3)

        ax.set_title(f"Frame {t}")
        ax.axis('off')
        plt.savefig(f"{outdir}/frame_{t:04d}.png", dpi=150)
        plt.close()

save_contour_images(data_cube, detection_mask)

def extract_lightcurves(data_cube, detection_mask, time_array):
    nt, ny, nx = data_cube.shape
    curves = []

    for y in range(ny):
        for x in range(nx):
            if detection_mask[:, y, x].any():
                lc = data_cube[:, y, x]
                curves.append(pd.DataFrame({
                    "time": time_array,
                    "intensity": lc,
                    "y": y,
                    "x": x
                }))

    return pd.concat(curves, ignore_index=True)


#------------------

def plot_wavelet_map(signal, scales, wavelet='gaus2'):
    coeffs, freqs = pywt.cwt(signal, scales, wavelet)
    power = np.abs(coeffs)**2

    plt.figure(figsize=(6,4))
    plt.imshow(power, aspect='auto', origin='lower',
               extent=[0, len(signal), scales[0], scales[-1]])
    plt.colorbar(label="Power")
    plt.xlabel("Time index")
    plt.ylabel("Scale")
    plt.title("Wavelet Power Spectrum")
    plt.show()
plot_wavelet_map(data_cube_dt[:, 50, 50], scales)
