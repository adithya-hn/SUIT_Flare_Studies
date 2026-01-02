import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------
transient_cube_file = "transient_cube.fits"
binary_file = "transient_features_binary.fits"

transient_cube = fits.getdata(transient_cube_file)       # (N,Y,X)
binary_mask = fits.getdata(binary_file)                  # (N,Y,X)

ntime = transient_cube.shape[0]

# ------------------------------------------------------------
# LIGHT CURVE: sum transient counts per frame
#   Option A: sum all transient-enhanced pixels
#   Option B: mean intensity of enhanced pixels
# ------------------------------------------------------------
sum_lc  = np.zeros(ntime)
mean_lc = np.zeros(ntime)
area    = np.zeros(ntime)    # number of enhanced pixels per frame

for i in range(ntime):
    frame = transient_cube[i]
    mask  = binary_mask[i].astype(bool)

    vals = frame[mask]
    if vals.size > 0:
        sum_lc[i]  = np.nansum(vals)
        mean_lc[i] = np.nanmean(vals)
        area[i]    = vals.size
    else:
        sum_lc[i]  = np.nan
        mean_lc[i] = np.nan
        area[i]    = 0

# ------------------------------------------------------------
# SAVE LIGHT CURVES
# ------------------------------------------------------------
hdu = fits.PrimaryHDU(np.vstack([sum_lc, mean_lc, area]))
hdu.header["COMMENT"] = "Row0=sum_lc, Row1=mean_lc, Row2=area"
hdu.writeto("transient_lightcurves.fits", overwrite=True)

# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------
plt.figure()
plt.plot(sum_lc, label="Sum transient counts")
plt.plot(mean_lc, label="Mean transient counts")
plt.plot(area, label="Transient area (pixels)")
plt.xlabel("Frame index")
plt.ylabel("Intensity / Pixel count")
plt.legend()
plt.tight_layout()
plt.savefig("transient_lightcurve.png", dpi=150)
plt.close()
