import numpy as np
from astropy.io import fits
import glob
from datetime import datetime
from astropy.time import Time
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

files = sorted(glob.glob("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/localize_brightenings/Diff_fits/NB04/*.fits"))

flux = []
time = []

x1,x2 = 100,200
y1,y2 = 450,550

# rect_flux = np.nanmean(data[y1:y2, x1:x2])

for f in files:

    hdu  = fits.open(f)
    data = hdu[0].data

    # spatial integration
    rect_flux = np.nansum(data[y1:y2, x1:x2])
    flux.append(rect_flux)

    # get observation time
    t = hdu[0].header['DATE-OBS']
    time.append(Time(t).datetime)
    fig,ax = plt.subplots()

    ax.imshow(data,origin='lower',vmin=-200,vmax=200)

    rect = Rectangle((x1,y1),
                     x2-x1,
                     y2-y1,
                     edgecolor='red',
                     facecolor='none',
                     linewidth=1.5)

    ax.add_patch(rect)

    plt.savefig(f"frame/frame_{t}.png",dpi=200)
    plt.close()

    hdu.close()

flux = np.array(flux)
import matplotlib.pyplot as plt

plt.plot(time, flux)
plt.ylabel("Integrated ROI Flux")
plt.xlabel("Time")
plt.savefig("roi_lightcurve.png",dpi=300)
plt.show()