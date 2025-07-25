from glob import glob
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import sunpy.map
import numpy as np
import astropy.units as u

channel= '94'

files = sorted(glob(f"/media/adithya/Adi_Disk3/SUIT_flare_work/June01_Flare/AIA_Data/aia_data/aia.lev1_euv_12s*{channel}.image_lev1.fits")) 
#/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/aia_fd/aia.lev1_euv_12s.2024-06-01T070237Z.94.image_lev1.fits
times = []
intensities = []

for file in tqdm(files):
    m = sunpy.map.Map(file)
    # Create coordinate grid (Helioprojective in arcsec)
    mask=np.zeros((4096,4096))
    cenX=int(m.meta.get("CRPIX1"))
    cenY=int(m.meta.get("CRPIX2"))
    #print(cenX,cenY)
    yy, xx = np.mgrid[0:4096, 0:4096] #row,col

    # Compute radial distance from disk center
    r = np.sqrt((xx-cenX)**2 + (yy-cenY)**2)  # Tx and Ty are in arcsec

    # 95% solar radius
    r_sun = int(m.meta.get('R_SUN'))
    mask_95 = r < (0.95 * r_sun)
    masked_data = np.where(mask_95, m.data, np.nan)
    total_intensity = np.nansum(masked_data)

    times.append(m.date.datetime)
    intensities.append(total_intensity)
np.savetxt(f'{channel}_lc.csv',np.c_[times,intensities],fmt='%s',header='Date,Intensity')
plt.figure(figsize=(10, 5))
plt.plot(times, intensities, marker='o',markersize=0.5)
plt.gcf().autofmt_xdate()
plt.xlabel("Time")
plt.ylabel("Total Intensity (within 95% R$_\odot$)")
plt.title(f"AIA {channel} Å Light Curve (Disk Only)")
plt.grid(True)
plt.show()