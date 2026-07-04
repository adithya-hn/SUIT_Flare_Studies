from glob import glob
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import sunpy.map
import numpy as np
import astropy.units as u

channel= 'MgII_2k'  # Channel name for the output file

files = sorted(glob(f"/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case10_Nov13/data/2k_data/*.fits")) 
#/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/aia_fd/aia.lev1_euv_12s.2024-06-01T070237Z.94.image_lev1.fits
times = []
intensities = []

for file in tqdm(files):
    try:
        m = sunpy.map.Map(file)
    except Exception as e:
        print(f"Error processing file {file}: {e}")
        continue
    # Create coordinate grid (Helioprojective in arcsec)
    mask=np.zeros((2048,2048))
    cenX=int(m.meta.get("CRPIX1"))
    cenY=int(m.meta.get("CRPIX2"))
    #print(cenX,cenY)
    yy, xx = np.mgrid[0:2048, 0:2048] #row,col

    # Compute radial distance from disk center
    r = np.sqrt((xx-cenX)**2 + (yy-cenY)**2)  # Tx and Ty are in arcsec

    # 95% solar radius
    r_sun = int(m.meta.get('R_SUN'))
    mask_95 = r < (0.98 * r_sun)
    masked_data = np.where(mask_95, m.data, np.nan)
    total_intensity = np.nansum(masked_data)
    #total_intensity = np.nansum(m.data)
    #m.plot()
    #plt.imshow(masked_data, cmap='gray', norm=plt.Normalize(vmin=np.nanmin(masked_data), vmax=np.nanmax(masked_data)))
    #plt.show()
    

    times.append(m.date.datetime)
    intensities.append(total_intensity)
np.savetxt(f'{channel}_lc.csv',np.c_[times,intensities],delimiter=',',fmt='%s',header='Date,Intensity')
plt.figure(figsize=(10, 5))
plt.plot(times, intensities, marker='o',markersize=0.5)
plt.gcf().autofmt_xdate()
plt.xlabel("Time")
plt.ylabel("Total Intensity (within 95% R$_\odot$)")
plt.title(f"SUIT {2796} Å Light Curve (Disk Only)")
plt.grid(True)
plt.savefig(f'MgII_2k_lc.png', dpi=300, bbox_inches='tight')
plt.show()