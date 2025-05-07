import numpy as np
import matplotlib.pyplot as plt
import sunpy.map
from skimage import measure
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib.path import Path
from astropy.io import fits
from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord

# Load two maps
iris_img = fits.open("../data/raw/flare_data/iris_l2_20240710_152851_3620108477_SJI_2796_t000.fits")
suit_map = sunpy.map.Map("/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/july10_11ut_to_18ut/SUT_T24_0956_000465_Lev1.0_2024-07-10T15.35.28.135_0983NB03.fits")

#SUIT plate scale correction
suit_pos = get_horizons_coord(-21, suit_map.date)
suit_map.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
        
sji_imgs=iris_img[0].data
iris_map=sunpy.map.Map(sji_imgs[21], iris_img[0].header) #SJI flare map

# QS region for threshold values
qs_coord=SkyCoord(Tx=(-10, 90) * u.arcsec, Ty=(-55, -110) * u.arcsec, frame=suit_map.coordinate_frame)
qs_map=suit_map.submap(qs_coord)
qs_data = qs_map.data * 1000 / qs_map.meta.get('CMD_EXPT')
print(np.median(qs_data),np.mean(qs_data),np.std(qs_data))
Thresh_val=np.median(qs_data)*3
print('Threshold: ',Thresh_val)

normalized_data = suit_map.data * 1000 / suit_map.meta.get('CMD_EXPT')
threshold = Thresh_val  # adjust as needed
contours_pix = measure.find_contours(normalized_data, threshold)
largest_contour = max(contours_pix, key=len) #finding the largest contour

# Converting to Helip-projective coordinates
hpc_coords2 = suit_map.pixel_to_world(largest_contour[:, 1]*u.pixel, largest_contour[:, 0]*u.pixel)
hpc_coords2= hpc_coords2.transform_to(iris_map.coordinate_frame) # not much chage

fig = plt.figure()
ax = plt.subplot(projection=suit_map)
#ax2 = fig.add_subplot(projection=iris_map)

suit_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent,alpha=0.8)
suit_map.draw_quadrangle(qs_coord,axes=ax,edgecolor="blue",linestyle="-",linewidth=1,label='QS region',alpha=0.5)
plt.colorbar()
ax.plot_coord(hpc_coords2, color='red', linewidth=1,label='Flare contour')
iris_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent,autoalign=True,alpha=0.7)

plt.title("Flare contour of SUIT overlaid on IRIS")
plt.colorbar()
#plt.legend() #taking so much time
plt.savefig("../results/suit_iris_contour.png", dpi=300)
plt.show()
