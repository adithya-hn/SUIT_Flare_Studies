
import sunpy.map

from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord


suit_map=sunpy.map.Map('/home/adithya/Research_Projects/Flare_studies/SUIT_Flares/May31_Flare/SUT_T24_0794_000394_Lev1.0_2024-05-31T20.00.46.964_0973NB03.fits')
suit_pos = get_horizons_coord(-21, suit_map.date)
suit_map.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))

hmi_map = sunpy.map.Map('/home/adithya/Research_Projects/Flare_studies/SUIT_Flares/May31_Flare/SHARP_Data/13697/20240531_200000/hmi.sharp_cea_720s.11297.20240531_200000_TAI.magnetogram.fits')
suit_map.peek()
hmi_map.peek()
out_aia =suit_map.reproject_to(hmi_map.wcs)


out_aia.peek()
out_aia.save('SUIT_cea.fits',overwrite='True')

'''fig,ax=plt.subplots(figsize=(10,5))
ax.imshow(hmi_map.data,cmap='gray')
ax.imshow(out_aia.data,alpha=0.6)
plt.show()'''