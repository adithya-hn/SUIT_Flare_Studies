
import sunpy.map

from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord


suit_map=sunpy.map.Map('/Analysis/Research_Projects/SUIT_work/Aligned_Fits/SUT_UNP_9999_999999_Lev1.0_2024-06-04T00.00.34.759_0972NB03.fits')
suit_pos = get_horizons_coord(-21, suit_map.date)
suit_map.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
suit_map.save('_SUT_UNP_9999_999999_Lev1.0_2024-06-04T00.00.34.759_0972NB03.fits')