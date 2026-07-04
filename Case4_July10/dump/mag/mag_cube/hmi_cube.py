"""
@Author      : Adithya H N
@Created On  : 2026-02-15
@Last Updated: 2026-02-15
@Project     : Project Name
@Version     : 1.0

@Description
-----------
Brief description
"""


import matplotlib
matplotlib.use('Agg')
import sunpy.map
from sunpy.physics.differential_rotation import solar_rotate_coordinate
from sunpy.physics.differential_rotation import differential_rotate
from astropy.time import Time
import glob
from sunpy.coordinates import Helioprojective, SphericalScreen, propagate_with_solar_surface
import warnings
import logging
from reproject import reproject_interp
import numpy as np
from astropy.io import fits


warnings.simplefilter('ignore')
log_ = logging.getLogger('sunpy')
log_.setLevel('WARNING')
files = sorted(glob.glob("/media/adithya/Adi_disk4/SUIT_flare_work/case4_jul10/data/HMI/HMI_cutouts/*.fits"))

maps = [sunpy.map.Map(f) for f in files]

# Choose reference time (first frame usually)
ref_time = maps[0].date
# derotated_maps = []

# for m in maps:
#     m_derot = differential_rotate(m, time=ref_time)
#     print((m_derot.data).shape,'----')
#     derotated_maps.append(m_derot)


# data_cube = np.array([m.data for m in derotated_maps])

# hdu = fits.PrimaryHDU(data_cube)
# hdu.writeto("HMI_derot_cube.fits", overwrite=True)

roi_maps=maps

ref_map = roi_maps[0]
ny, nx = ref_map.data.shape
# shapes = [m.data.shape for m in roi_maps]
# min_y = min(s[0] for s in shapes)
# min_x = min(s[1] for s in shapes)

# cube = np.array([m.data[:min_y, :min_x] for m in roi_maps])
cube = []

for m in roi_maps:
    data_reproj, _ = reproject_interp(m, ref_map.wcs, shape_out=(ny, nx))
    cube.append(data_reproj)
    print((data_reproj).shape,'----')

cube = np.array(cube).astype("float32")
hdu = fits.PrimaryHDU(cube)
hdu.writeto("HMI_derot_cube.fits", overwrite=True)