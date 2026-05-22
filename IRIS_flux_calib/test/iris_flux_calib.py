import matplotlib.pyplot as plt
import pooch

import astropy.units as u

from irispy.io import read_files
from irispy.utils.spectrograph import radiometric_calibration


raster_filename='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/IRIS_flux_calib/test/iris_l2_20240710_052822_3680201977_raster_t000_r00000.fits'
raster = read_files(raster_filename, memmap=False)
# print(raster)
si_iv_1403 = raster#["Mg II k 2796"][0]
calibrated_si_iv_1403 = radiometric_calibration(si_iv_1403)