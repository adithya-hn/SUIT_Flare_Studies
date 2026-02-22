# bad idea just to keep output clean
import warnings
warnings.filterwarnings('ignore')

from sunpy.net import Fido, attrs as a
from stixpy.net import client # Registers stixpy with Fido

from sunpy.timeseries import TimeSeries
from stixpy import timeseries # Registers stixpy TimeSeries with sunpy TimeSeries

from stixdcpy.quicklook import LightCurves
from astropy.table import QTable
from stixdcpy.quicklook import LightCurves
from stixdcpy.energylut import EnergyLUT
from stixdcpy import auxiliary as aux
from stixdcpy.net import FitsQuery as fq
from stixdcpy.net import Request as jreq
from stixdcpy import instrument as inst
from stixdcpy.science import PixelData, Spectrogram, spec_fits_crop, spec_fits_concatenate, fits_time_to_datetime
from stixdcpy.housekeeping import Housekeeping
from astropy.io import fits
from stixdcpy import detector_view as dv
from stixdcpy import spectrogram  as cspec
from stixdcpy.imgspec import ImgSpecArchive as isar
from datetime import datetime
from matplotlib import pyplot as plt

# Start_t="2024-11-13T15:00:30"
# End_t="2024-11-13T17:30:30"

lc = LightCurves.from_sdc(start_utc='2024-11-13T15:00:00', end_utc='2024-11-13T17:30:00', ltc=True)

lc.peek()
plt.show()
# spec_query = Fido.search(a.Time('2021-12-15T09:00', '2021-12-15T12:00:00'), a.Instrument.stix, a.stix.DataProduct.sci_xray_spec,ltc=True)
# spec_query