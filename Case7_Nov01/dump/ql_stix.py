import numpy as np
import pandas as pd

from pprint import pprint
from datetime import timedelta
from matplotlib import pyplot as plt
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
#from google.colab import data_table
#data_table.enable_dataframe_formatter()

lc = LightCurves.from_sdc(start_utc='2023-07-16T02:00:00', end_utc='2023-07-16T10:00:00', ltc=True)
lc.peek()

print(lc.data.keys())

fig=plt.figure()
plt.plot(lc.time, lc.counts[0,:])

hk=Housekeeping.from_sdc(start_utc='2021-09-06T12:00:00', end_utc='2021-09-06T14:00:00')
pprint(hk.param_names)

elut=EnergyLUT.request('2021-09-03T00:00:00')

l1_obj=PixelData.from_sdc(request_id=2208270251, level='L1')

hdul=l1_obj.hdul
hdul.info()