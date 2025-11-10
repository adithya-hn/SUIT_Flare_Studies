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

lc = LightCurves.from_sdc(start_utc='2024-11-01T15:00:00', end_utc='2024-11-01T17:15:00', ltc=True)

lc.peek()
plt.title('STIX_ql_lc')
plt.savefig('ql_lc.png',dpi=300)
plt.close()
# print(lc.data.keys())

fig=plt.figure(figsize=(12,6))
plt.plot(lc.time, lc.counts[0,:],label='4-10 keV')
plt.plot(lc.time, lc.counts[1,:],label='10-15 keV')
plt.plot(lc.time, lc.counts[2,:],label='15-25 keV')
plt.plot(lc.time, lc.counts[3,:],label='25-84 keV')
plt.ylim(1e1,1e5)
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('ql.png',dpi=300)
plt.show()

# hk=Housekeeping.from_sdc(start_utc='2021-09-06T12:00:00', end_utc='2021-09-06T14:00:00')
# pprint(hk.param_names)

# elut=EnergyLUT.request('2024-10-08T00:00:00')

# l1_obj=PixelData.from_sdc(request_id=2208270251, level='L1')

# hdul=l1_obj.hdul
# hdul.info()