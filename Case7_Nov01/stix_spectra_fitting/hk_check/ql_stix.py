import numpy as np
import pandas as pd

from pprint import pprint
from datetime import timedelta
from matplotlib import pyplot as plt
from astropy.table import QTable
import stixdcpy
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
start_tm='2024-11-01T00:16:00'
end_tm  ='2024-11-01T02:31:00'


lc = LightCurves.from_sdc(start_utc=start_tm, end_utc=end_tm, ltc=True)

lc.peek()
plt.title('STIX_ql_lc')
plt.savefig('ql_lc.png',dpi=300)
plt.show()

