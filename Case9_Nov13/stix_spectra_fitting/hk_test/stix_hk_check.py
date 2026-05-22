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
 
from sunpy.time import parse_time
#-------------------------------------

start_tm='2024-11-12T22:22:00'
end_tm  ='2024-11-13T00:37:00'

#---------------------------------------
lc = LightCurves.from_sdc(start_utc=start_tm, end_utc=end_tm, ltc=True)

lc.peek()
plt.savefig('lc.png',dpi=300)
plt.close()
# print(lc.data.keys())
# print(lc.rcr)
fig=plt.figure(figsize=(12,9))
plt.title('Attenuator motion')
plt.plot(lc.time, lc.rcr)
plt.savefig('attenuartor_motion.png',dpi=300)
plt.close()



hk=Housekeeping.from_sdc(start_utc=start_tm, end_utc=end_tm)
# pprint(hk.param_names)

# hk.plot('NIX00078,NIX00079,NIX00080,NIX00081')

# plt.step( hk.data['datetime'], hk.data['raw_values']['NIXD0023'])
# plt.xlabel('UTC')
# plt.ylabel('Mode')
# plt.close()

emph=aux.Ephemeris.from_sdc(start_utc=start_tm, end_utc=end_tm, steps=100)


emph.peek()
plt.close()
# print(parse_time(emph.data['orbit']['utc']))
fig=plt.figure(figsize=(12,9))
plt.title('Solo - earth viewing angle')
plt.plot(parse_time(emph.data['orbit']['utc']).datetime,emph.data['orbit']['earth_sun_solo_angle'])
plt.savefig('solo_earth_angle.png',dpi=300)
plt.show()
