

import numpy as np
import pandas as pd

from pprint import pprint
from datetime import timedelta
from matplotlib import pyplot as plt
from astropy.table import QTable
from stixdcpy.quicklook import LightCurves
from stixdcpy.energylut import EnergyLUT
from stixdcpy import auxiliary as aux
from stixdcpy.net import FitsQuery
from stixdcpy.net import Request as jreq
from stixdcpy import instrument as inst
#from stixdcpy.science import ScienceL1 , PixelData, Spectrogram, spec_fits_crop, spec_fits_concatenate, fits_time_to_datetime
from stixdcpy.housekeeping import Housekeeping
from astropy.io import fits
from stixdcpy import detector_view as dv
from stixdcpy import spectrogram  as cspec
from stixdcpy.imgspec import ImgSpecArchive as isar

from stixdcpy.science import  Spectrogram, spec_fits_crop, spec_fits_concatenate

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from datetime import datetime as dt
from sunpy.net import Fido
from sunpy.net import attrs as a
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
from stixpy.net import client
from sunpy.timeseries import TimeSeries
from stixpy import timeseries # Registers stixpy TimeSeries with sunpy TimeSeries


#-------------------------
start_time='2022-03-11T22:00:00'
end_time=  '2022-03-11T23:00:00'

start_dt = Time(start_time).to_datetime()
end_dt = Time(end_time).to_datetime()
print(f"Start time: {start_dt.date()}, End time: {end_dt}")

spec_res=FitsQuery.query(start_dt, end_dt, product_type='xray-spec')
print(len(spec_res))

FitsQuery.fetch(spec_res[:4]) #only fetch first 2 files



hk=Housekeeping.from_sdc(start_utc=start_dt, end_utc=end_dt)
print(hk.param_names)
emph=aux.Ephemeris.from_sdc(start_utc=start_dt, end_utc=end_dt, steps=100)
#hk.plot('NIX00078,NIX00079,NIX00080,NIX00081')