

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
start_time='2021-10-09T06:00:00'
end_time=  '2021-10-09T07:00:00'

start_dt = Time(start_time).to_datetime()
end_dt = Time(end_time).to_datetime()
print(f"Start time: {start_dt.date()}, End time: {end_dt}")

#-------------------------
lc = LightCurves.from_sdc(start_utc=start_time, end_utc=end_time, ltc=True)

lc.peek()
plt.savefig('raw_lc_plot.png',dpi=300)
plt.show()

ql_query = Fido.search(a.Time('2021-10-09', '2021-10-09T07:00:00'), a.Instrument.stix, a.stix.DataType.ql)
Path='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/STIX_Flares/case3/data'
files = Fido.fetch(ql_query,path=Path)