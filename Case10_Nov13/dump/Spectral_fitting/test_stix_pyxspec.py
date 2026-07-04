from xspec import *

#other packages 
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


file='solo_L1_stix-sci-xray-spec_20241113T110511-20241113T171001_V04_2411132589-50524.fits'
stixf= fits.open(file) 
stixf.info()

'''
ql_query = Fido.search(a.Time('2024-11-13', '2024-11-13T17:50:00'), a.Instrument.stix, a.stix.DataType.ql)

Path='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case10_Nov13/Spectral_fitting/spec_data'
files = Fido.fetch(ql_query,path=Path)
'''

ql_lc = TimeSeries('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case10_Nov13/Spectral_fitting/spec_data/solo_L1_stix-ql-lightcurve_20241113_V02.fits')

ql_lc.plot()
plt.show()