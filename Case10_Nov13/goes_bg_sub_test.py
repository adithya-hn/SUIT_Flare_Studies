import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
from sunpy.timeseries import TimeSeries
import pandas as pd
from sunpy import timeseries as ts
import astropy.units as u
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
# Create a new TimeSeries with background-subtracted flux
import numpy as np
from sunpy.time import TimeRange
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.metadata import MetaDict
from astropy.table import QTable
import pandas as pd
from sunkit_instruments.goes_xrs import calculate_temperature_em

from plot_goes_lightcurve import plot_goes_lightcurve


# 1. Define the time range
start_time = '2024-11-13 16:30'
end_time = '2024-11-13 17:30'

#result = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"),a.goes.SatelliteNumber(16),a.Resolution("flx1s"))
#print(result)

# 3. Download the data
#downloaded_files = Fido.fetch(result)

# 4. Load as TimeSeries
goes_ts = TimeSeries('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/goes_test/sci_xrsf-l2-flx1s_g16_d20241113_v2-2-0.nc', concatenate=True)
#print(goes_ts)

# Set background
ts_bg = goes_ts.truncate('2024-11-13T16:30:00', '2024-11-13T16:35:00')
goes_flare = goes_ts.truncate("2024-11-13 16:45", "2024-11-13 17:30")
#bg_flux = ts_bg.quantity.mean(axis=0)
#bg_flux = ts_bg.quantity.mean(axis=0)

bg_1_8 = ts_bg.quantity('xrsb').mean()
bg_0_5_4 = ts_bg.quantity('xrsa').mean()

int_a=np.copy(goes_flare.quantity('xrsa'))
int_b=np.copy(goes_flare.quantity('xrsb'))


corrected_ts = goes_flare

print(corrected_ts.quantity('xrsa')[0])
# Subtract background
corrected_ts.quantity('xrsb')[:] -= bg_1_8
corrected_ts.quantity('xrsa')[:] -= bg_0_5_4

print(corrected_ts.quantity('xrsa')[0])

#corrected_flux = goes_ts.quantity - bg_flux
corrected_ts = TimeSeries(data=corrected_ts.data, meta=goes_flare.meta, columns=goes_flare.columns, units=goes_flare.units, source='XRS')
# Get temperature and EM
#corrected_ts = goes_flare.copy()
#corrected_ts.quantity -= bg_flux

print(int_a[0],corrected_ts.quantity('xrsa')[0])

temp_em_ts = calculate_temperature_em(corrected_ts)

fig, ( ax2, ax3) = plt.subplots(2, sharex=True)
#temp_em_ts.plot(columns=["temperature"], axes=ax2)
#temp_em_ts.plot(columns=["emission_measure"], axes=ax3)
ax3.plot(int_a)
ax3.plot(corrected_ts.quantity('xrsa'))
ax3.set_yscale("log")
ax2.grid(True)
plt.savefig('GOES_temp_em1.png',dpi=300)
plt.show()