import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
from sunpy import timeseries as ts
from sunpy.data.sample import GOES_XRS_TIMESERIES
import numpy as np
from sunkit_instruments import goes_xrs


# 1. Define the time range
start_time = '2024-11-13 16:50'
end_time = '2024-11-13 17:30'

# 2. Search GOES data
#result = Fido.search(a.Time(start_time, end_time), a.Instrument('goes'))
#results = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"), a.Resolution("flx1s"))
#result_goes15 = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"), a.goes.SatelliteNumber(15), a.Resolution("flx1s"))
result = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"),a.goes.SatelliteNumber(16),a.Resolution("flx1s"))

# 3. Download the data
downloaded_files = Fido.fetch(result)

# 4. Load as TimeSeries
#goes_ts = TimeSeries(downloaded_files, concatenate=True)


goes_ts = ts.TimeSeries(downloaded_files)
#goes_ts.plot()

goes_flare = goes_ts.truncate("2024-11-13 16:45", "2024-11-13 17:30")
#background_tr = ts.TimeRange("2024-11-13 16:45", "2024-11-13 16:47")

# Download GOES XRS data for background interval
background_lc =goes_ts.truncate("2024-11-13 17:05",  "2024-11-13 17:10")
#goes_flare.update_background_times("2024-11-13 16:45", "2024-11-13 16:47")
# Estimate background flux as median of background interval for each channel

bg_1_8 = background_lc.quantity('xrsb').mean()
bg_0_5_4 = background_lc.quantity('xrsa').mean()

int_a=np.copy(goes_flare.quantity('xrsa'))
int_b=np.copy(goes_flare.quantity('xrsb'))

lc1 = goes_flare.quantity('xrsb') - bg_1_8
lc2 = goes_flare.quantity('xrsa') - bg_0_5_4
print(lc1,lc2)
goes_flare.quantity('xrsb').value=lc1
goes_flare.quantity('xrsa').value=lc2
goes_temp_em =goes_xrs.calculate_temperature_em(goes_flare) #goes_xrs.calculate_temperature_em(goes_flare)
print(goes_temp_em.columns)
data=goes_temp_em.to_dataframe()
print(data['temperature'])

start_time = "2024-11-13T07:45:00.000"  # ISOT start time
end_time   = "2024-11-13T08:20:00.000" 

plt.plot(int_b-goes_flare.quantity('xrsb'))
plt.plot(int_a-goes_flare.quantity('xrsa'))
plt.show()

#print(data['time'])

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
goes_flare.plot(axes=ax1)
goes_temp_em.plot(columns=["temperature"], axes=ax2)
goes_temp_em.plot(columns=["emission_measure"], axes=ax3)
ax3.set_yscale("log")
plt.savefig('GOES_temp_em.png',dpi=300)
plt.show()

fig, ( ax2, ax3) = plt.subplots(2, sharex=True)
goes_temp_em.plot(columns=["temperature"], axes=ax2)
goes_temp_em.plot(columns=["emission_measure"], axes=ax3)
ax3.set_yscale("log")
ax2.grid(True)
plt.savefig('GOES_temp_em1.png',dpi=300)
plt.show()