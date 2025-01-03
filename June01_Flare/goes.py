import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import time_support

from sunpy import timeseries as ts
from sunpy.net import Fido
from sunpy.net import attrs as a



#results = Fido.search(a.Time("2020-11-20 00:00", "2020-11-21 23:00"), a.Instrument("XRS"), a.Resolution("flx1s"))
#print(results)
results_16 = Fido.search(a.Time("2024-06-01 7:00", "2024-06-01 09:10"), a.Instrument("XRS"), a.goes.SatelliteNumber(16))
print(results_16)
files = Fido.fetch(results_16)
# We use the `concatenate=True` keyword argument in TimeSeries, as
# we have two files and want to create one timeseries from them.
goes_16 = ts.TimeSeries(files, concatenate=True)
#goes_16.peek()

df = goes_16.to_dataframe()
df = df[(df["xrsa_quality"] == 0) & (df["xrsb_quality"] == 0)]
#print(df['xrsb'])


goes_16_ = ts.TimeSeries(df, goes_16.meta, goes_16.units)
fig, ax = plt.subplots()
goes_16_.plot(axes=ax, columns=["xrsb"])
plt.close()

goes_flare = goes_16_.truncate("2024-06-01 7:00", "2024-06-01 09:10")
print(goes_flare.data['xrsb'])
#print(goes_flare.time)
#print(goes_flare["xrsb"])
time_support()
np.savetxt('Goes_data.csv',np.c_[goes_flare.time,goes_flare.data['xrsb']],delimiter=',',fmt='%s')
time_support()
fig, ax = plt.subplots()
ax.plot(goes_flare.time,goes_flare.quantity("xrsb"))
ax.set_ylabel("Flux (Wm$^{-2}$$s^{-1}$)")
ax.set_yscale('log')
fig.autofmt_xdate()
plt.show()


'''
tstart = "2015-06-21 01:00"
tend = "2015-06-21 23:00"
result = Fido.search(a.Time(tstart, tend), a.Instrument("XRS"))#
#print(result)

result_goes15 = Fido.search(a.Time(tstart, tend), a.Instrument("XRS"), a.goes.SatelliteNumber(15), a.Resolution("flx1s"))
#print(result_goes15)
file_goes15 = Fido.fetch(result_goes15)
goes_15 = ts.TimeSeries(file_goes15)
goes_15.peek()

df = goes_15.to_dataframe()
df = df[(df["xrsa_quality"] == 0) & (df["xrsb_quality"] == 0)]
goes_15 = ts.TimeSeries(df, goes_15.meta, goes_15.units)

fig, ax = plt.subplots()
goes_15.plot(axes=ax, columns=["xrsb"])
plt.show()

goes_flare = goes_15.truncate("2015-06-21 09:35", "2015-06-21 10:30")

time_support()
fig, ax = plt.subplots()
ax.plot(goes_flare.time, np.gradient(goes_flare.quantity("xrsb")))
ax.set_ylabel("Flux (Wm$^{-2}$$s^{-1}$)")
fig.autofmt_xdate()
plt.show()

'''