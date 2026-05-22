import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
from sunpy import timeseries as ts
from sunpy.data.sample import GOES_XRS_TIMESERIES
import numpy as np
from sunkit_instruments import goes_xrs
import matplotlib.dates as mdates
import datetime
from sunpy.timeseries.sources.goes import XRSTimeSeries
import astropy.units as u

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
# set_pub_style()


# 1. Define the time range
start_time = '2024-05-11 01:10'
end_time =   '2024-05-11 01:40'

# ost="2024-11-01 14:10"
# est="2024-11-01 14:35"

result = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"),a.goes.SatelliteNumber(16),a.Resolution("avg1m"))

# 3. Download the data
downloaded_files = Fido.fetch(result)

# 4. Load as TimeSeries
#goes_ts = TimeSeries(downloaded_files, concatenate=True)

goes_ts = ts.TimeSeries(downloaded_files,source='XRS')
goes_df = goes_ts.to_dataframe()
#print(goes_ts.quantity("xrsa").unit)
df_rd = goes_df.diff(periods=3)
# #running_diff_3min[running_diff_3min <= 0] = float("nan")
# running_diff_3min["xrsa"] = running_diff_3min["xrsa"].values * u.W / u.m**2
# running_diff_3min["xrsb"] = running_diff_3min["xrsb"].values * u.W / u.m**2
df_rd["xrsa"] = df_rd["xrsa"].values
df_rd["xrsb"] = df_rd["xrsb"].values
# df_rd["xrsa"] = df_rd["xrsa"].values * u.W / u.m**2
# df_rd["xrsb"] = df_rd["xrsb"].values * u.W / u.m**2


diff_ts = XRSTimeSeries(df_rd, meta=goes_ts.meta)
diff_ts.units = goes_ts.units
#diff_ts.meta = goes_ts.meta
# print(diff_ts.units)
# print(diff_ts)
print('-----------------------------')

goes_flare = goes_ts.truncate(start_time, end_time)
goes_temp_em =goes_xrs.calculate_temperature_em(goes_flare) #goes_xrs.calculate_temperature_em(goes_flare)

goes_diff_flare = diff_ts.truncate(start_time, end_time)
goes_diff_temp_em =goes_xrs.calculate_temperature_em(goes_diff_flare)
goes_flare_df = goes_diff_flare.to_dataframe()
neg_mask = (goes_flare_df["xrsa"].values <= 0) | (goes_flare_df["xrsb"].values <= 0)
df1 = goes_temp_em.to_dataframe()
df2 = goes_diff_temp_em.to_dataframe()
# print(len(df1['temperature']), len(df2['temperature']))
msk2= (df2["temperature"].values > 100)
df2.loc[neg_mask, "temperature"] = 0.0
df2.loc[msk2, "temperature"] = 0.0
df2.loc[neg_mask, "emission_measure"] = 0.0


# ts_modified = XRSTimeSeries(df2,meta=goes_ts.meta,units=goes_ts.units) 
# df2 = ts_modified.to_dataframe()
fig, ( ax2, ax3) = plt.subplots(2, sharex=True)
fig.suptitle("X5.8 Class Flare on 2024-05-11 01:10:00 UTC")

#ax2.plot(goes_flare_df.index, goes_flare_df["xrsa"], color="gray", label="GOES XRS-A")
ax2.plot(df1.index, df1["temperature"], color="k",markersize=2,marker='o', label="Temperature")
ax2.plot(df2.index, df2["temperature"], color="r",markersize=2,marker='o', label="Delta Temperature")
ax2.axhline(5, color="k", ls="--", lw=1)
ax2.axvline(datetime.datetime(2024, 5, 11, 1, 13), color="b", ls="-", lw=1)
ax2.set_ylabel("Temperature (MK)")
ax2.legend(loc="upper right")

#ax2 = ax2.twinx()
ax3.plot(df1.index, df1["emission_measure"]/1e50,markersize=2,marker='o', color="k",label="Emission Measure")
ax3.plot(df2.index, df2["emission_measure"]/1e50,markersize=2,marker='o', color="r",label="Delta Emission Measure")

ax3.set_ylabel("Emission Measure ($10^{50}$ cm$^{-3}$)")
ax3.axhline(5e-3, color="k", ls="--", lw=1)
ax3.axvline(datetime.datetime(2024, 5, 11, 1, 13), color="b", ls="-", lw=1)
ax3.legend(loc="upper right")
ax3.set_yscale("log")
ax3.set_ylim(1e-3,1e1)
ax2.grid(True)
ax3.grid(True)

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig('GOES_temp_em1.png',dpi=300)
plt.show()

df1.index.name = "date_time"
df1.to_csv('goes_temp_em.csv')
df2.index.name = "date_time"
df2.to_csv('goes_diff_temp_em.csv')