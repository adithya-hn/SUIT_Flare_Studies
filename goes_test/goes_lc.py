import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
from sunpy.timeseries import TimeSeries
import pandas as pd

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plot_goes_lightcurve import plot_goes_lightcurve

# 1. Define the time range
start_time = '2015-11-12 22:30'
end_time = '2015-11-13 00:30'

# 2. Search GOES data
#result = Fido.search(a.Time(start_time, end_time), a.Instrument('goes'))
#results = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"), a.Resolution("flx1s"))
#result_goes15 = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"), a.goes.SatelliteNumber(15), a.Resolution("flx1s"))
result = Fido.search(a.Time(start_time, end_time), a.Instrument("XRS"),a.goes.SatelliteNumber(15),a.Resolution("flx1s"))

print(result)

# 3. Download the data
downloaded_files = Fido.fetch(result,path='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/goes_test')

# 4. Load as TimeSeries
goes_ts = TimeSeries(downloaded_files, concatenate=True)

# 5. Convert to pandas DataFrame
df = goes_ts.to_dataframe()
df.index.name = "Time"  # Set index name
#df =df[(df["Time"] > start_time )& (df["Time"]<end_time)] 
df = df[(df["xrsa_quality"] == 0) & (df["xrsb_quality"] == 0)]

# 4. Filter data between start and end time
df_filtered = df.loc[start_time:end_time]

# 6. Save to CSV
csv_filename = "goes_xray_lightcurve_20241113_1.csv"
df_filtered.to_csv(csv_filename)

print(f"GOES light curve saved to {csv_filename}")


csv_file = csv_filename 
plot_goes_lightcurve(csv_file, save_plot=True, output_file="goes_lc_plot.png")
