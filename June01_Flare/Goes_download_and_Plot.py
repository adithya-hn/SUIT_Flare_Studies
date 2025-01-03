'''
from sunpy.net import Fido, attrs as a
from datetime import datetime

# Define the time range for the data
start_time = '2023-10-01T00:00:00'  # Change to your desired start time
end_time = '2023-10-01T23:59:59'    # Change to your desired end time

# Query GOES X-ray data using Fido
results = Fido.search(a.Time(start_time, end_time), a.Instrument('XRS'))

# Download and save the data files
downloaded_files = Fido.fetch(results, path='goes_data/{file}')

print("Files downloaded and saved to 'goes_data/' directory.")
'''
import matplotlib.pyplot as plt
from sunpy.timeseries import TimeSeries
import glob

# Load the GOES data files from the 'goes_data' directory
files = glob.glob('goes_data/*.nc')
goes_ts = TimeSeries(files, source='XRS')

# Plot the 1-8 Ångstrom data
fig, ax = plt.subplots()
goes_ts.plot(columns=['xrsb'], axes=ax, color='blue', label='1-8 Å')

# Customize the plot
ax.set_title('GOES X-ray 1-8 Ångstrom Channel')
ax.set_xlabel('Time')
ax.set_ylabel('Flux (W/m^2)')
plt.legend()
plt.show()
