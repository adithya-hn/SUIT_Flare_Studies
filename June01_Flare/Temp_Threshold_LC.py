import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
import datetime
from pylab import *

# Define the threshold temperature
threshold_temperature = 6_000_000  # 6 million Kelvin

# Directory containing the FITS files
fits_dir = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June01_Flare/Flare_temp/Temp_Fits'
fits_files = glob.glob(f'{fits_dir}/*.fits')

# Lists to store the time and count of pixels above the threshold for each file
times = []
pixels_above_threshold = []

# Process each FITS file
for fits_file in fits_files:
    hdul = fits.open(fits_file)
    data = hdul[0].data  # Assuming the temperature data is in the first extension
    header = hdul[0].header
    date_iso = header.get('DATE-OBS', fits_file)  # Replace 'DATE-OBS' with the appropriate keyword
    date = datetime.datetime.fromisoformat(date_iso)

    num_pixels_above_threshold = np.sum(data > threshold_temperature)
    times.append(date)
    pixels_above_threshold.append(num_pixels_above_threshold)
    hdul.close()

# Convert times to a numpy array for sorting if necessary
times = np.array(times)
pixels_above_threshold = np.array(pixels_above_threshold)

# Sort the data by time if necessary
sorted_indices = np.argsort(times)
times = times[sorted_indices]
pixels_above_threshold = pixels_above_threshold[sorted_indices]

# Plot the light curve
rc('axes', linewidth=1.2)
plt.rcParams["xtick.major.size"] = 10
fig,axs=plt.subplots(1,1, figsize=(10,5))
axs.xaxis.set_tick_params(size=0.5)
axs.yaxis.set_tick_params(size=0.5)
axs.tick_params(axis='both', direction='in', length=6, width=1)
axs.tick_params(which='minor', direction='in', length=3, width=1)
axs.yaxis.set_ticks_position('both')
axs.xaxis.set_ticks_position('both')
axs.minorticks_on()

axs.plot(times, pixels_above_threshold,'k')
plt.xlabel('Time')
plt.ylabel('Number of Pixels Above 6 Million Kelvin')
plt.title('Light Curve: Pixels Above 6 Million Kelvin Over Time')
#plt.grid(True)
plt.show()
