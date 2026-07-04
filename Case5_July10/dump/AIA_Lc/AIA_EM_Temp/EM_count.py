import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
from datetime import datetime
import pathlib
import os

# Get all FITS files in the directory
EM_files_path='/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/AIA/EM_fits'
fits_files = sorted(glob.glob(EM_files_path+"/*.fits"))  # Adjust the path if needed

pathlib.Path('EM_imgs').mkdir(parents=True, exist_ok=True)

total_counts = []
timestamps = []

for file in fits_files:
    with fits.open(file) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        date_obs = header.get('DATE-OBS', None)
        print(date_obs)
        '''
        plt.imshow(data, cmap='gray',vmin=0,vmax=3e28)
        plt.colorbar()
        plt.title(f'Total EM counts:{date_obs}')
        plt.savefig(f'EM_imgs/{os.path.basename(file)}.png',dpi=200)
        plt.close()'''

        if data is not None:
            total_counts.append(np.nansum(data))  # Sum all pixel values
            
            # Convert DATE-OBS to datetime object
            date_obs = header.get('DATE-OBS', None)
            if date_obs:
                timestamps.append(datetime.strptime(date_obs, "%Y-%m-%dT%H:%M:%S.%f"))
            else:
                timestamps.append(datetime.now())  # Use current time if DATE-OBS is missing

# Sort data based on timestamps
timestamps, total_counts = zip(*sorted(zip(timestamps, total_counts)))

np.savetxt("EM_total_counts.csv", np.c_[timestamps,total_counts],delimiter=',', fmt='%s')

# Plot the total counts over time
plt.figure(figsize=(10, 5))
plt.plot(timestamps, total_counts, marker='o', color='b',markersize=1)
plt.xlabel("Time")
plt.ylabel("Total Counts")
plt.title("Total Counts of AIA Cut-out Images Over Time")
plt.xticks(rotation=45)
plt.grid()

# Save the plot as a PNG file
plt.savefig("AIA_total_counts.png", dpi=300, bbox_inches='tight')
plt.show()

print("Plot saved as 'AIA_total_counts.png'")
