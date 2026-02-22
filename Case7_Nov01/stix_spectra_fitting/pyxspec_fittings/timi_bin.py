import numpy as np
from astropy.io import fits

file='stx_spectrum_2410315184.fits'
# hdul = fits.open(file)
# rate_data = hdul[1].data

import numpy as np
from astropy.io import fits
import xspec

# 1. Load Data
hdul = fits.open(file)
data = hdul['RATE'].data
header = hdul['RATE'].header

# 2. Define 1-minute window
# Using the absolute time calculation mentioned in your header comments
mjd_ref = header['MJDREF'] + header['TIMEZERO']
times = data['TIME'] # Seconds from start
mask = (times >= 0) & (times < 60) # First 60 seconds
binned_rows = data[mask]

# 3. Aggregate Rates into Counts
# Total Counts = Sum(Rate * Duration)
rates = binned_rows['RATE']        # Shape (N, 30)
durations = binned_rows['TIMEDEL'] # Shape (N,)
total_counts = np.sum(rates * durations[:, np.newaxis], axis=0)
total_exposure = np.sum(durations)

# 4. Fitting in pyXSPEC
# Note: You need to save total_counts to a temporary .pha file 
# or use a memory-based Spectrum object if your version supports it.

# s = xspec.Spectrum("binned_1min.pha")
# s.response = header['RESPFILE'] # 'stx_srm_2410315184.fits' from your header

# # Define Model: Absorption * (Thermal + Non-thermal)
# m = xspec.Model("phabs * (vapec + bknpower)")

# # Set Initial Parameters
# m.vapec.kT = 1.5           # Initial Temp ~17 MK
# m.bknpower.BreakE = 20.0   # Non-thermal cutoff at 20 keV
# m.bknpower.PhoIndex2 = 4.0 # Non-thermal spectral index (delta)

# # Fit
# xspec.Fit.perform()
# # 1. Define your 1-minute window (60 seconds)
# # Check the start time of the first row
# t_start_global = rate_data['TIME'][0] 
# t_end_target = t_start_global + 60.0

# # 2. Mask the rows that fall into this 1 minute
# mask = (rate_data['TIME'] >= t_start_global) & (rate_data['TIME'] < t_end_target)
# binned_rows = rate_data[mask]

# # 3. Sum the data
# # STIX 'RATE' table usually has a 'COUNTS' or 'DATA' column (dimensions 30 for channels)
# # We sum across the time axis (rows) but keep the 30 energy channels
# summed_counts = np.sum(binned_rows['RATE'], axis=0) 
# total_exposure = np.sum(binned_rows['TIMEDEL']) # Sum of integration times


# stixr= fits.open(file) #not sure why with open() as : syntax doesn't work here?
# stixr.info()
# print(stixr[1].header)
# # Assuming extension 2 is the SPECTRUM extension
# data = stixr[2].data
# header = stixr[2].header

# # 1. Identify rows for your 1-minute window
# # (Check if your file uses 'TIME' or 'TSTART/TSTOP')
# times = data['TIME'] 
# start_time = times[0]
# end_time = start_time + 60.0 # 60 seconds

# mask = (times >= start_time) & (times < end_time)
# binned_data = data[mask]

# # 2. Sum the counts and exposure
# # STIX data is often stored in 'COUNTS' or 'RATE'
# total_counts = np.sum(binned_data['COUNTS'], axis=0)
# total_exposure = np.sum(binned_data['EXPOSURE'])

# # 3. Create a new HDU for XSPEC
# # It is easiest to copy the existing structure and replace the data
# new_header = header.copy()
# # Update header for a single-row (binned) spectrum
# new_header['EXPOSURE'] = total_exposure
# # ... (save this to a new FITS file)