import xarray as xr
import numpy as np
from scipy.io import savemat

ds = xr.open_dataset("sci_xrsf-l2-flx1s_g16_d20241113_v2-2-0.nc")
print(ds)
# Extract arrays
time = ds['time'].values.astype('datetime64[s]')
long_flux = ds['xrsa_flux'].values
short_flux = ds['xrsb_flux'].values
dt=df.datetime
print(long_flux[time>dt])
# Save as MATLAB .sav format readable in IDL
savemat("goes_flux.mat", {"time": time.astype(str), 
                          "long_flux": long_flux, 
                          "short_flux": short_flux})
