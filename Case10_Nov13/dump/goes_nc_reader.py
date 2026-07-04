import xarray as xr
import numpy as np
from scipy.io import savemat

# Load your GOES-R NetCDF file
ds = xr.open_dataset('sci_xrsf-l2-flx1s_g16_d20241113_v2-2-0.nc')

# Extract flux and time
long_flux = ds['xrsb_flux'].values  # 1–8 Å
short_flux = ds['xrsa_flux'].values  # 0.5–4 Å
time = ds['time'].values.astype('datetime64[s]').astype(str)

# Save to a .sav file
savemat('goes_flux_for_idl.sav', {
    'time': time,
    'long_flux': long_flux,
    'short_flux': short_flux
})