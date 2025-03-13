import numpy as np
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord
import astropy.units as u
from sunpy.map import Map
from sunpy.data.sample import AIA_171_IMAGE  # Example Sun image

# Load a sample AIA image to get Sun's parameters (or use your own data)
aia_map = Map(AIA_171_IMAGE)  
R_sun_arcsec = aia_map.rsun_obs.to(u.arcsec).value  # Solar radius in arcsec
sun_center = aia_map.center  # Sun center in arcsec (HPC)

# Example flare locations in HPC (replace with your values)
flare_coords = [
    (-438, -291),  # (X, Y) in arcseconds
    (-263, 298),
    (-210, -285),
    (-148,197),
    (-65,-201),
    (109,119),
    (-413,225),
    (-324,234),
    (-43,-260),
    (152,-206),
    (264,200)
]

# Compute mu for each location
for x_arcsec, y_arcsec in flare_coords:
    r = np.sqrt((x_arcsec - sun_center.Tx.value)**2 + (y_arcsec - sun_center.Ty.value)**2) / R_sun_arcsec
    mu = np.sqrt(1 - r**2) if r <= 1 else 0  # Ensuring valid range
    print(f"Flare at (X={x_arcsec}, Y={y_arcsec}) -> mu = {mu:.3f}")
