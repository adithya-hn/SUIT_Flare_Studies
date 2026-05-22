from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord
import astropy.units as u

obstime="2024-06-01T07:29:00"
# Define L1 satellite position in Heliocentric coordinates
l1_position = SkyCoord(
    x=1.5e6 * u.km,  # Distance from Earth towards the Sun
    y=0 * u.km,
    z=0 * u.km,
    frame=frames.Heliocentric,
    obstime=obstime
)
# Example: Convert from heliographic to helioprojective
hg_coord = SkyCoord(lon=21*u.deg, lat=18*u.deg, frame=frames.HeliographicStonyhurst, observer=l1_position,obstime=obstime)


hp_coord = hg_coord.transform_to(frames.Helioprojective(observer="earth"))
print(hp_coord)