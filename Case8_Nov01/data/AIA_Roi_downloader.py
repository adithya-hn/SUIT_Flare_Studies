import os

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.coordinates  # NOQA
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
import drms
client = drms.Client()

Path='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/data/aia_dem_in'

start_time=Time('2024-11-01T11:00:00', scale='utc', format='isot')
end_time=Time('2024-11-01T15:00:00', scale='utc', format='isot')

x1=470
y1=50
x2=130
y2=400

bottom_left = SkyCoord(x1*u.arcsec, y1*u.arcsec, obstime=start_time, observer="earth", frame="helioprojective")
top_right = SkyCoord(x2*u.arcsec, y2 *u.arcsec, obstime=start_time, observer="earth", frame="helioprojective")

#start_time = Time('2012-09-24T14:56:03', scale='utc', format='isot')
#bottom_left = SkyCoord(-300*u.arcsec, -175*u.arcsec, obstime=start_time, observer="earth", frame="helioprojective")
#top_right = SkyCoord(50*u.arcsec, 200*u.arcsec, obstime=start_time, observer="earth", frame="helioprojective")

cutout = a.jsoc.Cutout(bottom_left, top_right=top_right, tracking=True)

#jsoc_email = os.environ["adithyabhattsringeri@gmail.com"]

series = a.jsoc.Series('aia.lev1_euv_12s')
wavelengths = [94,131,171,193,211,335]*u.angstrom
query = Fido.search(a.Time(start_time, end_time),series,a.AttrOr([a.Wavelength(wav) for wav in wavelengths]),a.jsoc.Notify("adithya1@atulbhats.com"),cutout,) #a.
print(query)

files = Fido.fetch(query,path=Path)



