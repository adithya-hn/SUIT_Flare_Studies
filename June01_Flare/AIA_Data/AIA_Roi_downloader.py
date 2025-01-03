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

Path='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June01_Flare/AIA_Data/Data_171'
#peaktime = Time('2012-01-30T21:48:00', scale='utc', format='isot') #2012-01-30T213538
peaktime = Time('2024-06-01T08:45:00', scale='utc', format='isot') 
#end_time=Time('2012-01-30T21:48:00', scale='utc', format='isot')
#star_time=Time('2012-01-30T21:35:38', scale='utc', format='isot')
xc=-467
yc=-313

x1=xc-200
y1=yc-200
x2=xc+200
y2=yc+200

bottom_left = SkyCoord(x1*u.arcsec, y1*u.arcsec, obstime=peaktime, observer="earth", frame="helioprojective")
top_right = SkyCoord(x2*u.arcsec, y2 *u.arcsec, obstime=peaktime, observer="earth", frame="helioprojective")
#start_time = Time('2012-09-24T14:56:03', scale='utc', format='isot')
#bottom_left = SkyCoord(-300*u.arcsec, -175*u.arcsec, obstime=start_time, observer="earth", frame="helioprojective")
#top_right = SkyCoord(50*u.arcsec, 200*u.arcsec, obstime=start_time, observer="earth", frame="helioprojective")

cutout = a.jsoc.Cutout(bottom_left, top_right=top_right, tracking=True)

#jsoc_email = os.environ["adithyabhattsringeri@gmail.com"]
series = a.jsoc.Series('aia.lev1_euv_12s')
query = Fido.search(a.Time(peaktime - 90*u.min, peaktime + 5*u.min),series,a.Wavelength(171*u.angstrom),a.jsoc.Notify("adithya1@atulbhats.com"),cutout,) #a.
print(query)

files = Fido.fetch(query,path=Path)
files.sort()


sequence = sunpy.map.Map(files, sequence=True)

fig = plt.figure()
ax = fig.add_subplot(projection=sequence.maps[0])
ani = sequence.plot(axes=ax, norm=ImageNormalize(vmin=0, vmax=5e3, stretch=SqrtStretch()))
#plt.savefig('movie.gif')
plt.show()