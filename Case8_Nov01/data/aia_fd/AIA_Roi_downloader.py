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

Path='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/data/aia_fd'
#peaktime = Time('2012-01-30T21:48:00', scale='utc', format='isot') #2012-01-30T213538
#peaktime = Time('2024-06-02T09:00:00', scale='utc', format='isot') 
start_time=Time('2024-11-01T12:00:00', scale='utc', format='isot')
end_time=Time('2024-11-01T15:00:00', scale='utc', format='isot')


#jsoc_email = os.environ["adithyabhattsringeri@gmail.com", adithya1@atulbhats.com
series = a.jsoc.Series('aia.lev1_euv_12s')
query = Fido.search(a.Time(start_time, end_time),series,a.Wavelength(171*u.angstrom),a.jsoc.Notify("adithya1@atulbhats.com"),a.jsoc.Segment("image")) #a.
print(query)

files = Fido.fetch(query,path=Path)
files.sort()


sequence = sunpy.map.Map(files, sequence=True)

fig = plt.figure()
ax = fig.add_subplot(projection=sequence.maps[0])
ani = sequence.plot(axes=ax, norm=ImageNormalize(vmin=0, vmax=5e3, stretch=SqrtStretch()))
#plt.savefig('movie.gif')
plt.show()