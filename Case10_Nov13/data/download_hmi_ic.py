import os
import drms
client = drms.Client()
import astropy.units as u
from astropy.time import Time
from sunpy.net import Fido
from sunpy.net import attrs as a


star_time=Time('2024-11-13T14:59:00', scale='utc', format='isot')
end_time= Time('2024-11-13T14:59:30', scale='utc', format='isot')
out_dir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case10_Nov13/data/aia'

series = a.jsoc.Series('hmi.ic_45s')
query = Fido.search(a.Time(star_time,end_time),series,a.jsoc.Notify("adithya1@atulbhats.com")) #adithya1@atulbhats.com adithyabhattsringeri@gmail.com
print(query)

files = Fido.fetch(query,path=out_dir)


#query = f"hmi.ic_45s[{start_time.strftime('%Y.%m.%d_%H:%M:%S_TAI')}-{end_time.strftime('%Y.%m.%d_%H:%M:%S_TAI')}]"
