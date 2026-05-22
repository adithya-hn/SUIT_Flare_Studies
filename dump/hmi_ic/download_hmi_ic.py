import os
import drms
client = drms.Client()
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.coordinates  # NOQA
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a


'''
#email_address = os.environ["adithya1@atulbhats.com"]  
client = drms.Client(email='adithya1@atulbhats.com') 
out_dir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Test_moduels/hmi_ic/data'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

export_request = client.export('hmi.M_45s[2024.06.02_TAI/1d@45s]{Magnetogram}', method='url', protocol='fits')
export_request.wait()
export_request.download(out_dir, index=0) 

'''

star_time=Time('2024-06-01T07:08:00', scale='utc', format='isot')
end_time=Time('2024-06-01T07:10:00', scale='utc', format='isot')
out_dir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Test_moduels/hmi_ic/data'

series = a.jsoc.Series('hmi.ic_45s')
query = Fido.search(a.Time(star_time,end_time),series,a.jsoc.Notify("adithya1@atulbhats.com")) #adithya1@atulbhats.com adithyabhattsringeri@gmail.com
print(query)

files = Fido.fetch(query,path=out_dir)


#query = f"hmi.ic_45s[{start_time.strftime('%Y.%m.%d_%H:%M:%S_TAI')}-{end_time.strftime('%Y.%m.%d_%H:%M:%S_TAI')}]"
