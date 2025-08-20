from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import datetime
from datetime import timedelta


import requests
#from datetime import datetime

# Example date
date = datetime.datetime(2024, 6, 1)


url = f"http://www2.nict.go.jp/aeri/swe/swx/observation/data/yamagawa/fits/yamagawa_{date.strftime('%Y%m%d')}.fits"

response = requests.get(url)
if response.status_code == 200:
    with open(f"yamagawa_{date.strftime('%Y%m%d')}.fits", 'wb') as f:
        f.write(response.content)
    print("FITS file downloaded.")
else:
    print("FITS file not found.")



data=fits.open('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/radio/YAMAGAWA_2024060108I.fits')

spectrum=data[0].data
header=data[0].header

date_str=str(header['DATE-OBS'])+str(header['TIME-OBS'])
#print(date_str)
start_time=datetime.datetime.strptime(date_str,'%Y-%m-%d%H:%M:%S')
num_steps=header['NAXIS1']
step_size=header['CDELT1']
sec_array=np.arange(0,num_steps)
time_array=[start_time+timedelta(seconds=int(ts)) for ts in sec_array]
freq_axis=np.arange(0,8,1)
#print(time_array)
plt.imshow(spectrum,aspect='auto', cmap='inferno', origin='lower',extent=[time_array[0], time_array[-1], freq_axis[0], freq_axis[-1]])
plt.show()