from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import datetime
from datetime import timedelta




data=fits.open('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/radio/YAMAGAWA_2024060107I.fits')

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
plt.imshow(spectrum,aspect='auto', cmap='inferno', origin='lower',extent=[time_array[0], time_array[-1], freq_axis[0], freq_axis[-1]],vmax=10)
plt.title('June-01-2024')
plt.savefig('Radio.png',dpi=300)
plt.show()