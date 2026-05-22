from astropy.time import Time
from astropy.io import fits
import numpy as np

hdul = fits.open('stx_spectrum_2410315184.fits')

rate = hdul['RATE'].data
energy = hdul['ENEBAND'].data
stix_spec_file='stx_spectrum_2410315184.fits'
tstart='2024-10-31T22:58:00'
tstop ='2024-10-31T22:58:30'


rate = hdul['RATE'].data
header = hdul['RATE'].header
print(header)
mjdref = header['MJDREF']
time_sec = rate['TIME']
tzero= header['timezero']
# ct=hdul['COUNT'].data
# print(ct)
print(hdul[1].header)
print('------------------')
i = 0
print('t del',rate['TIMEDEL'][i])
print('live time',rate['LIVETIME'][i])
print('exposure : ',rate['EXPOSURE'][i])

rate = hdul[1].data['RATE']
exposure = hdul[1].data['EXPOSURE']

counts = rate * exposure[:, None]

spec_counts = counts[mask].sum(axis=0)
total_exposure = exposure[mask].sum()


time_abs = Time(mjdref, format='mjd') #+ time_sec
t_start=Time(tstart)
t_stop=Time(tstop)


mask = (time_abs >= t_start) & (time_abs < t_stop)
print(t_start)
print('-----')
print(time_abs)
print(mask)
print(time_sec)
print(mjdref,tzero)