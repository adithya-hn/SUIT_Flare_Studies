import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS

from scipy.constants import c
from scipy.interpolate import interp1d

# Set up some default matplotlib options
#plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['image.origin'] = 'lower'
plt.rcParams['image.cmap'] = 'viridis'



iris_file = fits.open("iris_l2_20240710_152851_3620108477_raster_t000_r00000.fits",memmap=True, do_not_scale_image_data=True)
hd = iris_file[0].header
#print('Primary header:',hd)
print('Window. Name      : wave start - wave end\n')
for i in range(hd['NWIN']):
    win = str(i + 1)
    print('{0}. {1:15}: {2:.2f} - {3:.2f} Å' ''.format(win, hd['TDESC' + win], hd['TWMIN' + win], hd['TWMAX' + win]))

data = iris_file[8].data
print(data.shape)
wcs = WCS(iris_file[8].header)
m_to_nm = 1e9  # convert wavelength to nm
wave = wcs.all_pix2world(np.arange(wcs._naxis[0]), [0.], [0.], 1)[0] * m_to_nm
idx=(np.arange(wcs._naxis[0]))
fig = plt.figure(figsize=(10, 6))
#plt.plot(wave, data.mean((0, 1)))
plt.plot(idx, data.mean((0, 1)))
plt.show()



aux = iris_file[-2]
v_obs = aux.data[:, aux.header['OBS_VRIX']]
v_obs /= 1000.  # convert to km/s
plt.plot(v_obs)
plt.ylabel("Orbital velocity (km/s)")
plt.xlabel("Scan number")
plt.close()

fig = plt.figure(figsize=(6, 10))
plt.imshow(data[..., 64].T, vmin=-32000, vmax=-26600, aspect=0.5)
plt.show()


c_kms = c / 1000.
wave_shift = - v_obs * wave[64] / (c_kms)
# linear interpolation in wavelength, for each scan
for i in range(iris_file[0].header['NRASTERP']):
    tmp = interp1d(wave - wave_shift[i], data[i], bounds_error=False)
    data[i] = tmp(wave)
fig = plt.figure(figsize=(6, 10))
plt.imshow(data[..., 64].T, vmin=-32000, vmax=-26600, aspect=0.5)
plt.savefig('shifted.png')
plt.show()



mg_k_centre = 279.63521  # in nm
pos = 100  # in km/s around line centre
print(wave)

velocity =  (wave - mg_k_centre) * c_kms / mg_k_centre
index_p = np.argmin(np.abs(velocity - pos))
index_m = np.argmin(np.abs(velocity + pos))
doppl = data[..., index_m] - data[..., index_p]
fig = plt.figure(figsize=(6, 10))
plt.imshow(doppl.T, cmap='gist_gray', aspect=0.5,vmin=-700, vmax=700,)
plt.savefig('doppler.png')
plt.show()


'''
fig = plt.figure(figsize=(6, 10))
plt.imshow(data[..., 173].T, vmin=-32000, vmax=-31600, aspect=0.5)
plt.show()

plt.imshow(data[..., 350].T, vmin=-32000, vmax=-31600, aspect=0.5)
aux = iris_file[-2]
v_obs = aux.data[:, aux.header['OBS_VRIX']]
v_obs /= 1000.  # convert to km/s
plt.plot(v_obs)
plt.ylabel("Orbital velocity (km/s)")
plt.xlabel("Scan number")'''