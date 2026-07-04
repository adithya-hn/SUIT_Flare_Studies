import numpy as np
import pywt
import glob
from astropy.io import fits
from tqdm import tqdm
import matplotlib.pyplot as plt


fltr='NB04'
wavelet_name='db4'

cube = fits.getdata(f"{fltr}_aligned_cube.fits")

# choose pixel (y, x)
y, x = 300, 300
plt.figure(figsize=(4,4))
plt.title('SUIT Mg II h image')
plt.imshow(cube[10,: ,:],origin='lower')
plt.plot(y,x,'ro',markersize=3,label='chosen pixel')
plt.colorbar()
plt.legend()
plt.savefig('selected pixel.png',dpi=100)
plt.close()

ts = cube[:, y, x]    # intensity vs time
ts = ts#/ np.nanmedian(ts)   # normalisation (recommended)
coeffs = pywt.wavedec(ts, wavelet='haar')

recon = {}
def mad(x):
    return np.median(np.abs(x - np.median(x)))

for i in range(len(coeffs)):
    tmp = [np.zeros_like(c) for c in coeffs]
    # if coeffs[i]<0:
    #     coeffs[i]

    tmp[i] = coeffs[i]
    th=np.where(coeffs[i]>mad(coeffs[i]),coeffs[i],0)
    plt.plot(coeffs[i])
    plt.plot(th)
    plt.show()
    recon[i] = pywt.waverec(tmp, 'haar')
'''
fig,ax=plt.subplots(3,1,figsize=(14,8))
fig.subplots_adjust(hspace=0.5)
#ax2=ax.twinx()
ax[0].plot(ts, label='LC')
ax[0].set_ylabel('Intensity')
ax[0].set_title('1 Pixel Light curve')
#ax[1].plot(recon[0], label='approx')
ax[1].set_title('Haar wavelet co efficient reconstructed LC')
ax[1].plot(recon[1], label='C1')
ax[1].plot(recon[2], label='C2')
ax[1].plot(recon[3], label='C3')
ax[1].plot(recon[4], label='C4')
ax[1].plot(recon[5], label='C5')
ax[1].plot(recon[6], label='C6')
ax[2].set_title('Co efficient summed light curve')
ax[2].plot(recon[0]+recon[1]+recon[2]+recon[3]+recon[4]+recon[5]+recon[6],label='sum')
ax[1].set_ylabel('Co efficient values')
ax[1].legend()
ax[2].set_xlabel('Image frames (time)')
plt.legend()
plt.savefig('harr_wavelet.png',dpi=300)
plt.close()
'''

'''

wavelet = 'haar'
coeffs = pywt.wavedec(ts, wavelet)

# Threshold detail coefficients
sigma_factor = 1.0
coeffs_filt = [coeffs[0]]  # keep approximation

for c in coeffs[1:]:
    sigma = mad(c) / 0.6745
    c_f = pywt.threshold(c, sigma_factor * sigma, mode='hard')
    coeffs_filt.append(c_f)
    #plt.title()
    plt.plot(c_f)
    plt.plot(c,label='original')
    plt.legend()
    plt.show()

# Reconstruct
ts_rec = pywt.waverec(coeffs_filt, wavelet)
ts_rec = ts_rec[:len(ts)]  # safety trim

plt.figure(figsize=(12,6))
plt.plot(ts, color='k', alpha=0.6, label='Original LC')
plt.plot(ts_rec, color='r', lw=2, label='Wavelet reconstructed')
plt.legend()
plt.xlabel("Time index (85 s cadence)")
plt.ylabel("Intensity")
plt.tight_layout()
plt.savefig('Thresholded_haar_wavlet_recon_lc',dpi=300)
plt.show()
'''