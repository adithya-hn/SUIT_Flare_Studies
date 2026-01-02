
import numpy as np
from astropy.io import fits
import glob
import matplotlib.pyplot as plt
import pathlib

files = sorted(glob.glob("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case9_Nov13/data/raw/*NB02.fits"))
print('No. of files: ',len(files))

pathlib.Path('Spot_tracking').mkdir(parents=True, exist_ok=True) 
tracking_mode= 'patch' #patch /'ROI'
ref_roi_x, ref_roi_y = 440, 600   # example pixel location in reference ROI
rw= 50

with fits.open(files[0]) as hdul:
        data = hdul[0].data
        hdr  = hdul[0].header
        date = hdr['DATE-OBS']
        X = hdr['X1']
        Y = hdr['Y1']

i=0
for f in files:
    i+=1
    with fits.open(f) as hdul:
        data = hdul[0].data
        hdr  = hdul[0].header
        date = hdr['DATE-OBS']

    X1 = hdr['X1']
    Y1 = hdr['Y1']
    shift_x=X1-X
    shift_y=Y1-Y
    print('shift: ',shift_y,shift_x, 'pos: ',ref_roi_y-shift_y,ref_roi_x-shift_x)
    if tracking_mode=='patch':
        patch=data[(ref_roi_y-shift_y)-rw:(ref_roi_y-shift_y)+rw,(ref_roi_x-shift_x)-rw:(ref_roi_x-shift_x)+rw]
        plt.imshow(patch, origin='lower', cmap='gray')
        plt.plot(rw,rw,'or',markersize=10,alpha=0.2)
        plt.title(str(date)[:19]+' | shift pix: ('+str(shift_y)+' '+str(shift_x)+')')
        plt.savefig(f'Spot_tracking/tracking_patch_{i}.png',dpi=200)
        plt.close()
    if tracking_mode=='ROI':
        plt.imshow(data, origin='lower', cmap='gray')
        plt.plot(ref_roi_x-shift_x,ref_roi_y-shift_y,'or',markersize=10,alpha=0.2) #coordinates are fliped to accomodate origine lower
        plt.title(str(date)[:19]+' | shift pix: ('+str(shift_y)+' '+str(shift_x)+')')
        plt.savefig(f'Spot_tracking/spot_on_roi{i}.png',dpi=200)
        plt.close()
