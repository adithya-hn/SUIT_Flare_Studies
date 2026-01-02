from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt



masks=fits.getdata('NB04_db_transient_reg_cube.fits')
images=fits.getdata('NB04_aligned_cube.fits')
nt,ny,nx=images.shape
print(nt,ny,nx)
for i in range(62):
    img=images[i,:,:]
    mask=masks[i,:,:]
    fig,ax=plt.subplots()
    im=ax.imshow(img, origin='lower', cmap='gray',vmin=6000,vmax=25000)
    cbar = fig.colorbar(im, ax=ax)

    # if any transient pixels exist
    if mask.any():
        ax.contour(mask, levels=[0.5], colors='red', linewidths=1)
    plt.savefig(f'db2_wvlet_cont/{i}.png')
    plt.close()
    

