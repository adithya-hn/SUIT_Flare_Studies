from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pathlib

pathlib.Path('diff_img').mkdir(exist_ok=True,parents=True)
pathlib.Path('prop_diff_img').mkdir(exist_ok=True,parents=True)

images=fits.getdata('NB04_aligned_cube.fits')
nt,ny,nx=images.shape
print(nt,ny,nx)
prev_mask=[]
for i in range(78):
    img=images[i+1,:,:]-images[i,:,:]
    med = np.median(img)
    mad = np.median(np.abs(img - med))
    thresh = med + 4 * 1.4826 * mad

    fig,ax=plt.subplots()
    im=ax.imshow(img, origin='lower', cmap='gray')
    cbar = fig.colorbar(im, ax=ax)
    plt.savefig(f'diff_img/{i}.png')
    plt.close()

    mask=img>thresh
    ms=mask.copy()
    if i!=0:
        mask_=(img>0)*prev_mask
        mask=mask+mask_


    #fig,ax=plt.subplots()
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    im=axes[0].imshow(mask, origin='lower', cmap='gray')
    #cbar = fig.colorbar(im, ax=axes[0])
    axes[1].imshow(ms, origin='lower', cmap='gray')
    plt.savefig(f'prop_diff_img/{i}.png')
    plt.close()
    prev_mask=mask



    