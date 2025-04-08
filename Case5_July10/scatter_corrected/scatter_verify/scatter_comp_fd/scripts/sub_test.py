
import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from datetime import timedelta
import timeit
import pathlib
from astropy.coordinates import SkyCoord
#from colormap import filterColor
import numpy as np

cTx1=-350
cTy1=-50
cTx2=50
cTy2=-450

fd_img=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/scatter_corrected/scatter_verify/scatter_comp_fd/data/raw/SUT_T24_0956_000465_Lev1.0_2024-07-10T08.32.18.130_0971NB04.fits')
roi_img=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/scatter_corrected/scatter_verify/scatter_comp_fd/data/raw/SUT_T24_0956_000465_Lev1.0_2024-07-10T08.32.51.846_0973NB04.fits')

Imx=30000
Imn=1000
#fnm=(os.path.basename(fd_img))[:-5]
coords = SkyCoord(Tx=(cTx1,cTx2 ) * u.arcsec, Ty=(cTy1, cTy2) * u.arcsec, frame=fd_img.coordinate_frame)
coords1= SkyCoord(Tx=(cTx1,cTx2 ) * u.arcsec, Ty=(cTy1, cTy2) * u.arcsec, frame=roi_img.coordinate_frame)
cut_out=fd_img.submap(coords)
roi_cutout=roi_img.submap(coords1)
#roi_cutout.peek()
fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(projection=roi_img)
roi_img.plot(axes=ax,vmin=Imn,vmax=Imx)
roi_img.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
plt.colorbar()
plt.savefig('ROI.png',dpi=300)
plt.close()

fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(projection=cut_out)
cut_out.plot(axes=ax,vmin=Imn,vmax=Imx)
cut_out.draw_quadrangle(coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')
plt.colorbar()
plt.savefig('Cutout.png',dpi=300)
plt.close()

diff=roi_cutout.data[:633,:]-cut_out.data
fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(projection=cut_out)
diff_map=sunpy.map.Map(diff, cut_out.meta)
diff_map.plot(axes=ax,cmap='seismic')
plt.colorbar()
plt.savefig('Diff.png',dpi=300)
plt.close()
print(np.sum(abs(diff)))
