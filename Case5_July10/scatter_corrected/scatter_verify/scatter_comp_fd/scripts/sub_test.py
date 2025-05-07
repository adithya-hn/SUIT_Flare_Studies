
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

Tx_er1=-20
Ty_er1=-40
Tx_er2=40
Ty_er2=-120

fd_img=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/scatter_corrected/scatter_verify/scatter_comp_fd/data/raw/SUT_T24_0956_000465_Lev1.0_2024-07-10T08.32.18.130_0971NB04.fits')
roi_img=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/scatter_corrected/scatter_verify/scatter_comp_fd/data/raw/SUT_T24_0956_000465_Lev1.0_2024-07-10T08.32.51.846_0973NB04.fits')

Imx=30000
Imn=1000
#fnm=(os.path.basename(fd_img))[:-5]
coords = SkyCoord(Tx=(cTx1,cTx2 ) * u.arcsec, Ty=(cTy1, cTy2) * u.arcsec, frame=fd_img.coordinate_frame)
coords1= SkyCoord(Tx=(cTx1,cTx2 ) * u.arcsec, Ty=(cTy1, cTy2) * u.arcsec, frame=roi_img.coordinate_frame)
er_coords = SkyCoord(Tx=(Tx_er1,Tx_er2 ) * u.arcsec, Ty=(Ty_er1, Ty_er2) * u.arcsec, frame=fd_img.coordinate_frame)
er_coords1= SkyCoord(Tx=(Tx_er1,Tx_er2 ) * u.arcsec, Ty=(Ty_er1, Ty_er2) * u.arcsec, frame=roi_img.coordinate_frame)

cut_out=fd_img.submap(coords)
fd_qs=fd_img.submap(er_coords)
roi_cutout=roi_img.submap(coords)
roi_qs=roi_img.submap(er_coords)
print('--->',np.median(roi_qs.data),np.median(fd_qs.data))
print('--->',np.sum(roi_qs.data),np.sum(fd_qs.data))
#roi_cutout.peek()
fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(projection=roi_img)
roi_img.plot(axes=ax,vmin=Imn,vmax=Imx)
roi_img.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
roi_img.draw_quadrangle(er_coords,axes=ax,edgecolor="yellow",linestyle="-",linewidth=2,label='Region of interest')
plt.colorbar()
plt.savefig('ROI.png',dpi=300)
plt.close()

fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(projection=cut_out)
cut_out.plot(axes=ax,vmin=Imn,vmax=Imx)
cut_out.draw_quadrangle(coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')
cut_out.draw_quadrangle(er_coords1,axes=ax,edgecolor="green",linestyle="-",linewidth=2,label='Background')
plt.colorbar()
plt.savefig('Cutout.png',dpi=300)
plt.close()

diff=(roi_cutout.data[:633,:]*1000)/roi_cutout.meta.get('CMD_EXPT')-(cut_out.data*1000)/cut_out.meta.get('CMD_EXPT')
fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(projection=cut_out)
diff_map=sunpy.map.Map(diff, cut_out.meta)
diff_map.save('diff.fits',overwrite=True)
diff_map.plot(axes=ax,cmap='seismic')
plt.colorbar()
plt.savefig('Diff.png',dpi=300)
plt.close()
print(np.sum(abs(diff)))
