import matplotlib.pyplot as plt

import astropy.units as u
import numpy as np
import sunpy.data.sample
import sunpy.map
from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord
import os

File='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/May30_Flare/SUT_T24_0784_000393_Lev1.0_2024-05-30T05.00.23.248_0973NB03.fits'
fnm=os.path.basename(File)


suit_map = sunpy.map.Map(File)
suit_pos = get_horizons_coord(-21, suit_map.date)
suit_map.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
suit_map.save('original.fits',overwrite=True)

crx=suit_map.meta.get('CRPIX1')
cry=suit_map.meta.get('CRPIX2')
suit_map.meta.update({'CRPIX1':2558 })
suit_map.meta.update({'CRPIX2':1222})

dt=suit_map.data
image=dt
if np.isnan(image).any():
    print('ahhhh')

if not np.isfinite(image).all():
    print('hmm')
dt[np.isnan(dt)]=0
Sum_1=np.sum(np.array(suit_map.data))
print('Sum: ',Sum_1)
aia_rotated = suit_map.rotate(angle=45 * u.deg)
print('-->',suit_map.meta.get('CRPIX1'),suit_map.meta.get('CRPIX2'))
#aia_rotated.save(fnm,overwrite=True)
print(crx,cry)
#suit_map.meta.update({'CRPIX1':crx })
#suit_map.meta.update({'CRPIX2':cry})
old_header=suit_map.fits_header
old_header['CRPIX1']=crx
old_header['CRPIX2']=cry
rot_header=aia_rotated.fits_header
#old_header['CRPIX1']=2556 #rot_header['CRPIX1']
#old_header['CRPIX2']=1222 #rot_header['CRPIX2']

rt_new_map=sunpy.map.Map(aia_rotated.data,old_header)
rt_new_map.save(fnm,overwrite=True)
rt_new_map.peek()


'''r_aia_map = sunpy.map.Map(fnm)
r_aia_map.data[np.isnan (r_aia_map.data)]=0
Sum_=np.sum(np.array(r_aia_map.data))
fig = plt.figure()
ax = fig.add_subplot(projection=r_aia_map)
r_aia_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
#r_aia_map.draw_limb(axes=ax)
#r_aia_map.draw_grid(axes=ax)
plt.show()
print(Sum_-Sum_1)

err=(Sum_-Sum_1)*100/Sum_1

print('Error',err)'''

fig = plt.figure()
ax = fig.add_subplot(projection=suit_map)
suit_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
suit_map.draw_limb(axes=ax)
suit_map.draw_grid(axes=ax)
plt.show()


fig = plt.figure()
ax = fig.add_subplot(projection=aia_rotated)
aia_rotated.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
aia_rotated.draw_limb(axes=ax)
aia_rotated.draw_grid(axes=ax)
plt.show()