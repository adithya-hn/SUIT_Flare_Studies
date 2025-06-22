
#===============================
# NB03 based align option is there
#===============================

import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
from sunkit_image.coalignment import calculate_match_template_shift, apply_shifts
from datetime import timedelta
import timeit
import pathlib
from colormap import filterColor
from astropy.coordinates import SkyCoord
import numpy as np
from sunpy.map.maputils import all_coordinates_from_map, coordinate_is_on_solar_disk
import numpy.ma as ma
#---------------------------------------
#Input params

Ref_index=0
Ref_index_NB3=0
make_movie='yes'
Rf_tx=1
Rf_ty=1
Rf_width=1
rf_hight=1

#--------------------------------------

fol_nm='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10'

start = timeit.default_timer()
now = datetime.datetime.now()-timedelta(days=1)

jpg_fold=fol_nm+'/data/processed/'+'Coloured_ROI_imgs'
algn_dir=fol_nm+'/data/processed/'+'Aligned_images'

#Filter='NB03'
Filters=['NB03','NB04','NB08','NB02','NB05']
pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
search_fold=fol_nm+'/data/raw/' 
ref_img_=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/raw/SUT_T24_0956_000465_Lev1.0_2024-07-10T13.31.01.128_0983NB03.fits')
pixel_pos = np.argwhere(ref_img_.data == ref_img_.data.max()) * u.pixel
hpc_max = ref_img_.wcs.pixel_to_world(pixel_pos[:, 1], pixel_pos[:, 0])
fig = plt.figure()
ax = fig.add_subplot(projection=ref_img_)
ref_img_.plot(axes=ax)
ax.plot_coord(hpc_max, color='white', marker='x', markersize=15)
plt.show()

hpc_coords = all_coordinates_from_map(ref_img_)
r_mask = np.sqrt((hpc_coords.Tx - hpc_max.Tx) ** 2 +
                (hpc_coords.Ty - hpc_max.Ty) ** 2) / ref_img_.rsun_obs
mask = ma.masked_less_equal(r_mask, 0.06)
#print(r_mask)
scaled_map = sunpy.map.Map(ref_img_.data, ref_img_.meta)#, mask=mask.mask)
row,col=(scaled_map.data).shape
ll=scaled_map.wcs.pixel_to_world(0,0)
tr=scaled_map.wcs.pixel_to_world(row,col)
print(ll.Tx,ll.Ty,tr.Tx,tr.Ty.value)
coords = SkyCoord(Tx=(ll.Tx.value+20,tr.Tx.value-20) * u.arcsec, Ty=(ll.Ty.value+60,tr.Ty.value-60) * u.arcsec, frame=scaled_map.coordinate_frame)
ref_img=ref_img_.submap(coords)
scaled_map.plot()
scaled_map.draw_quadrangle(coords)
plt.show()
ref_img.peek()
for fltr in Filters:

    ###############################################################################
    # Let's plot the results.
    '''
    fig = plt.figure()
    ax = fig.add_subplot(projection=scaled_map)
    scaled_map.plot(axes=ax)
    plt.close()'''
    
    pathlib.Path(jpg_fold+'/'+fltr).mkdir(parents=True, exist_ok=True) 
    pathlib.Path(algn_dir+'/'+fltr).mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(search_fold + '*3'+fltr+'.fits'))
    print('Total ',fltr ,' files:',len(files) )
   
    aln_imgs=[]
    Sequence = sunpy.map.Map(files, sequence=True)    
    aligned_maps=mc_coalign(Sequence,layer_index=0,clip=False) #template=ref_img,layer=0,clip=False)
    ref_cdel=ref_img.meta['CDELT1']
    ref_img1=ref_img_.data[:-1,:-1]
    align_shift = calculate_match_template_shift(Sequence, template=ref_img1)
    shift_xPix = align_shift['x'].value / ref_cdel * -1
    shift_yPix = align_shift['y'].value / ref_cdel * -1
    #aligned_maps = apply_shifts(Sequence, yshift=shift_yPix * u.pixel, xshift=shift_xPix * u.pixel, clip=False)
    if make_movie== 'yes':
        mv_nm=fltr+'.mp4'
        fig = plt.figure()
        ax = fig.add_subplot(projection=ref_img_)
        ani = aligned_maps.plot(axes=ax,cmap=filterColor[fltr])
        #plt.axis('off')
        ani.save(mv_nm)
        plt.close()
    print('Aligned, Saving now...')
    for j in range(len(aligned_maps)):  #account for incerted image, j=0 is refference image
        aligned_img=sunpy.map.Map(aligned_maps[j].data,aligned_maps[j].fits_header)
        aligned_img.meta['CRPIX1']=(ref_img_.meta['CRPIX1'])
        aligned_img.meta['CRPIX2']=(ref_img_.meta['CRPIX2'])
        aligned_img.meta['CRVAL1']=(ref_img_.meta['CRVAL1'])
        aligned_img.meta['CRVAL2']=(ref_img_.meta['CRVAL2'])
        #aligned_img.meta['CROTA2']=0
        aligned_img.save(algn_dir+'/'+fltr+'/'+os.path.basename(files[j]),overwrite=True) #need this for alignement refference
        aln_imgs.append(algn_dir+'/'+fltr+'/'+os.path.basename(files[j]))
        fl_nm=jpg_fold+'/'+fltr+'/'+os.path.basename(files[j])[:-4]+'jpg'
        #Title_=aligned_img.meta.get('FTR_NAME') +' Filter: '+aligned_img.date
        fig=plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection=aligned_img)
        aligned_img.plot(title=False,cmap=filterColor[fltr],clip_interval=(2, 99.95) * u.percent)
        plt.figtext(0.04, 0.04,aligned_img.date, color='white',weight='bold', fontsize=18)
        plt.draw()
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(fl_nm,bbox_inches='tight', transparent="True", pad_inches=0)
        #plt.show()
        plt.close()
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 