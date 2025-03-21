
#===============================
# NB03 refference image added
#===============================

import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
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

start = timeit.default_timer()
now = datetime.datetime.now()-timedelta(days=1)
fol_nm='/home/adi/Desktop/Flare_rois/case12/nov_08_flare_data/Processed/'

jpg_fold=fol_nm+'/'+'Coloured_ROI_imgs'
algn_dir=fol_nm+'/'+'Aligned_images'

#Filter='NB03'
Filters=['NB03','NB04','NB08','BB01','BB02','BB03']
pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
search_fold='/home/adi/Desktop/Flare_rois/case12/nov_08_flare_data/' 
for fltr in Filters:
    if fltr=='NB03':
        Ref_idx=92

    if fltr=='NB04':
        Ref_idx=113

    if fltr=='NB08':
        Ref_idx=72

    ref_img_=sunpy.map.Map('/home/adi/Desktop/Flare_rois/case12/nov_08_flare_data/SUT_T24_1592_000646_Lev1.0_2024-11-08T02.59.56.250_0983NB03.fits')
    #ref_img_.peek()
    pixel_pos = np.argwhere(ref_img_.data == ref_img_.data.max()) * u.pixel
    hpc_max = ref_img_.wcs.pixel_to_world(pixel_pos[:, 1], pixel_pos[:, 0])
    #m_coord= SkyCoord(Tx=(50, 150) * u.arcsec, Ty=(-125, -200) * u.arcsec, frame=ref_img_.coordinate_frame)
    #coords = SkyCoord(Tx=(-100, 300) * u.arcsec, Ty=(0, -400) * u.arcsec, frame=ref_img_.coordinate_frame)

    fig = plt.figure()
    ax = fig.add_subplot(projection=ref_img_)
    ref_img_.plot(axes=ax)
    ax.plot_coord(hpc_max, color='white', marker='x', markersize=15)
    plt.show()

    hpc_coords = all_coordinates_from_map(ref_img_)
    r_mask = np.sqrt((hpc_coords.Tx - hpc_max.Tx) ** 2 +
                    (hpc_coords.Ty - hpc_max.Ty) ** 2) / ref_img_.rsun_obs
    mask = ma.masked_less_equal(r_mask, 0.03)
    #print(r_mask)
    scaled_map = sunpy.map.Map(ref_img_.data, ref_img_.meta, mask=mask.mask)
    coords = SkyCoord(Tx=(-100, 300) * u.arcsec, Ty=(0, -400) * u.arcsec, frame=scaled_map.coordinate_frame)
    ref_img=scaled_map.submap(coords)
    ref_img.peek()

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
    aligned_maps=mc_coalign(Sequence,template=ref_img,clip=False)
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
        aligned_img.meta['CROTA2']=0
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