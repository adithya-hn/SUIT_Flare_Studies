
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

search_fold='/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Flare_data_pCorr/' 
suit_map=sunpy.map.Map('/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Flare_data_pCorr/SUT_T24_0785_000396_Lev1.0_2024-06-01T09.14.59.374_0983NB03.fits')
coords = SkyCoord(Tx=(-650, -450) * u.arcsec, Ty=(-400, -100) * u.arcsec, frame=suit_map.coordinate_frame)
suit_box=suit_map.submap(coords)
jpg_fold=search_fold+'processed/roi_imgs'
algn_dir=search_fold+'processed/aligned_images'

#Filter='NB03'
Filters=['NB03','NB04','NB08'] 
pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)

fdir =search_fold 
for fltr in Filters:
    
    pathlib.Path(jpg_fold+'/'+fltr).mkdir(parents=True, exist_ok=True)
    pathlib.Path(algn_dir+'/'+fltr).mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(fdir + '*3'+fltr+'.fits'))
    print('Total ',fltr ,' files:',len(files))
   
    aln_imgs=[]
    Sequence = sunpy.map.Map(files, sequence=True)    
    aligned_maps=mc_coalign(Sequence,template=suit_box,clip=False)
    if make_movie== 'yes':
        mv_nm=fltr+'.mp4'
        fig = plt.figure()
        ax = fig.add_subplot(projection=suit_box)
        ani = aligned_maps.plot(axes=ax,cmap=filterColor[fltr])
        #plt.axis('off')
        ani.save(mv_nm)
        plt.close()
    print('Aligned, Saving now...')
    for j in range(len(aligned_maps)):  #account for incerted image, j=0 is refference image
        aligned_img=sunpy.map.Map(aligned_maps[j].data,aligned_maps[j].fits_header)
        aligned_img.meta['CRPIX1']=(suit_box.meta['CRPIX1'])
        aligned_img.meta['CRPIX2']=(suit_box.meta['CRPIX2'])
        aligned_img.meta['CRVAL1']=(suit_box.meta['CRVAL1'])
        aligned_img.meta['CRVAL2']=(suit_box.meta['CRVAL2'])
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
        plt.close()
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 