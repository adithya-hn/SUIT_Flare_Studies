
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
import warnings
def ignore_sunpy_warnings():
    # Filter out specific SunPy warning messages
    warnings.filterwarnings("ignore", category=UserWarning, module="sunpy")
    #warnings.filterwarnings("ignore", category=UserInfo, module="sunpy")

ignore_sunpy_warnings()
warnings.simplefilter('ignore')


#---------------------------------------
#Input params

Ref_index=0
Ref_index_NB3=0
make_movie='no'
Rf_tx=1
Rf_ty=1
Rf_width=1
rf_hight=1
ref_fltr='NB08'

#--------------------------------------

fol_nm='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case10_Nov13'

start = timeit.default_timer()
now = datetime.datetime.now()-timedelta(days=1)

jpg_fold=fol_nm+'/data/processed/'+'coloured_imgs'
algn_dir=fol_nm+'/data/processed/'+'aligned_fits'

#Filter='NB03'
Filters=['NB02','NB04','NB08','NB03','NB05'] #'NB03',
pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
search_fold=fol_nm+'/data/raw/' 
nb3_files = sorted(glob.glob(search_fold + f'*3{ref_fltr}.fits'))
nb3_files =sorted(nb3_files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
#ref_img_=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/raw/SUT_T24_0956_000465_Lev1.0_2024-07-10T13.31.01.128_0983NB03.fits')
ref_img_=sunpy.map.Map(nb3_files[0])
pixel_pos = np.argwhere(ref_img_.data == ref_img_.data.max()) * u.pixel
hpc_max = ref_img_.wcs.pixel_to_world(pixel_pos[:, 1], pixel_pos[:, 0])
fig = plt.figure()
ax = fig.add_subplot(projection=ref_img_)
ref_img_.plot(axes=ax)
ax.plot_coord(hpc_max, color='white', marker='x', markersize=15)
plt.show()
print(ref_img_.date)
Ref_idx=0
for fltr in Filters:
    
    pathlib.Path(jpg_fold+'/'+fltr).mkdir(parents=True, exist_ok=True) 
    pathlib.Path(algn_dir+'/'+fltr).mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(search_fold + '*3'+fltr+'.fits'))
    print('Total ',fltr ,' files:',len(files) )
    
    aln_imgs=[]
    if fltr == ref_fltr:
        Sequence = sunpy.map.Map(files, sequence=True)    
        aligned_maps=mc_coalign(Sequence,layer_index=0,clip=False) #template=ref_img,layer=0,clip=False)
        im_mx=13000
        im_mn= 1000
    
    else:
        s_maps=[]
        map_dates=[]
        flt_ref=sunpy.map.Map(files[0])
        aln_sq=sunpy.map.MapSequence([ref_img_,flt_ref])
        if aln_sq[0].meta.get('FTR_NAME') == ref_fltr:
            ref_map=aln_sq[1]
            aln_idx=0
            flt_idx=1
            print('its a second img',aln_sq[1].meta.get('FTR_NAME'),aln_sq[1].meta.get('DATE-OBS'))
        else:
            ref_map=aln_sq[0]
            aln_idx=1
            flt_idx=0
            print('its a first img',aln_sq[0].meta.get('FTR_NAME'),aln_sq[0].meta.get('DATE-OBS'))

        aln_fltrs=mc_coalign(aln_sq,layer_index=aln_idx,clip=False)

        s_maps.append(aln_fltrs[flt_idx])
        ref_date=ref_img_.meta.get('DATE-OBS')
        for i in range (1,len(files)):
            s_maps.append(sunpy.map.Map(files[i]))
        print('Total images in sequence:',len(s_maps),'frist img',s_maps[0].meta.get('FTR_NAME'),s_maps[0].meta.get('DATE-OBS'))    
        Sequence = sunpy.map.MapSequence(s_maps)
        '''#print('Total images in sequence:',len(Sequence))
        for maps in Sequence:
            map_dates.append(maps.meta.get('DATE-OBS'))
        print('1st img', Sequence[0].date)
        map_dates=np.array(map_dates)
        Ref_idx=np.where(map_dates==ref_date)[0][0]
        print('NB08 ref idx:',Ref_idx)'''
        aligned_maps=mc_coalign(Sequence,layer_index=0,clip=False)

        if fltr == 'NB03':
            im_mx=25000
            im_mn= 0
        elif fltr == 'NB04':
            im_mx=25000
            im_mn= 1000
        elif fltr == 'NB02':
            im_mx=30000
            im_mn= 1000
        elif fltr == 'NB05':
            im_mx=30000
            im_mn= 1000
        elif fltr == 'NB08':
            im_mx=13000
            im_mn= 0
            
    
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
        #print('Tot images',len(aligned_maps))
        '''if j ==Ref_idx :
            if fltr != 'NB08':
                print('Skipping NB08 reference image...')
                continue'''
        aligned_img=sunpy.map.Map(aligned_maps[j].data,aligned_maps[j].fits_header)
        aligned_img.meta['CRPIX1']=(ref_img_.meta['CRPIX1'])
        aligned_img.meta['CRPIX2']=(ref_img_.meta['CRPIX2'])
        aligned_img.meta['CRVAL1']=(ref_img_.meta['CRVAL1'])
        aligned_img.meta['CRVAL2']=(ref_img_.meta['CRVAL2'])
        #aligned_img.meta['CROTA2']=0
        fname= aligned_maps[j].meta.get('F_NAME')
        aligned_img.save(algn_dir+'/'+fltr+'/'+fname,overwrite=True) #need this for alignement refference
        aln_imgs.append(algn_dir+'/'+fltr+'/'+fname)
        fl_nm=jpg_fold+'/'+fltr+'/'+fname[:-4]+'jpg'
        #Title_=aligned_img.meta.get('FTR_NAME') +' Filter: '+aligned_img.date
        fig=plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection=aligned_img)
        aligned_img.plot(title=False,cmap=filterColor[fltr],vmin=im_mn,vmax=im_mx,axes=ax)
        plt.figtext(0.04, 0.04,aligned_img.date, color='white',weight='bold', fontsize=18)
        plt.draw()
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(fl_nm,bbox_inches='tight', transparent="True", pad_inches=0)
        #plt.show()
        plt.close()
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 