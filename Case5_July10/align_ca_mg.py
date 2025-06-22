
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
import numpy as np
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
fol_nm='/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed/'

jpg_fold=fol_nm+'/'+'Coloured_ROI_imgs'
algn_dir=fol_nm+'/'+'Aligned_images'

ca_img='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/raw/SUT_T24_0956_000465_Lev1.0_2024-07-10T13.30.12.420_0983NB08.fits'
mg_img='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/raw/SUT_T24_0956_000465_Lev1.0_2024-07-10T13.31.01.128_0983NB03.fits'

files=np.array([mg_img,ca_img])
ca_map=sunpy.map.Map(ca_img)
mg_map=sunpy.map.Map(mg_img)

Ref_idx=0

map_sq=sunpy.map.MapSequence(mg_map,ca_map)
print(map_sq[1].meta.get('FTR_NAME'))

aligned_maps=mc_coalign(map_sq,layer_index=1,clip=False)
for j in range(len(aligned_maps)):  #account for incerted image, j=0 is refference image
        aligned_img=sunpy.map.Map(aligned_maps[j].data,aligned_maps[j].fits_header)
        aligned_img.meta['CRPIX1']=(aligned_maps[Ref_idx].meta['CRPIX1'])
        aligned_img.meta['CRPIX2']=(aligned_maps[Ref_idx].meta['CRPIX2'])
        aligned_img.meta['CRVAL1']=(aligned_maps[Ref_idx].meta['CRVAL1'])
        aligned_img.meta['CRVAL2']=(aligned_maps[Ref_idx].meta['CRVAL2'])
        aligned_img.save(os.path.basename(files[j-1]),overwrite=True) #need this for alignement refference
        
        #aln_imgs.append(algn_dir+'/'+fltr+'/'+os.path.basename(files[j]))
        fl_nm=os.path.basename(files[j])[:-4]+'jpg'
        #Title_=aligned_img.meta.get('FTR_NAME') +' Filter: '+aligned_img.date
        fig=plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection=aligned_img)
        aligned_img.plot(title=False,clip_interval=(2, 99.95) * u.percent)
        plt.figtext(0.04, 0.04,aligned_img.date, color='white',weight='bold', fontsize=18)
        plt.draw()
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(fl_nm,bbox_inches='tight', transparent="True", pad_inches=0)
        #plt.show()
        plt.close()

'''
#Filter='NB03'
Filters=['NB03','NB04','NB08'] #,'BB01','BB02','BB03']
pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
search_fold='/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/' #'/scratch/suit_data/level1fits/2024/'+str(now.month).zfill(2)+'/'+str(now.day).zfill(2)+'/'+'normal_2k'+'/'
fdir =search_fold 
for fltr in Filters:
    if fltr=='NB03':
        Ref_idx=92

    if fltr=='NB04':
        Ref_idx=113

    if fltr=='NB08':
        Ref_idx=72

   
    pathlib.Path(jpg_fold+'/'+fltr).mkdir(parents=True, exist_ok=True)
    pathlib.Path(algn_dir+'/'+fltr).mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(fdir + '*3'+fltr+'.fits'))
    print('Total ',fltr ,' files:',len(files))
   
    aln_imgs=[]
    Sequence = sunpy.map.Map(files, sequence=True)    
    aligned_maps=mc_coalign(Sequence,layer_index=Ref_idx,clip=False)
   
    print('Aligned, Saving now...')
    for j in range(len(aligned_maps)):  #account for incerted image, j=0 is refference image
        aligned_img=sunpy.map.Map(aligned_maps[j].data,aligned_maps[j].fits_header)
        aligned_img.meta['CRPIX1']=(aligned_maps[Ref_idx].meta['CRPIX1'])
        aligned_img.meta['CRPIX2']=(aligned_maps[Ref_idx].meta['CRPIX2'])
        aligned_img.meta['CRVAL1']=(aligned_maps[Ref_idx].meta['CRVAL1'])
        aligned_img.meta['CRVAL2']=(aligned_maps[Ref_idx].meta['CRVAL2'])
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
'''