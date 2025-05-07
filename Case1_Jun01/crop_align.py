
#-----------------


# code works but looks like we have an issue in june 1st images itself in WCS, WCS cropping not working well,m need to crop manulally using pixel coordinates

#----------------


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
import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame


#---------------------------------------
#Input params

Ref_index=1
Ref_index_NB3=1
make_movie='yes'
Rf_tx=1
Rf_ty=1
Rf_width=1
rf_hight=1

#--------------------------------------

start = timeit.default_timer()
now = datetime.datetime.now()-timedelta(days=1)
fol_nm='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data'

jpg_crop=fol_nm+'/cropped/'+'coloured_rois'
crop_dir=fol_nm+'/cropped/'+'crop_fits'
jpg_fold=fol_nm+'/processed/'+'coloured_rois'
algn_dir=fol_nm+'/processed/'+'aligned_fits'



#Filter='NB03'
Filters=['NB03','NB04','NB08'] #,'BB01','BB02','BB03']
pathlib.Path(crop_dir).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_crop).mkdir(parents=True, exist_ok=True)

fdir =fol_nm+'/raw/'


for fltr in Filters:
    if fltr=='NB03' or fltr=='NB04':
        Ref_idx=Ref_index_NB3
        vminT=1000
        vmaxT=30000
    else:
        Ref_idx=Ref_index
        vminT=500
        vmaxT=8000
    pathlib.Path(jpg_crop+'/'+fltr).mkdir(parents=True, exist_ok=True)
    pathlib.Path(crop_dir+'/'+fltr).mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(fdir + '*3'+fltr+'.fits'))
    print('Total ',fltr ,' files:',len(files))
   
 
    crop_maps=[]
    for k in range(len(files)):
        suit_map=sunpy.map.Map(files[k])
        rotation_angle=suit_map.meta["CROTA2"]
        #print('Angle rotated',rotation_angle)
        center_coord4 = SkyCoord(-386 * u.arcsec, -354 * u.arcsec, frame=suit_map.coordinate_frame)
        offset_frame4 = SkyOffsetFrame(origin=center_coord4, rotation=-rotation_angle*u.deg)
        width4 = 350 * u.arcsec
        height4 =350 * u.arcsec
        rectangle4=SkyCoord(lon=[-1/2, 1/2] * width4, lat=[-1/2, 1/2] * height4, frame=offset_frame4)
        crop_suit=suit_map.submap(rectangle4)
        crop_maps.append(crop_suit)



        fig = plt.figure()
        ax = fig.add_subplot(projection=crop_suit)
        suit_map.plot(axes=ax, vmin=vminT, vmax=vmaxT)
        
        ax.plot_coord(center_coord4, "o", color="red")

        suit_map.draw_quadrangle(rectangle4,axes=ax,edgecolor="yellow",linestyle="--",linewidth=2,)
        fl_nm=jpg_crop+'/'+fltr+'/'+os.path.basename(files[k])[:-4]+'png'
        plt.savefig(fl_nm)
        #crop_suit.save(algn_dir+'/'+fltr+'/'+os.path.basename(files[k]),overwrite=True)
        plt.close()
        
        #suit_map.peek()

        

    Sequence = sunpy.map.MapSequence(crop_maps)#sunpy.map.Map(files, sequence=True)    
    Sequence.peek()
    aligned_maps=mc_coalign(Sequence,layer_index=Ref_idx,clip=False)
    print('Aligned, Saving now...')
    for j in range(len(aligned_maps)):  #account for incerted image, j=0 is refference image
        aligned_img=sunpy.map.Map(aligned_maps[j].data,aligned_maps[j].fits_header)
        aligned_img.meta['CRPIX1']=(aligned_maps[Ref_idx].meta['CRPIX1'])
        aligned_img.meta['CRPIX2']=(aligned_maps[Ref_idx].meta['CRPIX2'])
        aligned_img.meta['CRVAL1']=(aligned_maps[Ref_idx].meta['CRVAL1'])
        aligned_img.meta['CRVAL2']=(aligned_maps[Ref_idx].meta['CRVAL2'])
        aligned_img.save(crop_dir+'/'+fltr+'/'+os.path.basename(files[j]),overwrite=True) #need this for alignement refference
        
        fl_nm=jpg_fold+'/'+fltr+'/'+os.path.basename(files[j])[:-4]+'png'
        #Title_=aligned_img.meta.get('FTR_NAME') +' Filter: '+aligned_img.date
        fig=plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection=aligned_img)
        aligned_img.plot(title=False,cmap=filterColor[fltr], vmin=vminT, vmax=vmaxT)
        plt.figtext(0.04, 0.04,aligned_img.date, color='white',weight='bold', fontsize=18)
        plt.draw()
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(fl_nm,bbox_inches='tight', transparent="True", pad_inches=0)
        #plt.show()
        plt.close()
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 
 
    
