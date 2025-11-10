import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
from sunkit_image.coalignment import calculate_match_template_shift,apply_shifts
from datetime import timedelta
import timeit
import pathlib
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
#import ImagesToMovie_pkg
import matplotlib.image as mpimg
from PIL import Image
import pandas as pd
from scipy.ndimage import shift
from scipy.ndimage import gaussian_filter
from sunpy.time import parse_time
from astropy.time import Time
from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord
from tqdm import tqdm
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from sunpy.visualization.animator import MapSequenceAnimator
from sunpy.coordinates import Helioprojective, SphericalScreen, propagate_with_solar_surface

import warnings
warnings.simplefilter('ignore')

import logging
#sunpy.log.set_level("WARNING")
log_ = logging.getLogger('sunpy')
log_.setLevel('WARNING')


#-----------------Intial paths and params--------------------------

Filters=['171','1600']
suit_raw_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case10_Nov13/data/raw/'
aia_imgs_pth='/media/adithya/Adi_disk4/SUIT_flare_work/case10_nov13/data/aia/cut_outs/'
hmi_imgs_pth='/media/adithya/Adi_disk4/SUIT_flare_work/case10_nov13/data/hmi/HMI_cutouts/'
suit_filters=['NB03','NB08','NB04']
ref_1600='/media/adithya/Adi_disk4/SUIT_flare_work/case10_nov13/data/aia/cut_outs/1600_cutouts/aia.lev1_uv_24s.2024-11-13T145904Z.1600.image_lev1.fits'
tx1,ty1=-100,-90
tx2,ty2=320,-350
save_aligned='yes'
save_pngs='no'
save_aligned_pth='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case10_Nov13/data/1600_aligned/'
alin_fltr='NB04' #Filter to align other SUIT filters to
# Threshold levels

th_lvs2=[1800,3000,3250]  # Ca II H
th_lvs=[2000,10000,14000] # Mg II h/k

pathlib.Path(save_aligned_pth).mkdir(parents=True, exist_ok=True)


#-------------------------------------------

for fltr2 in Filters:
    
    c1_data=[]
    c2_data=[]
    dates=[]

    search_fold=suit_raw_files #Custom Folder
    product_fold=save_aligned_pth
    if fltr2=='HMI':
        base_fold=hmi_imgs_pth

    else:
        base_fold=aia_imgs_pth+f'{fltr2}_cutouts/'
        #base_fold=f'/media/adithya/Adi_disk4/SUIT_flare_work/case8_nov01/data/aia/cut_outs/{fltr2}_cutouts_suit/'
        #base_fold ='/media/adithya/Adi_disk4/SUIT_flare_work/case8_nov01/data/aia/aia_fd_data/'    
    
    print(f'Searching for {fltr2} images in {search_fold} folder')
    
    files2 = glob.glob(search_fold + '*3'+'NB08.fits')
    files2=sorted(files2, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    files = glob.glob(search_fold + '*3'+alin_fltr+'.fits')
    files =sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
  
    b_files=glob.glob(base_fold + '*.fits')

    print('Total SUIT NB08 files:',len(files2))
    print(f'Total SUIT {alin_fltr} files:',len(files))
    print('Total AIA files:',len(b_files))

    ca_map_time=[]
    base_time_array=[]
    mg_results=[]
    ca_results=[]
    for f in range(len(b_files)):##    @
        if fltr2=='HMI':
            #hmi.m_45s.20240602_023000_TAI.2.magnetogram
            base_time_array.append(datetime.datetime.strptime(os.path.basename(b_files[f])[10:25], "%Y%m%d_%H%M%S"))
        elif fltr2=='GONG':
            base_time_array.append(datetime.datetime.strptime(os.path.basename(b_files[f])[:11], "%Y%m%d%H%M%S"))
        
        elif fltr2=='1600':
            base_time_array.append(datetime.datetime.strptime(os.path.basename(b_files[f])[16:33], "%Y-%m-%dT%H%M%S"))

        else:
            base_time_array.append(datetime.datetime.strptime(os.path.basename(b_files[f])[17:33], "%Y-%m-%dT%H%M%S"))
    base_time_array=Time(parse_time(base_time_array))

    

    #---------------aia alignemnt--------
    ref_mg_map=sunpy.map.Map(files[0]) # or specific map
    ref_time=Time(parse_time(ref_mg_map.date))
    suit_pos = get_horizons_coord(-21, ref_mg_map.date)
    ref_mg_map.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
    rf_idx=np.argmin(np.abs(base_time_array - ref_time))

    ref_baseMap=sunpy.map.Map(ref_1600)

    rect_template=SkyCoord(Tx=(tx1,tx2 )* u.arcsec, Ty=(ty1,ty2 )* u.arcsec, frame=ref_baseMap.coordinate_frame)
    aia_template=ref_baseMap.submap(rect_template)
    #aia_template.peek()

    ref_sq=sunpy.map.MapSequence(ref_mg_map,ref_mg_map)
    aln_mgMap=mc_coalign(ref_sq,template=aia_template)
    fnm_=aln_mgMap[0].meta.get('F_NAME')
    #aln_mgMap[0].save(f'sigle_{fnm_}',overwrite=True)
    mg_map_set=[]
    for i in range(len(files)):
        if i==0:
            mg_map_set.append(aln_mgMap[0])
        mg_map_set.append(sunpy.map.Map(files[i]))
    mg_map_sq=sunpy.map.MapSequence(mg_map_set)

    #print('Ref img',ref_mg_map.date,ref_mg_map.meta.get('WAVELNTH'))
    ca_map_set=[]
    for i in range(len(files2)):
        if i==0:
            ca_map_set.append(aln_mgMap[0])
        ca_map=sunpy.map.Map(files2[i])
        ca_map_set.append(ca_map)
        
    ca_map_sq=sunpy.map.MapSequence(ca_map_set) #contains Mg ref_map

    #print('Ca 1 img',ca_map_sq[0].date)

    # Get Mg ref map layer index and update plate scale
    for i in range(len(ca_map_sq)):
        ca_map_sq[i].meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
        if ca_map_sq[i].meta['FTR_NAME']==alin_fltr:
            mg_ca_ref_idx=i
    
    
    for maps in mg_map_sq:
        maps.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))

    # Co-align maps
    print('Aligning images...')
    mg_aln_maps=mc_coalign(mg_map_sq,layer_index=0,clip=False)
    ca_aln_maps=mc_coalign(ca_map_sq,layer_index=mg_ca_ref_idx,clip=False) #includes referce mg map
    print('Done')

    Map_sq=[]
    caMap_sq=[]
    for maps in mg_aln_maps:
        alnMap=sunpy.map.Map(maps.data,maps.meta)
        alnMap.meta['CRPIX1']=mg_aln_maps[0].meta['CRPIX1']
        alnMap.meta['CRPIX2']=mg_aln_maps[0].meta['CRPIX2']
        alnMap.meta['CRVAL1']=mg_aln_maps[0].meta['CRVAL1']
        alnMap.meta['CRVAL2']=mg_aln_maps[0].meta['CRVAL2']
        
        Map_sq.append(alnMap)
        fname=str(maps.meta.get('F_NAME'))[:-5]+'.png'
        #print(fname)
        if save_aligned=='yes':
            alnMap.save(save_aligned_pth+maps.meta.get('F_NAME'),overwrite=True)
        if save_pngs =='yes':
            alnMap.plot()
            plt.savefig(fname)
            plt.close()
    
    for i in range(len(ca_aln_maps)):
        alnMap=sunpy.map.Map(ca_aln_maps[i].data,ca_aln_maps[i].meta)
        alnMap.meta['CRPIX1']=mg_aln_maps[0].meta['CRPIX1'] #Important to update all to same reference
        alnMap.meta['CRPIX2']=mg_aln_maps[0].meta['CRPIX2']
        alnMap.meta['CRVAL1']=mg_aln_maps[0].meta['CRVAL1']
        alnMap.meta['CRVAL2']=mg_aln_maps[0].meta['CRVAL2']
        
        if i != mg_ca_ref_idx:
            caMap_sq.append(alnMap)
        fname=str(alnMap.meta.get('F_NAME'))[:-5]+'.png'
        #print(fname)
        if save_aligned=='yes':
            alnMap.save(save_aligned_pth+alnMap.meta.get('F_NAME'),overwrite=True)
        if save_pngs =='yes':
            alnMap.plot(vmin=3500,vmax=7000)
            plt.savefig(fname)
            plt.close()

    for maps in caMap_sq:
        ca_map_time.append(maps.date)
    ca_map_time_array=Time(parse_time(ca_map_time)) 
    
    #--------------------

    fol_nm=os.getcwd()
    print(fol_nm)

    jpg_fold=fol_nm+'/'+'Contour_imgs'

    #pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)
    pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
    pathlib.Path(jpg_fold+f'/{fltr2}').mkdir(parents=True, exist_ok=True)

    for i in tqdm (range(len(files))):
        #suitMap=sunpy.map.Map(files[i]) #mg IIk  image
        suitMap=Map_sq[i]
        base_time=Time(parse_time(suitMap.date))
        idx=np.argmin(np.abs(base_time_array - base_time))
        idx2=np.argmin(np.abs(ca_map_time_array - base_time))

        CaII_Map=caMap_sq[idx2] #aligned Ca map   
        BaseMap=sunpy.map.Map(b_files[idx])
   
        #-- Derotate aia maps --
        with propagate_with_solar_surface():
            Base_img=BaseMap.reproject_to(ref_baseMap.wcs) 

        #-- Make expo-normalised map--

        #img_head=suitMap.fits_header
        norm_data=suitMap.data*1000/int(suitMap.meta.get('CMD_EXPT'))
        CaII_data=CaII_Map.data*1000/int(CaII_Map.meta.get('CMD_EXPT'))
        norm_ca_Map=sunpy.map.Map(gaussian_filter(CaII_data, sigma=1),CaII_Map.fits_header)
        norm_mg_Map=sunpy.map.Map(gaussian_filter(norm_data, sigma=1),suitMap.fits_header)


        #th_lvs2=[2000,3400,3800]
        #th_lvs=[3000,8000]

        fl_nm=jpg_fold+f'/{fltr2}'+'/'+os.path.basename(files[i])[:-4]+'jpg'
        fig=plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection=Base_img)
        Base_img.plot(cmap='gray',autoalign=True,title=str(BaseMap.date))
       
        
        norm_mg_Map.draw_contours(axes=ax, levels=th_lvs,lws=0.5,colors=['blue','pink','green'])
        norm_ca_Map.draw_contours(axes=ax, levels=th_lvs2,lws=0.5,colors=['red','skyblue','yellow'],alpha=0.7)

        plot_str='Mg II k: '+str(suitMap.date) +'\n'+ 'Ca II h: '+str(CaII_Map.date) 
        ax.text(50,50, plot_str, color='white', fontsize=10)
        plt.draw()
        #plt.colorbar()
        plt.savefig(fl_nm,dpi=300)
        plt.close()    
        # Append results to the list
        mg_th1=np.where(norm_mg_Map.data>th_lvs[1],norm_mg_Map.data,0)
        mg_th1_area=np.count_nonzero(mg_th1)
        mg_th1_count=np.sum(mg_th1)

        mg_th2=np.where(norm_mg_Map.data>th_lvs[2],norm_mg_Map.data,0)
        mg_th2_area=np.count_nonzero(mg_th2)
        mg_th2_count=np.sum(mg_th2)

        ca_th1=np.where(norm_ca_Map.data>th_lvs2[1],norm_ca_Map.data,0)
        ca_th1_area=np.count_nonzero(ca_th1)
        ca_th1_count=np.sum(ca_th1)

        ca_th2=np.where(norm_ca_Map.data>th_lvs2[2],norm_ca_Map.data,0)
        ca_th2_area=np.count_nonzero(ca_th2)
        ca_th2_count=np.sum(ca_th2)

        mg_results.append({'Date time': norm_mg_Map.date,
                     'Threshold 1 count': mg_th1_count,
                     'Threshold 1 area':mg_th1_area,
                     'Threshold 2 count':mg_th2_count,
                     'Threshold 2 area':mg_th2_area})
        
        ca_results.append({'Date time': norm_ca_Map.date,
                     'Threshold 1 count': ca_th1_count,
                     'Threshold 1 area':ca_th1_area,
                     'Threshold 2 count':ca_th2_count,
                     'Threshold 2 area':ca_th2_area})

mg_results_df= pd.DataFrame(mg_results)
ca_results_df = pd.DataFrame(ca_results)
mg_results_df.to_csv(('mgIIh_contours.csv'), index=False)
ca_results_df.to_csv(('caIIh_contours.csv'), index=False)
  

