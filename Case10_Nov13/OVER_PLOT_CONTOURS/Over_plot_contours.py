import os
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

#Threshold values:

nb3T=11000
nb4T=11500
nb8T=3900
nb6T=86000
nb7T=295000
nb3Mx=14000
nb4Mx=15000
nb8Mx=4300
nb6Mx=95000
nb7Mx=310000

Filters=['1600']

for fltr in Filters:
    fltr2=fltr
    c1_data=[]
    c2_data=[]
    dates=[]

    search_fold=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case10_Nov13/data/raw/' #Custom Folder
    #search_fold2=f'/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/P_corr_data/'
    if fltr2=='HMI':
        base_fold=f'/media/adithya/Adi_disk4/SUIT_flare_work/case8_nov01/data/HMI/{fltr2}_cutouts/'
    elif fltr2=='GONG':
        base_fold=f'/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/{fltr2}/{fltr2}_cutouts/'
    
    else:
        base_fold=f'/media/adithya/Adi_disk4/SUIT_flare_work/case10_nov13/data/aia/cut_outs/{fltr2}_cutouts/'
        #base_fold=f'/media/adithya/Adi_disk4/SUIT_flare_work/case8_nov01/data/aia/cut_outs/{fltr2}_cutouts_suit/'
        #base_fold ='/media/adithya/Adi_disk4/SUIT_flare_work/case8_nov01/data/aia/aia_fd_data/'    
    
    print(f'Searching for {fltr} images in {search_fold} folder')
    
    files2 = glob.glob(search_fold + '*3'+'NB08.fits')
    files2=sorted(files2, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    files = glob.glob(search_fold + '*3'+'NB04'+'.fits')
    files =sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

    #b_files=glob.glob(base_fold + '*.fits')
    b_files=glob.glob(base_fold + '*.fits')
    #b_files=sorted(b_files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

    print('Total SUIT NB08 files:',len(files2))
    print('Total SUIT NB04 files:',len(files))
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

    #ref_baseMap=sunpy.map.Map(b_files[rf_idx]) #or specific map
    ref_baseMap=sunpy.map.Map('/media/adithya/Adi_disk4/SUIT_flare_work/case10_nov13/data/aia/cut_outs/1600_cutouts/aia.lev1_uv_24s.2024-11-13T150016Z.1600.image_lev1.fits')

    '''
    target_scale = 0.698
    current_scale =0.6093729 #BaseMap.scale  # Tuple in arcsec/pixel
    factor_x = current_scale / target_scale
    factor_y = current_scale / target_scale
    new_shape = (int(ref_baseMap.data.shape[1] * factor_y),int(ref_baseMap.data.shape[0] * factor_x))
    aia_resampled = ref_baseMap.resample(new_shape * u.pixel)
    smoothed_data=gaussian_filter(aia_resampled.data, sigma=5)
    smoothed_map=sunpy.map.Map(smoothed_data,aia_resampled.fits_header)
    '''

    rect_template=SkyCoord(Tx=(-100,320 )* u.arcsec, Ty=(-90,-350 )* u.arcsec, frame=ref_baseMap.coordinate_frame)
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

    print('Ref img',ref_mg_map.date,ref_mg_map.meta.get('WAVELNTH'))
    ca_map_set=[]
    for i in range(len(files2)):
        if i==0:
            ca_map_set.append(aln_mgMap[0])
        ca_map=sunpy.map.Map(files2[i])
        ca_map_set.append(ca_map)
        
    ca_map_sq=sunpy.map.MapSequence(ca_map_set) #contains Mg ref_map

    print('Ca 1 img',ca_map_sq[0].date)

    # Get Mg ref map layer index and update plate scale
    for i in range(len(ca_map_sq)):
        ca_map_sq[i].meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
        if ca_map_sq[i].meta['FTR_NAME']=='NB04':
            mg_ca_ref_idx=i
    
    
    for maps in mg_map_sq:
        maps.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))

    # Co-align maps
    mg_aln_maps=mc_coalign(mg_map_sq,layer_index=0,clip=False)
    ca_aln_maps=mc_coalign(ca_map_sq,layer_index=mg_ca_ref_idx,clip=False) #includes referce mg map

    Map_sq=[]
    caMap_sq=[]
    for maps in mg_aln_maps:
        alnMap=sunpy.map.Map(maps.data,maps.meta)
        alnMap.meta['CRPIX1']=mg_aln_maps[0].meta['CRPIX1']
        alnMap.meta['CRPIX2']=mg_aln_maps[0].meta['CRPIX2']
        alnMap.meta['CRVAL1']=mg_aln_maps[0].meta['CRVAL1']
        alnMap.meta['CRVAL2']=mg_aln_maps[0].meta['CRVAL2']
        alnMap.plot(vmin=5000,vmax=25000)
        Map_sq.append(alnMap)
        fname=str(maps.meta.get('F_NAME'))[:-5]+'.png'
        #print(fname)
        plt.savefig(fname)
        plt.close()
    
    for i in range(len(ca_aln_maps)):
        alnMap=sunpy.map.Map(ca_aln_maps[i].data,ca_aln_maps[i].meta)
        alnMap.meta['CRPIX1']=mg_aln_maps[0].meta['CRPIX1']
        alnMap.meta['CRPIX2']=mg_aln_maps[0].meta['CRPIX2']
        alnMap.meta['CRVAL1']=mg_aln_maps[0].meta['CRVAL1']
        alnMap.meta['CRVAL2']=mg_aln_maps[0].meta['CRVAL2']
        alnMap.plot()
        if i != mg_ca_ref_idx:
            caMap_sq.append(alnMap)
        fname=str(maps.meta.get('F_NAME'))[:-5]+'.png'
        #print(fname)
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
   
        '''
        target_scale = 0.698
        current_scale =0.6093729 #BaseMap.scale  # Tuple in arcsec/pixel
        factor_x = current_scale / target_scale
        factor_y = current_scale / target_scale
        new_shape = (int(BaseMap.data.shape[1] * factor_y),int(BaseMap.data.shape[0] * factor_x))
        aia_resampled = BaseMap.resample(new_shape * u.pixel)
        smoothed_data=gaussian_filter(aia_resampled.data, sigma=10)
        smoothed_map=sunpy.map.Map(smoothed_data,aia_resampled.fits_header)

        #rect_template=SkyCoord(Tx=(-450,-100 )* u.arcsec, Ty=(44,336 )* u.arcsec, frame=ref_baseMap.coordinate_frame)
        #aia_template=ref_baseMap.submap(rect_template)
        #Sequence = sunpy.map.MapSequence(suitMap,suitMap)     

        #aligned_maps=mc_coalign(Sequence,template=aia_template)
        #algn_dir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/OVER_PLOT_CONTOURS/aligned'
        
        #for j in range(len(aligned_maps)):
        #fname=aligned_maps[0].meta.get('F_NAME')
        #aligned_maps[0].save(fname,overwrite=True)
        #print('0th_img in_sequence: ',aligned_maps[0].meta.get('WAVELNTH'))
        '''
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
        th_lvs2=[1500,2900,3100]
        th_lvs=[3500,10000,11000]
        #th_lvs=[3000,8000]

        fl_nm=jpg_fold+f'/{fltr2}'+'/'+os.path.basename(files[i])[:-4]+'jpg'
        fig=plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection=Base_img)
        Base_img.plot(cmap='gray',autoalign=True,title=str(BaseMap.date))
       
        
        norm_mg_Map.draw_contours(axes=ax, levels=th_lvs,lws=0.5,colors=['blue','pink','green'])
        norm_ca_Map.draw_contours(axes=ax, levels=th_lvs2,lws=0.5,colors=['red','skyblue','yellow'],alpha=0.7)

        plot_str='Mg II h: '+str(suitMap.date) +'\n'+ 'Ca II h: '+str(CaII_Map.date) 
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
  

