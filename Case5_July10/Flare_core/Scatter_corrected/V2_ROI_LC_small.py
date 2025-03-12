


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
from astropy.coordinates import SkyCoord
#from colormap import filterColor
import numpy as np


start = timeit.default_timer()

fol_nm=os.getcwd()+'/Light_curve_images/'
Filters=['NB03','NB04','NB08']

fdir='/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/july10_11ut_to_18ut/Processed/Aligned_images/'
#dat_file='Flare_files_Jun2_M1.2.dat'
#Filters=['NB08']

for fltr in Filters:
    plot_data=[]
    Plot_data_er=[]
    dates=[]
    box_pth=fol_nm+'/'+fltr+'/'+'/Box'
    pathlib.Path(fol_nm+'/'+fltr).mkdir(parents=True, exist_ok=True)
    pathlib.Path(box_pth).mkdir(parents=True, exist_ok=True)
    #print(fdir+fltr+'/'+'*3'+fltr+'.fits')

    #files = sorted(glob.glob(fdir+'/'+fltr+'/'+'*3'+fltr+'.fits'))
    files=[]
  
    print('Searching in: ',fdir+fltr )
    files = sorted(glob.glob(fdir+fltr + '/*3'+fltr+'.fits'))
   
    fltr_count=[]
    date_array=[]
    fltr_count_err=[]
    bx_area=[]
    Sequence = sunpy.map.Map(files, sequence=True)

    Align_LC=''#True
    temp_lyr=0


    if Align_LC==True:
        print('Aligning image')
        aligned_maps=mc_coalign(Sequence,layer_index=temp_lyr,clip=False)
    else:
        print('Considering its already aligned...')
        aligned_maps=Sequence

    for i in range(len(Sequence)):
        suit_map=Sequence[i]
        fnm=(os.path.basename(files[i]))[:-5]
        F_name=f'{fol_nm}/{fltr}/{fnm}.jpg'
        Box_fnm=f'{box_pth}/{fnm}.jpg'
        Headr_data=Sequence[0].fits_header
        #print(Headr_data,aligned_maps[i].date)
        Headr_data['DATE-OBS']=str(aligned_maps[i].date)
        if Align_LC==True:
            suit_map=sunpy.map.Map(aligned_maps[i].data,Headr_data)
            #suit_map.meta.update({'DATE-OBS':aligned_maps[i].date})
            F_name=f'{fol_nm}/{fltr}/algn_{fnm}.jpg'
            Box_fnm=f'{box_pth}/algn_{fnm}.jpg'

        
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection=suit_map)
        suit_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
        coords = SkyCoord(Tx=(-175, -75) * u.arcsec, Ty=(-250, -150) * u.arcsec, frame=suit_map.coordinate_frame)
        #coords = SkyCoord(Tx=(-580, -390) * u.arcsec, Ty=(-320, -160) * u.arcsec, frame=suit_map.coordinate_frame)
        suit_map.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
        er_coords = SkyCoord(Tx=(-10, 90) * u.arcsec, Ty=(0, -100) * u.arcsec, frame=suit_map.coordinate_frame)
        suit_map.draw_quadrangle(er_coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')
        
        #er_coords = SkyCoord(Tx=(-290, -390) * u.arcsec, Ty=(-380, -400) * u.arcsec, frame=suit_map.coordinate_frame)
        #suit_map.draw_quadrangle(er_coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')
        
        #print(fnm)
        plt.savefig(F_name,dpi=300)
        plt.close()
        fig = plt.figure(figsize=(5, 5))
        suit_box=suit_map.submap(coords)
        er_box=suit_map.submap(er_coords)
        l,h=np.shape(er_box.data)
        er_area=l*h
        L,H=np.shape(suit_box.data)
        bx_area.append(L*H)
        ax = fig.add_subplot(projection=suit_box)
        suit_box.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
        plt.savefig(Box_fnm,dpi=300)
        plt.close()
        fltr_count.append(np.sum(suit_box.data*1000/Sequence[i].meta.get('MEAS_EXP')))
        #fltr_count_err.append(np.sqrt(np.sum(suit_box.data))*1000/Sequence[i].meta.get('MEAS_EXP')) #poisonian error
        #fltr_count_err.append(np.sqrt(er_area)*(np.std(er_box.data))*1000/Sequence[i].meta.get('MEAS_EXP'))
        fltr_count_err.append(np.sum(er_box.data*1000/Sequence[i].meta.get('MEAS_EXP'))/er_area)
        date_array.append(Sequence[i].date)
    fltr_count=np.array(fltr_count)
    fltr_count_err=np.array(fltr_count_err)
    plot_data=np.array(plot_data)
    Plot_data_er=np.array(Plot_data_er)
    dates=np.array(dates)
    bx_area=np.array(bx_area)
    
   
    np.savetxt(f'{fltr}_M1.0_Light_curve_data.csv',np.c_[date_array,fltr_count,fltr_count_err,bx_area],delimiter=',',fmt='%s')
    #np.savetxt(f'{fltr}_X1.4_date_data.dat',dates,fmt='%s')
    
    
    
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 