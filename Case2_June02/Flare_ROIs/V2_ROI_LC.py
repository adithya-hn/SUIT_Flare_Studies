


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

fol_nm=os.getcwd()#'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June01_Flare/Image_data/'
Filters=['NB03','NB04','NB08']


dat_file='Flare_files_Jun2_M2.1.dat'

for fltr in Filters:
    plot_data=[]
    dates=[]
    box_pth=fol_nm+'/'+fltr+'/'+'/img2.1/2.1Box'
    pathlib.Path(fol_nm+'/'+fltr+'/img2.1').mkdir(parents=True, exist_ok=True)
    pathlib.Path(box_pth).mkdir(parents=True, exist_ok=True)
    #print(fdir+fltr+'/'+'*3'+fltr+'.fits')

    #files = sorted(glob.glob(fdir+'/'+fltr+'/'+'*3'+fltr+'.fits'))
    files=[]
    with open(dat_file, 'r') as file:
        files_ = file.read().splitlines()
    for fls in files_:
         #print(f'{fltr}.fits')
         if fls.endswith(f'{fltr}.fits'):
            files.append(fls)

    print('Total ',fltr ,' files:',len(files))
    if len(files)==0:
        print('Sorry.. No images')
        continue

    fltr_count=[]
    date_array=[]
    Sequence = sunpy.map.Map(files, sequence=True)

    Align_LC=True
    temp_lyr=0

    if Align_LC==True:
        print('Aligning image')
        aligned_maps=mc_coalign(Sequence,layer_index=temp_lyr,clip=False)
        
    for i in range(len(Sequence)):
        suit_map=Sequence[i]
        fnm=(os.path.basename(files[i]))[:-5]
        F_name=f'{fol_nm}/{fltr}/img2.1/{fnm}.jpg'
        Box_fnm=f'{box_pth}/{fnm}.jpg'
        if Align_LC==True:
            suit_map=sunpy.map.Map(aligned_maps[i].data,Sequence[0].fits_header)
            F_name=f'{fol_nm}/{fltr}/img2.1/algn_{fnm}.jpg'
            Box_fnm=f'{box_pth}/algn_{fnm}.jpg'

        
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection=suit_map)
        suit_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
        coords = SkyCoord(Tx=(-390, -200) * u.arcsec, Ty=(-230, -330) * u.arcsec, frame=suit_map.coordinate_frame)
        suit_map.draw_quadrangle(coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='2-element SkyCoord')
        
        #print(fnm)
        plt.savefig(F_name,dpi=300)
        plt.close()
        fig = plt.figure(figsize=(5, 5))
        suit_box=suit_map.submap(coords)
        ax = fig.add_subplot(projection=suit_box)
        suit_box.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
        plt.savefig(Box_fnm,dpi=300)
        plt.close()
        fltr_count.append(np.sum(suit_box.data/Sequence[i].meta.get('MEAS_EXP')))
        date_array.append(Sequence[i].date)
    fltr_count=np.array(fltr_count)
    plot_data.append(fltr_count)
    dates.append(date_array)
    plot_data=np.array(plot_data)
    dates=np.array(dates)
    print(plot_data.shape)
    np.savetxt(f'{fltr}_M2.1_Light_curve_data.dat',plot_data)
    np.savetxt(f'{fltr}_M2.1_date_data.dat',dates,fmt='%s')
    
    
    
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 