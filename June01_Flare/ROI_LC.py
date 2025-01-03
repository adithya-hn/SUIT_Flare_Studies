


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
now = datetime.datetime.now()-timedelta(days=1)
fol_nm='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June01_Flare/Image_data/'#str(now.day)+'_'+str(now.month).zfill(2)+'_'+str(now.year).zfill(2)


Filters=['NB08']

search_fold='/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Processed_Data/Aligned_images/'#'/scratch/suit_data/level1fits/2024/'+str(now.month).zfill(2)+'/'+str(now.day).zfill(2)+'/'+'normal_2k'+'/'
fdir =search_fold 
plot_data=[]
dates=[]

for fltr in Filters:
    box_pth=fol_nm+'/'+fltr+'/'+'Box'
    pathlib.Path(fol_nm+'/'+fltr).mkdir(parents=True, exist_ok=True)
    pathlib.Path(box_pth).mkdir(parents=True, exist_ok=True)
    #print(fdir+fltr+'/'+'*3'+fltr+'.fits')
    files = sorted(glob.glob(fdir+'/'+fltr+'/'+'*3'+fltr+'.fits'))
    print('Total ',fltr ,' files:',len(files))
    if len(files)==0:
        print('Sorry.. No images')
        continue

    fltr_count=[]
    date_array=[]
    Sequence = sunpy.map.Map(files, sequence=True)
    for i in range(len(Sequence)):
        suit_map=Sequence[i]
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection=suit_map)
        suit_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
        
        coords = SkyCoord(Tx=(-580, -390) * u.arcsec,Ty=(-160, -320) * u.arcsec,frame=suit_map.coordinate_frame,)
        suit_map.draw_quadrangle(coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='2-element SkyCoord')
        fnm=(os.path.basename(files[i]))[:-5]
        #print(fnm)
        plt.savefig(f'{fol_nm}{fltr}/{fnm}.jpg',dpi=300)
        plt.close()
        fig = plt.figure(figsize=(5, 5))
        suit_box=suit_map.submap(coords)
        ax = fig.add_subplot(projection=suit_box)
        suit_box.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
        plt.savefig(f'{box_pth}/{fnm}.jpg',dpi=300)
        plt.close()
        fltr_count.append(np.sum(suit_box.data/Sequence[i].meta.get('MEAS_EXP')))
        date_array.append(Sequence[i].date)
    fltr_count=np.array(fltr_count)
    #print(len(fltr_count))
    plot_data.append(fltr_count)
    dates.append(date_array)
    #print(plot_data.shape)
plot_data=np.array(plot_data)
dates=np.array(dates)
#plot_data=plot_data.reshape(8,189)
#plot_data=np.vstack([plot_data,date_array])
print(plot_data.shape)
np.savetxt('NB03_Light_curve_data.dat',plot_data)
np.savetxt('NB03_date_data.dat',dates,fmt='%s')
    
    
    
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 