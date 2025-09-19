


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
fol_nm='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June02_Flare/Case2'
fol_nm=fol_nm+'/'+'SHARP_imgs'
pathlib.Path(fol_nm).mkdir(parents=True, exist_ok=True)


search_fold='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June02_Flare/SHARP_data_2/13697/'

fdir =os.listdir(search_fold)
fdir.sort()
#print(fdir)
key=['magnetogram','Shear','TOTPOT','USJH']
#param=key[0]

plot_data=[]
dates=[]
for param in key:
    fltr_count=[]
    date_array=[]
    for tstp in fdir:
        box_pth=fol_nm+'/'+param
        pathlib.Path(box_pth+'/SHARP_box').mkdir(parents=True, exist_ok=True)
        #pathlib.Path(box_pth).mkdir(parents=True, exist_ok=True)
        #print(fdir+fltr+'/'+'*3'+fltr+'.fits')
        files = sorted(glob.glob(f'{search_fold}/{tstp}/*{param}.fits'))
        print('Total files:',len(files))
        #print(files)
        if len(files)==0:
            print('Sorry.. No images')
            continue

        
        hmi_map=sunpy.map.Map(files[0])
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(projection=hmi_map)
        hmi_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
        #hmi_map.peek()
        coords = SkyCoord(lat=(-21,-16)  * u.deg,lon=(348, 354)* u.deg,frame=hmi_map.coordinate_frame,)
        hmi_map.draw_quadrangle(coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='2-element SkyCoord')
        fnm=(os.path.basename(files[0]))[:-5]

        #print('-------',f'{fol_nm}/{fnm}.jpg')
        plt.savefig(f'{box_pth}/{fnm}.jpg',dpi=300)
        plt.close()
        fig = plt.figure(figsize=(5, 5))
        suit_box=hmi_map.submap(coords)
        ax = fig.add_subplot(projection=suit_box)
        suit_box.plot(axes=ax ,clip_interval=(1, 99.99)*u.percent)
        plt.savefig(f'{box_pth}/SHARP_box/{fnm}.jpg',dpi=300)
        plt.close()
        fltr_count.append(np.sum(abs(suit_box.data)))
        date_array.append(hmi_map.date)
        
    fltr_count=np.array(fltr_count)
    #Plot_data_er=np.array(Plot_data_er)
    date_array=np.array(date_array)
    np.savetxt(f'{param}_M1.2_Light_curve_data.csv',np.c_[date_array,fltr_count],delimiter=',',fmt='%s')

    
    
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 