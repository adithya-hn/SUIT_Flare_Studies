
#==========================================
#Differential rotation implimented
#==========================================

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
from sunpy.coordinates import RotatedSunFrame
from tqdm import tqdm
import warnings
warnings.simplefilter('ignore')


start = timeit.default_timer()
now = datetime.datetime.now()-timedelta(days=1)
fol_nm='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/GONG'
fol_nm=fol_nm+'/'+'Imgs'
pathlib.Path(fol_nm).mkdir(parents=True, exist_ok=True)

search_fold='/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/GONG/GONG_cutouts/'
tx1=-475
tx2=-400
ty1=-270
ty2=-150
'''x1=-320
tx2=-190
ty1=-360
ty2=-250'''

cadence=20
Thresh_val=0
th1=3500
th2=3700
th3=4200

key=['H_alpha']
param=key[0]

plot_data=[]
dates=[]


box_pth=fol_nm+'/'+param
pathlib.Path(box_pth+'/Box').mkdir(parents=True, exist_ok=True)
files = sorted(glob.glob(f'{search_fold}*.fits'))
print('Total files:',len(files))
#print(files)
Sequence = sunpy.map.Map(files, sequence=True)
fltr_count=[]
date_array=[]
fltr_count1=[]
fltr_count2=[]
fltr_count3=[]
for i in tqdm(range (len(Sequence))):
    hmi_map=Sequence[i]
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(projection=hmi_map)
    hmi_map.plot(axes=ax)
    
    ct=hmi_map.date
    if i ==0:
        ref_tm=hmi_map.date
        
    delta_t=ct-ref_tm
 
    point1 = SkyCoord(tx1*u.arcsec, ty1*u.arcsec, frame=hmi_map.coordinate_frame) #ca bright box
    diffrot_point1 = SkyCoord(RotatedSunFrame(base=point1, duration=delta_t.sec*u.second))
    transformed_diffrot_point1 = diffrot_point1.transform_to(hmi_map.coordinate_frame)

    point2 = SkyCoord(tx2*u.arcsec, ty2*u.arcsec, frame=hmi_map.coordinate_frame) #
    diffrot_point2 = SkyCoord(RotatedSunFrame(base=point2, duration=delta_t.sec*u.second))
    transformed_diffrot_point2 = diffrot_point2.transform_to(hmi_map.coordinate_frame)

    coords = SkyCoord(Tx=(transformed_diffrot_point1.Tx.value, transformed_diffrot_point2.Tx.value) * u.arcsec,Ty=(transformed_diffrot_point1.Ty.value, transformed_diffrot_point2.Ty.value) * u.arcsec,frame=hmi_map.coordinate_frame,)
    hmi_map.draw_quadrangle(coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='2-element SkyCoord')
    fnm=(os.path.basename(files[i]))[:-5]

    plt.savefig(f'{box_pth}/{fnm}.jpg',dpi=300)
    plt.close()
    fig = plt.figure(figsize=(5, 5))
    suit_box=hmi_map.submap(coords)
    Img_data=abs(suit_box.data)
    #Thresh_data=np.where(Img_data< Thresh_val, 0, Img_data)#
    ThMap=sunpy.map.Map(Img_data,hmi_map.fits_header)

    Thresh_data1=np.where(Img_data< th1, 0, Img_data)
    Thresh_data2=np.where(Img_data< th2, 0, Img_data)
    Thresh_data3=np.where(Img_data< th3, 0, Img_data)

    ax = fig.add_subplot(projection=suit_box)
    ThMap.plot(axes=ax)
    plt.savefig(f'{box_pth}/Box/{fnm}.jpg',dpi=300)
    plt.close()
    
    fltr_count.append(np.sum(abs(Img_data)))
    date_array.append(hmi_map.date)
    fltr_count1.append(np.sum(abs(Thresh_data1)))
    fltr_count2.append(np.sum(abs(Thresh_data2)))
    fltr_count3.append(np.sum(abs(Thresh_data3)))
    

fltr_count=np.array(fltr_count)
date_array=np.array(date_array)
np.savetxt(f'cabox_{param}_X1.4_Light_curve_data.csv',np.c_[date_array,fltr_count,fltr_count1,fltr_count2,fltr_count3],delimiter=',',fmt='%s')

stop = timeit.default_timer()
print('Run Time: ', (stop - start)/60,'Mins') 