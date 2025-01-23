
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
fol_nm='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case3_June02/Mag/Full_AR'
fol_nm=fol_nm+'/'+'Imgs'
pathlib.Path(fol_nm).mkdir(parents=True, exist_ok=True)

search_fold='/Analysis/Projects_Data/Flare_Data/June02_Flare_Data2/HMI/HMI_cutouts/'
tx1=-470 #AR
tx2=-30
ty1=-500
ty2=-150

'''
tx1=-420 #AR
tx2=-380
ty1=-225
ty2=-320

tx1=-470 #AR
tx2=-30
ty1=-500
ty2=-150

tx1=-325 #ca box
tx2=-200
ty1=-250
ty2=-170
x1=-320 #core
tx2=-190
ty1=-360
ty2=-250'''

cadence=45
Thresh_val=100
th1=100
th2=500
th3=1000

key=['magnetogram']
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
fltr_count4=[]
fltr_count5=[]
for i in tqdm(range (len(Sequence))):
    hmi_map_=Sequence[i]
    if hmi_map_.meta.get('CROTA2')>10:
        hmi_map = hmi_map_.rotate(order=3)
    else:
        hmi_map=hmi_map_

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(projection=hmi_map)
    hmi_map.plot(axes=ax)
 
    point1 = SkyCoord(tx1*u.arcsec, ty1*u.arcsec, frame=hmi_map.coordinate_frame) #ca bright box
    diffrot_point1 = SkyCoord(RotatedSunFrame(base=point1, duration=cadence*(i+1)*u.second))
    transformed_diffrot_point1 = diffrot_point1.transform_to(hmi_map.coordinate_frame)

    point2 = SkyCoord(tx2*u.arcsec, ty2*u.arcsec, frame=hmi_map.coordinate_frame) #
    diffrot_point2 = SkyCoord(RotatedSunFrame(base=point2, duration=cadence*(i+1)*u.second))
    transformed_diffrot_point2 = diffrot_point2.transform_to(hmi_map.coordinate_frame)

    coords = SkyCoord(Tx=(transformed_diffrot_point1.Tx.value, transformed_diffrot_point2.Tx.value) * u.arcsec,Ty=(transformed_diffrot_point1.Ty.value, transformed_diffrot_point2.Ty.value) * u.arcsec,frame=hmi_map.coordinate_frame,)
    hmi_map.draw_quadrangle(coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='2-element SkyCoord')
    fnm=(os.path.basename(files[i]))[:-5]

    plt.savefig(f'{box_pth}/{fnm}.jpg',dpi=300)
    plt.close()
    fig = plt.figure(figsize=(5, 5))
    suit_box=hmi_map.submap(coords)
    Img_data=suit_box.data
    l,h=Img_data.shape
    area=l*h
    Thresh_data=np.where(abs(Img_data)< Thresh_val, 0, Img_data)#
    ThMap=sunpy.map.Map(Thresh_data,hmi_map.fits_header)

    Thresh_data1=np.where(abs(Img_data)< th1, 0, abs(Img_data))
    Thresh_data2=np.where(abs(Img_data)< th2, 0, abs(Img_data))
    Thresh_data3=np.where(abs(Img_data)< th3, 0, abs(Img_data))
    pos_th=np.where(Img_data< th1, 0, Img_data)
    neg_th=np.where(Img_data> -th1, 0, Img_data)

    ax = fig.add_subplot(projection=suit_box)
    ThMap.plot(axes=ax)
    plt.savefig(f'{box_pth}/Box/{fnm}.jpg',dpi=300)
    plt.close()
    
    fltr_count.append(np.sum(abs(Thresh_data))*area*1.33e15)
    date_array.append(hmi_map.date)
    fltr_count1.append(np.sum(abs(Thresh_data1))*area*1.33e15)
    fltr_count2.append(np.sum(abs(Thresh_data2))*area*1.33e15)
    fltr_count3.append(np.sum(abs(Thresh_data3))*area*1.33e15)
    fltr_count4.append(np.sum(abs(pos_th))*area*1.33e15)
    fltr_count5.append(np.sum(abs(neg_th))*area*1.33e15)

print(area)
fltr_count=np.array(fltr_count)
date_array=np.array(date_array)
np.savetxt(f'fullar_900th{param}_M2.1_Light_curve_data.csv',np.c_[date_array,fltr_count,fltr_count1,fltr_count2,fltr_count3],delimiter=',',fmt='%s')

stop = timeit.default_timer()
print('Run Time: ', (stop - start)/60,'Mins') 