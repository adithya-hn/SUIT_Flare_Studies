import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import glob
import datetime

import pathlib
from astropy.coordinates import SkyCoord
import numpy as np
from scipy.ndimage import shift
from scipy.ndimage import gaussian_filter

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
Filters=['NB03','NB08']

cTx1=-370
cTy1=-500
cTx2=0
cTy2=-200

Tx_er1=-175
Ty_er1=-515
Tx_er2=-130
Ty_er2=-475


for fltr in Filters:
    plot_data=[]
    tot_count=[]
    thres_val_array=[]
    dates=[]

    search_fold=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case3_June02/data/1600_aligned/' #Custom Folder

    
    print(f'Searching for {fltr} images in {search_fold} folder')
    fdir =search_fold 
    files = glob.glob(fdir + '*3'+fltr+'.fits')
    files=sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    files=files
    print('Total files:',len(files))

    ref_img=sunpy.map.Map(files[0])
    fl_date = datetime.datetime.fromisoformat(str(ref_img.date))
    fol_nm=os.getcwd() 

    jpg_fold=fol_nm+'/'+'Threshold_imgs'
    pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
    pathlib.Path(jpg_fold+f'/{fltr}').mkdir(parents=True, exist_ok=True)

    sigma=1
    if fltr=='NB03':
        Thresh_val=nb3T
        Tmax=nb3Mx
    
    if fltr=='NB04':
        Thresh_val=nb4T
        Tmax=nb4Mx

    if fltr=='NB08':
        Thresh_val=nb8T
        Tmax=nb8Mx

    if fltr=='NB06':
        Thresh_val=nb6T
        Tmax=nb6Mx
    
    if fltr=='NB07':
        Thresh_val=nb7T
        Tmax=nb7Mx

    for i in range(len(files)):
        suitMap=sunpy.map.Map(files[i])
        img_head=suitMap.fits_header
        suit_data=suitMap.data#gaussian_filter(suitMap.data,sigma=sigma)
        #print('before norm:',suit_data[100,100])
        norm_img=suit_data*1000/int(suitMap.meta.get('CMD_EXPT'))
        suit_map=sunpy.map.Map(norm_img,img_head)

        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection=suit_map)
        #suit_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
        coords = SkyCoord(Tx=(-580, -390) * u.arcsec, Ty=(-320, -160) * u.arcsec, frame=suit_map.coordinate_frame)
        suit_map.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
        coords = SkyCoord(Tx=(cTx1, cTx2) * u.arcsec, Ty=(cTy1, cTy2) * u.arcsec, frame=suit_map.coordinate_frame)
        er_coords = SkyCoord(Tx=(Tx_er1, Tx_er2) * u.arcsec, Ty=(Ty_er1, Ty_er2) * u.arcsec, frame=suit_map.coordinate_frame)
        suit_map.draw_quadrangle(er_coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')
        #suit_map.peek()
        er_box=suit_map.submap(er_coords)

        print('Test thresh val: ', np.mean(er_box.data)*3)
        Thresh_val= 10500 #np.mean(er_box.data)*3
        
        Thresh_alned_data=np.where(norm_img>Thresh_val,norm_img,0)
        alignedMap=sunpy.map.Map(Thresh_alned_data,img_head)

        fl_nm=jpg_fold+f'/{fltr}'+'/'+os.path.basename(files[i])[:-4]+'jpg'
        fig=plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection=alignedMap)
        alignedMap.plot(cmap='gray', vmin=Thresh_val-1, vmax=Tmax)
        plot_data.append(np.count_nonzero(Thresh_alned_data))
        tot_count.append(np.sum(Thresh_alned_data))
        dates.append(suitMap.date)
        thres_val_array.append(Thresh_val)
        #alignedMap.draw_limb(axes=ax)
        #alignedMap.draw_grid(axes=ax)
        plot_str=str(alignedMap.date)+' - '+str(suitMap.date)
        ax.text(50,50, plot_str, color='white', fontsize=10)
        plt.draw()
        plt.colorbar()
        plt.savefig(fl_nm)
        plt.close()    
        #print(i,' / ',len(files))
    plot_data=np.array(plot_data)
    tot_count=np.array(tot_count)
    thres_val_array=np.array(thres_val_array)
    np.savetxt(f'{fltr}_threshold_count.csv',np.c_[dates,plot_data,tot_count,thres_val_array],header='Date,area,Total,qs thresh',delimiter=',',fmt='%s')



