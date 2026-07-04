import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import glob
import datetime
from scipy import stats
import pathlib
from astropy.coordinates import SkyCoord
import numpy as np
from scipy.ndimage import shift
from scipy.ndimage import gaussian_filter
from astropy.coordinates import SkyCoord, SkyOffsetFrame

#Threshold values:


nb3T=11000
nb4T=11500
nb8T=3900
nb6T=86000
nb7T=295000

nb1Mx=10000
nb2Mx=20000
nb3Mx=14000
nb4Mx=15000
nb5Mx=45000
nb8Mx=4300
nb6Mx=95000
nb7Mx=230000

Filters=['NB04',]
search_fold=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/aligned_crop/' #Custom Folder

for fltr in Filters:
    plot_data=[]
    tot_count=[]
    qs=[]
    dates=[]
    shut_angl=[]
    expo=[]
    
    print(f'Searching for {fltr} images in {search_fold} folder')
    fdir =search_fold 
    files = glob.glob(fdir + '*3'+fltr+'.fits')
    files=sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    files=files
    print('Total files:',len(files))

    ref_img=sunpy.map.Map(files[0])
    norm_data=ref_img.data*1000/ref_img.meta.get("CMD_EXPT")
    data=(norm_data.flatten()).astype(int)
    mode_val = stats.mode(data, keepdims=True).mode[0]
    Thresh_val=4764*2.5#mode_val*1.05
    print('Frist img mode',mode_val)
    
    fol_nm=os.getcwd()
    print(fol_nm)
    jpg_fold=fol_nm+'/Threshold_imgs'

    pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
    pathlib.Path(jpg_fold+f'/{fltr}').mkdir(parents=True, exist_ok=True)

    sigma=1
    if fltr=='NB01':
        #Thresh_val=nb3T
        Tmax=nb1Mx

    if fltr=='NB02':
        #Thresh_val=nb3T
        Tmax=nb2Mx

    if fltr=='NB03':
        #Thresh_val=nb3T
        Tmax=nb3Mx

    if fltr=='NB05':
        #Thresh_val=nb3T
        Tmax=nb5Mx

    
    if fltr=='NB04':
        #Thresh_val=nb4T
        Tmax=nb4Mx

    if fltr=='NB08':
        #Thresh_val=nb8T
        Tmax=nb8Mx

    if fltr=='NB06':
        #Thresh_val=nb6T
        Tmax=nb6Mx
    
    if fltr=='NB07':
        #Thresh_val=nb7T
        Tmax=nb7Mx

    for i in range(len(files)):
        suitMap=sunpy.map.Map(files[i])
        img_head=suitMap.fits_header
        suit_data=suitMap.data#gaussian_filter(suitMap.data,sigma=sigma)
        norm_data=suit_data*1000/int(suitMap.meta.get('CMD_EXPT'))
        data=(norm_data.flatten()).astype(int)
        mode_val = stats.mode(data, keepdims=True).mode[0]
        
        Thresh_norm_data=np.where(norm_data>Thresh_val,norm_data,0)
        suit_map=sunpy.map.Map(Thresh_norm_data,img_head)

        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection=suit_map)
        suit_map.plot(axes=ax, vmin=Thresh_val-1,vmax=Tmax)
        plt.colorbar()

        print('QS val (mode): ',mode_val)
        print('Threshold: ',Thresh_val,'|',mode_val*1.02)
    
        fl_nm=jpg_fold+f'/{fltr}'+'/'+os.path.basename(files[i])[:-4]+'jpg'       
        plt.savefig(fl_nm)
        plt.close()

        plot_data.append(np.count_nonzero(Thresh_norm_data))
        tot_count.append(np.sum(Thresh_norm_data))
        qs.append(mode_val)
        dates.append(suitMap.date.datetime)
        shut_angl.append(suitMap.meta.get("SHTR_STR"))
        expo.append(suitMap.meta.get("CMD_EXPT"))
        
    dates=np.array(dates)
    qs=np.array(qs)
    np.savetxt(f'{fltr}_thresh_data.csv',np.c_[dates,qs,plot_data,tot_count,expo,shut_angl],comments='',header='time,mode,thresh_area,thresho_tot,exposure,shutter_angle',delimiter=',',fmt='%s')
    plt.figure(figsize=(10, 5))
    ax=plt.subplot(111)
    ax1=ax.twinx()
    ax.plot(dates, tot_count, marker='o',markersize=0.5)
    ax1.plot(dates, qs,'r', marker='o',markersize=0.5)
    plt.gcf().autofmt_xdate()
    plt.xlabel("Time")
    ax1.set_ylabel("QS (Mode) value")
    ax.set_ylabel("Threshold total")
    ax.legend('upper right')
    ax.legend('upper left')
    plt.title(f"Threshold count variation- {fltr}")
    plt.savefig(f'Threshold_{fltr}.png', dpi=300, bbox_inches='tight')
    plt.close()



