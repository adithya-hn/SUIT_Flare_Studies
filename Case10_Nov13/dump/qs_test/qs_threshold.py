
#Date created: 15 September 2025
#Author: @adithya-hn
#Purpose: To find proper QS region and set threshold based on that

#-------------------------------------------------------------------

import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import glob
import datetime
import matplotlib.dates as mdates
import pathlib
from astropy.coordinates import SkyCoord
import numpy as np
from scipy.ndimage import shift
from scipy.ndimage import gaussian_filter
from astropy.coordinates import SkyCoord, SkyOffsetFrame

#---------Input parameters-----------------

ArTx= 100
ArTy= -235

arW=420
arH=240

qs1Tx=-50
qs1Ty=-370

qs2Tx=0
qs2Ty=-65

qs3Tx=300
qs3Ty=-70

Filters=['NB04','NB08']
search_fold=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case10_Nov13/data/1600_aligned/' #Custom Folder

#------------------------------------------


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

for fltr in Filters:
    plot_data=[]
    tot_count=[]
    qs1=[]
    qs2=[]
    qs3=[]
    dates=[]
    test_point=[]
    shut_angle=[]
    exposure=[]
    qs1_er=[]
    qs2_er=[]
    qs3_er=[]
    ar_er=[]
    
    print(f'Searching for {fltr} images in {search_fold} folder')
    fdir =search_fold 
    files = glob.glob(fdir + '*3'+fltr+'.fits')
    files=sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    files=files
    print('Total files:',len(files))

    ref_img=sunpy.map.Map(files[0])
    fl_date = datetime.datetime.fromisoformat(str(ref_img.date))
    fol_nm=os.getcwd()
    print(fol_nm)

    jpg_fold=fol_nm+'/'+'box_images'
    pathlib.Path(f'Threshold_box/{fltr}').mkdir(parents=True, exist_ok=True)
    pathlib.Path(jpg_fold+f'/{fltr}').mkdir(parents=True, exist_ok=True) # this can create required parent folders too

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
        suit_data=suitMap.data#gaussian_filter(suitMap.data,sigma=sigma)f
        alned_data=suit_data*1000/int(suitMap.meta.get('CMD_EXPT'))
        suit_map=sunpy.map.Map(alned_data,img_head)

        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection=suit_map)

        rotation_angle=suit_map.meta["CROTA2"]

        #AR sub box
        cen_cord      = SkyCoord(ArTx* u.arcsec,ArTy* u.arcsec, frame=suit_map.coordinate_frame)
        offset_frame = SkyOffsetFrame(origin=cen_cord, rotation=-rotation_angle*u.deg)
        width1  = arW * u.arcsec
        height1 = arH * u.arcsec
        coords = SkyCoord(lon=[-1/2, 1/2] * width1, lat=[-1/2, 1/2] * height1, frame=offset_frame)
        suit_map.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')

        #Quiet sub box
        qs_cen1 = SkyCoord(qs1Tx * u.arcsec, qs1Ty * u.arcsec, frame=suit_map.coordinate_frame)
        qs_cen2 = SkyCoord(qs2Tx * u.arcsec, qs2Ty * u.arcsec, frame=suit_map.coordinate_frame) # point left, belwoe th ar box
        qs_cen3 = SkyCoord(qs3Tx * u.arcsec, qs3Ty * u.arcsec, frame=suit_map.coordinate_frame)
        offset_frame1 = SkyOffsetFrame(origin=qs_cen1, rotation=-rotation_angle*u.deg)
        offset_frame2 = SkyOffsetFrame(origin=qs_cen2, rotation=-rotation_angle*u.deg)
        offset_frame3 = SkyOffsetFrame(origin=qs_cen3, rotation=-rotation_angle*u.deg)
        width4  = 40 * u.arcsec
        height4 = 40 * u.arcsec
        qs_coords1=SkyCoord(lon=[-1/2, 1/2] * width4, lat=[-1/2, 1/2] * height4, frame=offset_frame1)
        qs_coords2=SkyCoord(lon=[-1/2, 1/2] * width4, lat=[-1/2, 1/2] * height4, frame=offset_frame2)
        qs_coords3=SkyCoord(lon=[-1/2, 1/2] * width4, lat=[-1/2, 1/2] * height4, frame=offset_frame3)

        suit_map.draw_quadrangle(qs_coords1,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')
        suit_map.draw_quadrangle(qs_coords2,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')
        suit_map.draw_quadrangle(qs_coords3,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')
        
        
        er_box1=suit_map.submap(qs_coords1)
        er_box2=suit_map.submap(qs_coords2)
        er_box3=suit_map.submap(qs_coords3)

        test_box=suit_map.submap(coords)
        print('QS val: ', np.mean(er_box1.data))
        print('File',suit_map.meta['DATE-OBS'])
        Thresh_val= np.mean(er_box1.data)*2
        if fltr== 'NB08':
            Thresh_val= np.mean(er_box1.data)*1.25

        Thresh_alned_data=np.where(alned_data>Thresh_val,alned_data,0)
        alignedMap=sunpy.map.Map(Thresh_alned_data,img_head)
        '''
        suit_map.plot(axes=ax, vmin=1000,vmax=16000)
        fl_nm=jpg_fold+f'/{fltr}'+'/'+os.path.basename(files[i])[:-4]+'jpg'
        plt.savefig(fl_nm,dpi=300)
        plt.close()

        '''
        fl_nm='Threshold_box/'+fltr+'/'+os.path.basename(files[i])[:-4]+'jpg'
        alignedMap.plot(cmap='gray', vmin=Thresh_val-1, vmax=Tmax)
        plot_str=str(alignedMap.date)+' - '+str(suitMap.date)
        ax.text(50,50, plot_str, color='white', fontsize=10)
        plt.draw()
        plt.savefig(fl_nm,dpi=300)
        plt.close()
        
        

        plot_data.append(np.count_nonzero(Thresh_alned_data))
        tot_count.append(np.sum(Thresh_alned_data))
        qs1.append(np.mean(er_box1.data))
        qs2.append(np.mean(er_box2.data))
        qs3.append(np.mean(er_box3.data))
        dates.append(suitMap.date.datetime)
        test_point.append(np.mean(test_box.data)) #test_box.meta.get("CMD_EXPT") test_box.meta.get("SHTR_STR"))#
        shut_angle.append(suit_map.meta.get("SHTR_STR"))
        exposure.append(suit_map.meta.get("CMD_EXPT"))
        ar_er.append (np.std(test_box.data))#np.sqrt((np.sum(test_box.data))*suit_map.meta.get("CMD_EXPT")/1000)/(test_box.data).size)
        qs1_er.append(np.std(er_box1.data))#np.sqrt(np.sum(er_box1.data*suit_map.meta.get("CMD_EXPT")/1000))/(er_box1.data).size)
        qs2_er.append(np.std(er_box1.data))#np.sqrt((np.sum(er_box2.data))*suit_map.meta.get("CMD_EXPT")/1000)/(er_box2.data).size)
        qs3_er.append(np.std(er_box1.data))#np.sqrt((np.sum(er_box3.data))*suit_map.meta.get("CMD_EXPT")/1000)/(er_box3.data).size)

    dates=np.array(dates)
    qs1=np.array(qs1)
    test_point=np.array(test_point)
    np.savetxt(f'{fltr}_ar_qs_th_count.csv',np.c_[dates,qs1,qs1_er,qs2,qs2_er,qs3,qs3_er,test_point,ar_er,shut_angle,exposure,tot_count,plot_data],comments='',header='Time,QS1mean,qs1_se,QS2mean,qs2_se,QS3mean,qs3_se,AR_mean,AR_se,Shutter_Ang,Exposure,Threshold_count,Thresh_area',delimiter=',',fmt='%s')
    plt.figure(figsize=(10, 5))
    ax=plt.subplot(111)
    ax1=ax.twinx()

    date_stamp=dates[0].strftime('%Y-%m-%d')
    plt.title(f'AR and QS Box Light curve - {fltr} ({date_stamp})')
    #markers, caps, bars=ax.errorbar(dates, test_point,yerr=ar_er, marker='o',markersize=0.5,label='AR box Intensity')
    ax.plot(dates, test_point,marker='o',markersize=0.5,label='AR box Intensity')
    markers1, caps1, bars1=ax1.errorbar(dates, qs1,yerr=qs1_er,fmt='r-', marker='o',markersize=0.5,label='QS1 box Intensity')
    markers2, caps2, bars2=ax1.errorbar(dates, qs2,yerr=qs2_er,fmt='c-', marker='o',markersize=0.5,label='QS2 box Intensity')
    markers3, caps3, bars3=ax1.errorbar(dates, qs3,yerr=qs3_er,fmt='m-', marker='o',markersize=0.5,label='QS3 box Intensity')
    #[bar.set_alpha(0.3) for bar in bars]
    [bar.set_alpha(0.3) for bar in bars1]
    [bar.set_alpha(0.3) for bar in bars2]
    [bar.set_alpha(0.3) for bar in bars3]

    plt.gcf().autofmt_xdate()
    plt.xlabel("Time")
    ax1.set_ylabel("QS Intensity")
    ax.set_ylabel("AR Intensity")
    plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
    time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
    plt.gca().xaxis.set_major_formatter(time_formatter)
    plt.savefig(f'QS_boxes_intensity_{fltr}.png', dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(10, 5))
    ax=plt.subplot(111)
    ax1=ax.twinx()
    ax.plot(dates, test_point, marker='o',markersize=0.5,label='AR box Intensity')
    #markers, caps, bars=ax.errorbar(dates, test_point,yerr=ar_er, marker='o',markersize=0.5,label='AR box Intensity')
    markers1, caps1, bars1=ax1.errorbar(dates, qs1,yerr=qs1_er,fmt='r-', marker='o',markersize=0.5,label='QS1 box Intensity')
    #[bar.set_alpha(0.3) for bar in bars]
    [bar.set_alpha(0.3) for bar in bars1]
    plt.gcf().autofmt_xdate()
    plt.xlabel("Time")
    ax1.set_ylabel("QS Intensity")
    ax.set_ylabel("AR Intensity")
    plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
    time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
    plt.gca().xaxis.set_major_formatter(time_formatter)
    plt.savefig(f'box_QS1_intensity_{fltr}.png', dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(10, 5))
    plt.title(f'AR Threshold Light curve - {fltr} ({date_stamp})')
    ax=plt.subplot(111)
    ax1=ax.twinx()
    ax.plot(dates, tot_count,color='k',marker='o',markersize=0.5,label='Threshold count')
    ax1.plot(dates, plot_data,color='b',marker='o',markersize=0.5,label='Threshold area')
    plt.gcf().autofmt_xdate()
    plt.xlabel("Time")
    ax.set_ylabel("Threshold Intensity")
    plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
    time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
    plt.gca().xaxis.set_major_formatter(time_formatter)
    plt.savefig(f'Threshold_intensity_{fltr}.png', dpi=300, bbox_inches='tight')
    plt.close()




