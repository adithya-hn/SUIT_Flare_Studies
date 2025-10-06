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

for fltr in Filters:
    plot_data=[]
    tot_count=[]
    qs=[]
    dates=[]
    test_point=[]
    qsstd=[]
    

    search_fold=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case2_June02/data/processed/aligned_fits/{fltr}/' #Custom Folder

    
    print(f'Searching for {fltr} images in {search_fold} folder')
    fdir =search_fold 
    files = glob.glob(fdir + '*3'+fltr+'.fits')
    files=sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    files=files
    print('Total files:',len(files))

    ref_img=sunpy.map.Map(files[0])
    #base_fold='/home1/Data/Adithya/POC_Works/Jitter/'
    fl_date = datetime.datetime.fromisoformat(str(ref_img.date))
    fol_nm=os.getcwd() #str(fl_date.day).zfill(2)+'_'+str(fl_date.month).zfill(2)+'_'+str(fl_date.year).zfill(2)
    print(fol_nm)

    jpg_fold=fol_nm+'/'+'box_images'
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
        alned_data=suit_data*1000/int(suitMap.meta.get('CMD_EXPT'))
        suit_map=sunpy.map.Map(alned_data,img_head)

        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection=suit_map)
        #suit_map.plot(axes=ax, vmin=1000,vmax=16000)

        rotation_angle=suit_map.meta["CROTA2"]

        #AR sub box
        cen_cord      = SkyCoord(-255 * u.arcsec, -310 * u.arcsec, frame=suit_map.coordinate_frame)
        offset_frame1 = SkyOffsetFrame(origin=cen_cord, rotation=-rotation_angle*u.deg)
        width1  = 365 * u.arcsec
        height1 = 240 * u.arcsec
        coords = SkyCoord(lon=[-1/2, 1/2] * width1, lat=[-1/2, 1/2] * height1, frame=offset_frame1)

        #Quiet sub box
        center_coord4 = SkyCoord(-170 * u.arcsec, -500 * u.arcsec, frame=suit_map.coordinate_frame)
        #center_coord4 = SkyCoord(-275 * u.arcsec, -540 * u.arcsec, frame=suit_map.coordinate_frame) # point left, belwoe th ar box
        #center_coord4 = SkyCoord(-350 * u.arcsec, -570 * u.arcsec, frame=suit_map.coordinate_frame)
        offset_frame4 = SkyOffsetFrame(origin=center_coord4, rotation=-rotation_angle*u.deg)
        width4  = 40 * u.arcsec
        height4 = 40 * u.arcsec
        er_coords=SkyCoord(lon=[-1/2, 1/2] * width4, lat=[-1/2, 1/2] * height4, frame=offset_frame4)
        suit_map.draw_quadrangle(er_coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')

        suit_map.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
        
        er_box=suit_map.submap(er_coords)
        test_box=suit_map.submap(coords)
        print('QS val: ', np.mean(er_box.data))
        print('File',suit_map.meta['DATE-OBS'])
        Thresh_val= 3700*4 #np.mean(er_box.data)*3
        if fltr== 'NB08':
            Thresh_val= np.mean(er_box.data)*1.5

        Thresh_alned_data=np.where(alned_data>Thresh_val,alned_data,0)
        alignedMap=sunpy.map.Map(Thresh_alned_data,img_head)

        fl_nm=jpg_fold+f'/{fltr}'+'/Th'+os.path.basename(files[i])[:-4]+'jpg'
        
        alignedMap.plot(cmap='gray', vmin=Thresh_val-1, vmax=Tmax)
        plot_data.append(np.count_nonzero(Thresh_alned_data))
        tot_count.append(np.sum(Thresh_alned_data))
        #plot_str=str(alignedMap.date)+' - '+str(suitMap.date)
        #ax.text(50,50, plot_str, color='white', fontsize=10)
        #plt.draw()
        plt.savefig(fl_nm,dpi=300)
        plt.close()
        qs.append(np.mean(er_box.data))
        dates.append(suitMap.date.datetime)
        test_point.append(np.mean(test_box.data)) #test_box.meta.get("CMD_EXPT") test_box.meta.get("SHTR_STR"))#
        qsstd.append(np.std(er_box.data))

    dates=np.array(dates)
    qs=np.array(qs)
    test_point=np.array(test_point)
    np.savetxt(f'{fltr}_qs1_count.csv',np.c_[dates,qs,test_point,qsstd],comments='',header='Time,QS_mean,AR_mean,qs_std',delimiter=',',fmt='%s')
    plt.figure(figsize=(10, 5))
    ax=plt.subplot(111)
    ax1=ax.twinx()

    date_stamp=dates[0].strftime('%Y-%m-%d')
    plt.title(f'AR and QS Box Light curve - {fltr} ({date_stamp})')
    ax.plot(dates, test_point, marker='o',markersize=0.5,label='AR box Intensity')
    #ax1.errorbar(dates, qs,'r',yerr=qsstd, marker='o',markersize=0.5,label='QS1 box Intensity')
    markers1, caps1, bars1=ax1.errorbar(dates, qs,yerr=qsstd,fmt='r-', marker='o',markersize=0.5,label='QS1 box Intensity')
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



