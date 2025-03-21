


import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from datetime import timedelta
import timeit
import pathlib
from astropy.coordinates import SkyCoord
#from colormap import filterColor
import numpy as np


start = timeit.default_timer()

fol_nm=os.getcwd()+'/lc_images/'
Filters=['NB03','NB04','NB08','BB01','BB02','BB03']

fdir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Flare_stat/case14/data/raw/'
cTx1=-460
cTy1=-260
cTx2=-380
cTy2=-120

Tx_er1=-360
Ty_er1=-80
Tx_er2=-280
Ty_er2=0
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
    files = sorted(glob.glob(fdir + '/*3'+fltr+'.fits'))
   
    fltr_count=[]
    date_array=[]
    fltr_count_err=[]
    bx_area=[]
    Sequence = sunpy.map.Map(files, sequence=True)

    aligned_maps=Sequence

    for i in range(len(Sequence)):
        suit_map=Sequence[i]
        fnm=(os.path.basename(files[i]))[:-5]
        F_name=f'{fol_nm}/{fltr}/{fnm}.jpg'
        Box_fnm=f'{box_pth}/{fnm}.jpg'
        Headr_data=Sequence[0].fits_header
        #print(Headr_data,aligned_maps[i].date)
        Headr_data['DATE-OBS']=str(aligned_maps[i].date)
  

        if fltr=='NB03':
            Imx=30000
            Imn=1000
        if fltr=='NB03':
            Imx=30000
            Imn=1000
        if fltr=='NB08':
            Imx=15000
            Imn=1000
        if fltr=='BB01':
            Imx=21000
            Imn=1000
        if fltr=='BB02':
            Imx=26000
            Imn=1000
        if fltr=='BB03':
            Imx=26000
            Imn=1000
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection=suit_map)
        suit_map.plot(axes=ax,vmin=Imn,vmax=Imx)
        coords = SkyCoord(Tx=(cTx1,cTx2 ) * u.arcsec, Ty=(cTy1, cTy2) * u.arcsec, frame=suit_map.coordinate_frame)
        er_coords = SkyCoord(Tx=(Tx_er1,Tx_er2) * u.arcsec, Ty=(Ty_er1,Ty_er2) * u.arcsec, frame=suit_map.coordinate_frame)

        suit_map.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
        suit_map.draw_quadrangle(er_coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')
        plt.colorbar()
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

        fltr_count.append(np.sum(suit_box.data*1000/Sequence[i].meta.get('CMD_EXPT')))
        fltr_count_err.append(np.sum(er_box.data*1000/Sequence[i].meta.get('CMD_EXPT'))/er_area)
        date_array.append(Sequence[i].date)

    fltr_count=np.array(fltr_count)
    fltr_count_err=np.array(fltr_count_err)
    plot_data=np.array(plot_data)
    Plot_data_er=np.array(Plot_data_er)
    dates=np.array(dates)
    bx_area=np.array(bx_area)
    
   
    np.savetxt(f'{fltr}_c14_lc_data.csv',np.c_[date_array,fltr_count,fltr_count_err,bx_area],delimiter=',',fmt='%s')
    #np.savetxt(f'{fltr}_X1.4_date_data.dat',dates,fmt='%s')
    
    
    
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 