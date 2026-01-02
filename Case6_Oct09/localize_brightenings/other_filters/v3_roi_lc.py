#Data created on 11th June 2025
#Author : Adithya HN

#Script ti create light curves for different filters in a specified ROI and background box
#The ROI and background box coordinates can be modified as per requirement

#Implementing ROI counter rotation to keep account for submap issues.



import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import pandas as pd
import datetime
import timeit
import pathlib
from astropy.coordinates import SkyCoord
import numpy as np
import glob
from astropy.coordinates import SkyCoord, SkyOffsetFrame

start = timeit.default_timer()


# ----------Input Parameters------

fdir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/aligned_crop/'
fol_nm = os.getcwd() + '/lc_images/'
Filters = ['NB02','NB05','NB06','NB07']

# cTx1=-200
# cTy1=-220
# arW=250
# arH=195

cTx1=140
cTy1=130
arW=70
arH=50

Tx_er1=-50
Ty_er1=-60
qsH=40
qsW=40

#------------------------




for fltr in Filters:

    files = sorted(glob.glob(fdir+ '*3'+fltr+'.fits'))
    if not files:
        print(f"No files found for filter {fltr}")
        continue

    print('Processing filter:', fltr)
    plot_data = []
    Plot_data_er = []
    dates = []
    box_pth = f'{fol_nm}/{fltr}/Box'
    pathlib.Path(f'{fol_nm}/{fltr}').mkdir(parents=True, exist_ok=True)
    pathlib.Path(box_pth).mkdir(parents=True, exist_ok=True)
    pathlib.Path('csv_files').mkdir(parents=True, exist_ok=True)

    fltr_count = []
    date_array = []
    fltr_count_err = []
    qs_box=[]
    qs_box_err=[]
    bx_area = []
    er_bx_area = []

    Sequence = sunpy.map.Map(files, sequence=True)
    aligned_maps = Sequence

    for i, suit_map in enumerate(Sequence):
        fnm = os.path.basename(files[i])[:-5]
        F_name = f'{fol_nm}/{fltr}/{fnm}.jpg'
        Box_fnm = f'{box_pth}/{fnm}.jpg'
        Headr_data = Sequence[0].fits_header
        Headr_data['DATE-OBS'] = str(aligned_maps[i].date)

        # Intensity scaling
        Imn, Imx = 5000, 45000
        if fltr == 'NB08': Imx = 12000
        elif fltr == 'NB04': Imx = 33000
        elif fltr == 'BB01': Imx = 21000
        elif fltr in ['BB02', 'BB03']: Imx = 26000

        # Plot full map
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection=suit_map)
        suit_map.plot(axes=ax, vmin=Imn, vmax=Imx)
        rotation_angle=suit_map.meta["p_angle"]
        
        cen_cord      = SkyCoord(cTx1 * u.arcsec, cTy1 * u.arcsec, frame=suit_map.coordinate_frame)
        offset_frame1 = SkyOffsetFrame(origin=cen_cord, rotation=-rotation_angle*u.deg)
        width1  = arW * u.arcsec
        height1 = arH * u.arcsec
        coords = SkyCoord(lon=[-1/2, 1/2] * width1, lat=[-1/2, 1/2] * height1, frame=offset_frame1)

        center_coord4 = SkyCoord(Tx_er1 * u.arcsec, Ty_er1 * u.arcsec, frame=suit_map.coordinate_frame)
        offset_frame4 = SkyOffsetFrame(origin=center_coord4, rotation=-rotation_angle*u.deg)
        width4  = qsW * u.arcsec
        height4 = qsH * u.arcsec
        er_coords=SkyCoord(lon=[-1/2, 1/2] * width4, lat=[-1/2, 1/2] * height4, frame=offset_frame4)

        # bl = base_drot.pixel_to_world(bottom_left[0], bottom_left[1] )
        # tr = base_drot.pixel_to_world(top_right[0] , top_right[1] )

        #suit_map.draw_quadrangle(er_coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=1,label='Background')
        suit_map.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=1,label='Region of interest')
        
        plt.colorbar()
        plt.savefig(F_name, dpi=300)
        plt.close()

        # Plot submaps
        suit_box = suit_map.submap(coords)
        er_box = suit_map.submap(er_coords)
        er_area = np.prod(er_box.data.shape)
        L, H = suit_box.data.shape
        
        

        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(projection=suit_box)
        suit_box.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
        plt.savefig(Box_fnm, dpi=300)
        plt.close()

        # Light curve values
        exposure = Sequence[i].meta.get('CMD_EXPT')
        fltr_count.append(np.sum(suit_box.data * 1000 / exposure))
        fltr_count_err.append(np.sqrt(np.sum(suit_box.data))*1000/exposure)

        qs_box.append(np.sum(er_box.data * 1000 / exposure))
        qs_box_err.append(np.sqrt(np.sum(er_box.data))*1000/exposure)

        #fltr_count_err.append(np.sum(er_box.data * 1000 / exposure) / er_area)
        date_array.append(Sequence[i].date)
        bx_area.append(L * H)
        er_bx_area.append(er_area)
   

    # Save light curve
    np.savetxt(f'csv_files/{fltr}_c4_lc_data.csv',
               np.c_[date_array, fltr_count, fltr_count_err,qs_box,qs_box_err, bx_area,er_bx_area],
               delimiter=',',header='Time,AR_total,AR_count_Er,QS_total,QS_count_Er,AR_area,QS_area',comments='' ,fmt='%s')

stop = timeit.default_timer()
print('Run Time: ', (stop - start)/60, 'Mins')
