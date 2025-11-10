
#version in use
import os
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

# Parameters

#csv_path = 'Flare_files_Nov11_M1.4_case28f.dat'  # <- change this to your actual file path
fdir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case2_June02/data/processed/aligned_fits/'
fol_nm = os.getcwd() + '/lc_images/'
Filters = ['NB03','NB04','NB08'] #,'NB02','NB05'

# ROI and Background box coordinates
#cTx1, cTy1, cTx2, cTy2 = -413, -178, -285, -93
#Tx_er1, Ty_er1, Tx_er2, Ty_er2 = -174, -85, -100, -11

cTx1=-255
cTy1=-305
arH=365
arW=240

Tx_er1=-150
Ty_er1=-500
qsH=50
qsW=50

# Loop over filters
for fltr in Filters:

    files = sorted(glob.glob(fdir+fltr + '/*3'+fltr+'.fits'))
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

    fltr_count = []
    date_array = []
    fltr_count_err = []
    bx_area = []
    er_bx_area = []
    qs_box=[]
    qs_box_err=[]

    Sequence = sunpy.map.Map(files, sequence=True)
    aligned_maps = Sequence

    for i, suit_map in enumerate(Sequence):
        fnm = os.path.basename(files[i])[:-5]
        F_name = f'{fol_nm}/{fltr}/{fnm}.jpg'
        Box_fnm = f'{box_pth}/{fnm}.jpg'
        Headr_data = Sequence[0].fits_header
        Headr_data['DATE-OBS'] = str(aligned_maps[i].date)

        # Intensity scaling
        Imn, Imx = 10000, 45000
        if fltr == 'NB08': Imx = 12000
        elif fltr == 'NB04': Imx = 33000
        elif fltr == 'BB01': Imx = 21000
        elif fltr in ['BB02', 'BB03']: Imx = 26000

        # Plot full map
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection=suit_map)
        suit_map.plot(axes=ax, vmin=Imn, vmax=Imx)
        #coords = SkyCoord(Tx=(cTx1, cTx2) * u.arcsec, Ty=(cTy1, cTy2) * u.arcsec, frame=suit_map.coordinate_frame)

        rotation_angle=suit_map.meta["CROTA2"]
        
        cen_cord      = SkyCoord(-400 * u.arcsec, -300 * u.arcsec, frame=suit_map.coordinate_frame)
        offset_frame1 = SkyOffsetFrame(origin=cen_cord, rotation=-rotation_angle*u.deg)
        width1  = 280 * u.arcsec
        height1 = 225 * u.arcsec
        coords = SkyCoord(lon=[-1/2, 1/2] * width1, lat=[-1/2, 1/2] * height1, frame=offset_frame1)

        center_coord4 = SkyCoord(-325 * u.arcsec, -485 * u.arcsec, frame=suit_map.coordinate_frame)
        offset_frame4 = SkyOffsetFrame(origin=center_coord4, rotation=-rotation_angle*u.deg)
        width4  = 40 * u.arcsec
        height4 = 40 * u.arcsec
        er_coords=SkyCoord(lon=[-1/2, 1/2] * width4, lat=[-1/2, 1/2] * height4, frame=offset_frame4)
        suit_map.draw_quadrangle(er_coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')
        suit_map.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
        
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
    np.savetxt(f'{fltr}_c2_lc_data.csv',
               np.c_[date_array, fltr_count, fltr_count_err,qs_box,qs_box_err, bx_area,er_bx_area],
               delimiter=',',header='Time,AR_total,AR_count_Er,QS_total,QS_count_Er,AR_area,QS_area' ,fmt='%s')

stop = timeit.default_timer()
print('Run Time: ', (stop - start)/60, 'Mins')
