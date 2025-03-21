


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


fol_nm=os.getcwd()+'/Light_curve_images/'
Filters=['NB03','NB04','NB08','BB01','BB02','BB03']

fdir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Flare_stat/case12/data/processed/'
f_disc=''
roi_patch=''

for fltr in Filters:

    files=[]
  
    print('Searching in: ',fdir+fltr )
    files = sorted(glob.glob(fdir+fltr + '/*3'+fltr+'.fits'))
   
    fltr_count=[]
    date_array=[]
    fltr_count_err=[]
    bx_area=[]
    Sequence = sunpy.map.Map(files, sequence=True)

    Align_LC=''#True
    temp_lyr=0
    
    suit_map=sunpy.map.Map()
    suit_disc=sunpy.map.Map()

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(projection=suit_map)
    suit_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
    coords = SkyCoord(Tx=(-175, -75) * u.arcsec, Ty=(-250, -150) * u.arcsec, frame=suit_map.coordinate_frame)
    suit_map.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
    er_coords = SkyCoord(Tx=(-10, 90) * u.arcsec, Ty=(0, -100) * u.arcsec, frame=suit_map.coordinate_frame)
    suit_map.draw_quadrangle(er_coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')

    suit_box=suit_map.submap(coords)
    er_box=suit_map.submap(er_coords)

    fltr_count.append(np.sum(suit_box.data*1000/Sequence[i].meta.get('MEAS_EXP')))
    fltr_count_err.append(np.mean(er_box.data*1000/Sequence[i].meta.get('MEAS_EXP')))
    date_array.append(Sequence[i].date)
    fltr_count=np.array(fltr_count)
    fltr_count_err=np.array(fltr_count_err)
    plot_data=np.array(plot_data)
    Plot_data_er=np.array(Plot_data_er)
    dates=np.array(dates)
    bx_area=np.array(bx_area)
    
   
    np.savetxt(f'{fltr}_M1.0_Light_curve_data.csv',np.c_[date_array,fltr_count,fltr_count_err,bx_area],delimiter=',',fmt='%s')
    #np.savetxt(f'{fltr}_X1.4_date_data.dat',dates,fmt='%s')
    
stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 