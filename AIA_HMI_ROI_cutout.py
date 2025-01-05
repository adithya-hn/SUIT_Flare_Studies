


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
Filters=['171','1600','131','HMI']


for fltr in Filters:
    
    fdir=f'/Analysis/Projects_Data/Flare_Data/June02_Flare_Data2/AIA/{fltr}/'
    
    if fltr=='HMI':
        fdir='/Analysis/Projects_Data/Flare_Data/June02_Flare_Data2/HMI'
    Cutouts=fdir+f'/{fltr}_cutouts'
    Cutouts_imgs=fdir+f'/{fltr}_cutout_imgs'

    pathlib.Path(Cutouts).mkdir(parents=True, exist_ok=True)
    pathlib.Path(Cutouts_imgs).mkdir(parents=True, exist_ok=True)
    print('Searching in: ',fdir+fltr )
    files = sorted(glob.glob(fdir+'*.fits'))

    for i in range(len(files)):
        suit_map=sunpy.map.Map(files[i])
        fnm=(os.path.basename(files[i]))
        print(fnm,fnm[:-5])
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(projection=suit_map)
        suit_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
        
        coords = SkyCoord(Tx=(-400, -220) * u.arcsec, Ty=(-180, -350) * u.arcsec, frame=suit_map.coordinate_frame)

        suit_map.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
        suit_box=suit_map.submap(coords)
        suit_box.save((f'{Cutouts}/{fnm}'),overwrite=True)
        plt.savefig((f'{Cutouts_imgs}/{fnm[:-4]}'),dpi=300)
        plt.show()
        
    
    
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 