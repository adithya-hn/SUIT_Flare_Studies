


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
import warnings
warnings.simplefilter('ignore')


start = timeit.default_timer()


Filters=['GONG']
File_Source='/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2'

tx1=-300
tx2=300
ty1=-450
ty2=0


for fltr in Filters:
    
    fdir=f'{File_Source}/{fltr}/'
    
    if fltr=='HMI':
        fdir=f'{File_Source}/HMI/'
    Cutouts=fdir+f'{fltr}_cutouts'
    Cutouts_imgs=fdir+f'{fltr}_cutout_imgs'

    pathlib.Path(Cutouts).mkdir(parents=True, exist_ok=True)
    pathlib.Path(Cutouts_imgs).mkdir(parents=True, exist_ok=True)
    print('Searching in: ',fdir, fltr )
    files = sorted(glob.glob(fdir+'*.fits.fz'))
    print('No. of files found ',len(files))
    count=0

    for i in range(len(files)):
        try:
            suit_map=sunpy.map.Map(files[i],allow_errors=True)

            fnm=(os.path.basename(files[i])[:-3])
            
            fig = plt.figure(figsize=(5, 5))
            ax = fig.add_subplot(projection=suit_map)
            if fltr=='HMI':
                suit_map.plot(axes=ax)
            else:
                suit_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
            
            coords = SkyCoord(Tx=(tx1, tx2) * u.arcsec, Ty=(ty1,ty2) * u.arcsec, frame=suit_map.coordinate_frame)

            suit_map.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
            suit_box=suit_map.submap(coords)
            suit_map.draw_limb(axes=ax)
            suit_map.draw_grid(axes=ax)
            
            suit_box.save(f'{Cutouts}/{fnm}',overwrite=True)
            plt.savefig((f'{Cutouts_imgs}/{fnm[:-4]}'),dpi=300)
            print('[',i,'/',len(files),']',fnm)
            if i==0:
                plt.close()
            else:
                plt.close()

        except RuntimeError:
            print('Could not read image :',os.path.basename(files[i]))
            count+=1
            pass

        
    print('Images with issues:',count)    
    
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 