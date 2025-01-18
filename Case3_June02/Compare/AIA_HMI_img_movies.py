import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
import pathlib
import timeit
from matplotlib.colors import LogNorm

Filters=['171','1600','131','HMI']

File_Source='/Analysis/Projects_Data/Flare_Data/June02_Flare_Data2'
IMG_dir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case3_June02/SDO_cutouts'

for fltr in Filters:

    fdir=f'{File_Source}/AIA/{fltr}/{fltr}_cutouts/'
    if fltr=='HMI':
        fdir=f'{File_Source}/HMI/HMI_cutouts/'

    Cutouts=IMG_dir+f'/{fltr}'

    pathlib.Path(Cutouts).mkdir(parents=True, exist_ok=True)
    print('Searching in: ',fdir, fltr )
    files = sorted(glob.glob(fdir+'*.fits'))
    print('No. of files found ',len(files))
    count=0
    
    Sequence = sunpy.map.Map(files, sequence=True)  
    fig = plt.figure()
    ax = fig.add_subplot(projection=Sequence.maps[0])
    ani = Sequence.plot(axes=ax)
    ani.save(f'{IMG_dir}/{fltr}.mp4')

    fig.clear()

    #'''
    for i in range(len(files)):
        try:
            suit_map=sunpy.map.Map(files[i],allow_errors=True)

            fnm=(os.path.basename(files[i]))
            fig = plt.figure(figsize=(5, 5))
            ax = fig.add_subplot(projection=suit_map)
            if fltr=='HMI':
                suit_map.plot(axes=ax)
            elif fltr=='131':
                suit_map.plot(axes=ax,vmin=0,vmax=1000)
            elif fltr=='171':
                suit_map.plot(axes=ax,vmin=0,vmax=5400)
            else:
                suit_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
            
            plt.colorbar(   ) 
            plt.savefig((f'{Cutouts}/{fnm[:-4]}'),dpi=300)
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
    #'''
    
 

