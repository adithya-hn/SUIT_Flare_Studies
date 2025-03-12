


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

Filters=['NB03','NB04','NB08']

fdir='/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed/Aligned_images/'

for fltr in Filters:  
    print('Searching in: ',fdir+fltr )
    files = sorted(glob.glob(fdir+fltr + '/*3'+fltr+'.fits'))
   
    exp_count=[]
    date_array=[]
    mes_count=[]
    Sequence = sunpy.map.Map(files, sequence=True)
    for i in range(len(Sequence)):
        suit_map=Sequence[i]
        Headr_data=Sequence[0].fits_header
        exp_count.append(Sequence[i].meta.get('CMD_EXPT'))
        mes_count.append(Sequence[i].meta.get('MEAS_EXP'))
        date_array.append(Sequence[i].date)

    exp_count=np.array(exp_count)
    mes_count=np.array(mes_count) 
    date_array=np.array(date_array)
    '''
    plt.plot(date_array,exp_count)
    plt.plot(date_array,mes_count)
    plt.xlabel('Date')
    plt.ylabel('Exposure Time')
    plt.title(f'{fltr} Filter Exposure Time')
    plt.legend(['CMD_EXPT','MEAS_EXP'])
    plt.savefig(f'{fltr}_M1.0_exp_time.png')
    plt.show()'''
    np.savetxt(f'{fltr}_M1.0_exp_time.csv',np.c_[date_array,exp_count,mes_count],delimiter=',',fmt='%s')

    
    
    
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 