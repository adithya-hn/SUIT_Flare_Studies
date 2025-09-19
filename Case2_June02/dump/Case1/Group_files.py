import numpy as np
import matplotlib.pyplot as plt
import glob
import pathlib
import shutil
import os

MainFolder='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June02_Flare/SHARP_data_2'
mag_files=sorted(glob.glob(MainFolder+'/*magnetogram.fits'))

#all_files=glob(MainFolder+'.fits')

for i in range(len(mag_files)):
    mag_f=os.path.basename(mag_files[i])
    fold_nm=mag_f[25:39]
    
    #print('--',fold_nm,mag_f)
    dest=MainFolder+'/'+fold_nm
    pathlib.Path(MainFolder+'/'+fold_nm).mkdir(parents=True, exist_ok=True) 
    batch_files=glob.glob(MainFolder+'/hmi.sharp_cea_720s.11297.'+fold_nm+'*')
    print(batch_files)
    for fls in batch_files:
        print(fls)
        shutil.move(fls,dest)
    
#hmi.sharp_cea_720s.11297.20240602_063600_TAI.magnetogram.fits
    
