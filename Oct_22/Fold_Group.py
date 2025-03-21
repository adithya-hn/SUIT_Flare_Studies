#Author: Adithya HN
#Date: 23/02/24
#Purpose: Regroup the files in to folder for DEM pairs, wavlength of 94A is used as folder name

#----------------------------------------------
# Modified for idoc-mdeoc downloaded images - 94 AIA key is at starting of the name
#----------------------------------------------


import sunpy.data.sample
import sunpy.map
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
from datetime import datetime, timedelta
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from sunpy.coordinates import Helioprojective
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import astropy.time as time
import pathlib
import shutil
from itertools import groupby
#from astropy import time
#import time


xrt_fnames=[]
aia_fnames=[]
xrt_dates=[]
aia_all_dates=[]
aia_dates=[]
#listing all XRT-Be thin Images

#'''
aia_path='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Oct_22/source/' #Source: path to where AIA cut outs are downloaded.

for root, d, f_names in os.walk(aia_path):
    for f in f_names:
        if f.endswith(".fits"):
            aia_fnames.append(os.path.join(root, f))  # file name along with path
            #date_obj=datetime.strptime(fits.open(os.path.join(root, f))[1].header['DATE-OBS'],'%Y-%m-%dT%H:%M:%S.%f')
    
 
          #aia_all_dates.append(str(date_obj))
ref_fnames=[]
for root, d, f_names in os.walk(aia_path):
    for f in f_names:
        if f.startswith("aia.lev1.94"):
            ref_fnames.append(os.path.join(root, f))


''' 
# To check missing file and re download

start_time = time.Time('2012-01-23T03:12:48.620', scale='utc', format='isot')
t=start_time
end_time=time.Time('2012-01-23T04:00:00.344', scale='utc', format='isot')
#T_start='2012-01-23T03:12:48.620','2012-01-23T03:13:00.635','2012-01-23T03:13:12.626' , '2012-01-23T03:13:24.621' 
print(end_time)
Time_bin=[]
while t < end_time:
    Time_bin.append(t)
    t=t+timedelta(seconds=12)
    #print(t)

print(len(Time_bin))

#for i in range(len(Time_bin)-1):
    #print(Time_bin[i],'|',Time_bin[i],Time_bin[i+1])
    #pathlib.Path(aia_path+Time_bin[i]).mkdir(parents=True, exist_ok=True)
print(len(aia_fnames))
#test_list=[]
#res = [list(i) for j, i in groupby(test_list,lambda a: a.split('_')[0])]
'''

aia_fnames.sort()
Length=int((len(ref_fnames)))
print(Length)
folder_list=[]
dest_fold='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Oct_22/AIA_data/' #Destination loaction, Files will be moved to this folder
for j in range(Length):
    Fname=(os.path.splitext((ref_fnames[j].split(os.sep))[-1]))[0]
    fld_nm=(Fname[23:39])
    pathlib.Path(dest_fold+fld_nm).mkdir(parents=True, exist_ok=True)
    
    t1=j*7
    t2=t1+7
    #print(fld_nm)
    
    folder_list.append(dest_fold+fld_nm)
    #aia.lev1.131.1798872856.2024-10-22T00:00:54.623Z.image_lev1.fits
    for k in range(len(aia_fnames)):
        f=(os.path.splitext((aia_fnames[k].split(os.sep))[-1]))[0]
        if f[8:10]=='94':
            print(f)
        fl=((os.path.splitext((aia_fnames[k].split(os.sep))[-1]))[0])[24:40]
        
        if fld_nm==fl:
            #print('can be moved')
            #print('moving',aia_fnames[k], 'to',fld_nm)
            shutil.move(aia_fnames[k],(dest_fold+fld_nm))
        
np.savetxt(dest_fold+'List_flare_frames.dat',folder_list,fmt='%s')
    #print(aia_fnames[t1:t2])
        


