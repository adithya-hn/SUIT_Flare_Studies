
import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
from datetime import datetime
from datetime import timedelta
import timeit
import pathlib
from astropy.time import Time
import numpy as np
from sunpy.time import parse_time

start = timeit.default_timer()
now = datetime.now()-timedelta(days=1)
fol_nm='/Analysis/Projects_Data/Flare_Data/Oct09_Flare_Data/'
jpg_fold=fol_nm+'/'+'P_corr_data'

search_fold='/Analysis/Projects_Data/Flare_Data/Oct09_Flare_Data/'
fdir =search_fold 
pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)

files = sorted(glob.glob(fdir + '*.fits'))


P_angle=-3.23
for i in range (len(files)):
    suitMap=sunpy.map.Map(files[i])
    obsDate=suitMap.date
    ref_time=Time(datetime.fromisoformat(suitMap.meta.get('CRTIME')))
    cord_dif=obsDate-ref_time

    print(abs(suitMap.meta.get('CRVAL1')/0.698),suitMap.meta.get('CRVAL1'))
#NAXIS1
    suitMap.meta['CRPIX1']=-1*(suitMap.meta.get('CRVAL1')/0.698)+(suitMap.meta.get('NAXIS1')/2)
    suitMap.meta['CRPIX2']=-1*(suitMap.meta.get('CRVAL2')/0.698)+(suitMap.meta.get('NAXIS1')/2)
    suitMap.meta['CRVAL1']=0
    suitMap.meta['CRVAL2']=0
    suitMap.meta['CROTA2']=P_angle
    #print(crpix1,crpix2)
    suitMap.save((jpg_fold+'/'+os.path.basename(files[i])),overwrite=True)


    '''
    idx2=np.where(ref_time==Pdate)[0]
    if len(idx2)==0:
        print('not found finding the closest')
        idx2=np.argmin(np.abs(Pdate - ref_time))
    else:
         idx2=idx2[0]
    print(ref_time,obsDate,Pdate[idx2])'''



    


