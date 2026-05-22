
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
fol_nm='/Analysis/Projects_Data/Flare_Data/Nov01_Flare_Data2/'
save_fold=fol_nm+'/'+'P_corr_data'

search_fold=fol_nm
fdir =search_fold 
pathlib.Path(save_fold).mkdir(parents=True, exist_ok=True)

files = sorted(glob.glob(fdir + '*.fits'))
print('No. of images found: ',len(files))

P_angle=6.27

for i in range (len(files)):
    suitMap=sunpy.map.Map(files[i])
    obsDate=suitMap.date
    ref_time=Time(datetime.fromisoformat(suitMap.meta.get('CRTIME')))
    cord_dif=obsDate-ref_time # to see coordinates assigned from are correct ref img or not

    #print(abs(suitMap.meta.get('CRVAL1')/0.698),suitMap.meta.get('CRVAL1'))

    suitMap.meta['CRPIX1']=-1*(suitMap.meta.get('CRVAL1')/0.698)+(suitMap.meta.get('NAXIS1')/2)
    suitMap.meta['CRPIX2']=-1*(suitMap.meta.get('CRVAL2')/0.698)+(suitMap.meta.get('NAXIS1')/2)
    suitMap.meta['CRVAL1']=0
    suitMap.meta['CRVAL2']=0
    suitMap.meta['CROTA2']=P_angle
    #print(crpix1,crpix2)
    suitMap.save((save_fold+'/'+os.path.basename(files[i])),overwrite=True)




    


