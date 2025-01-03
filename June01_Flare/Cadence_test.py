
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



start = timeit.default_timer()
now = datetime.datetime.now()-timedelta(days=1)
fol_nm='/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Processed_Data'

jpg_fold=fol_nm+'/'+'Coloured_ROI_imgs'
algn_dir=fol_nm+'/'+'Aligned_images'

#Filter='NB03'
Filters=['NB03','NB04','NB08']#,'BB01','BB02','BB03']
pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
search_fold='/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Flare_data/'
fdir =search_fold 


for fltr in Filters:
    pathlib.Path(jpg_fold+'/'+fltr).mkdir(parents=True, exist_ok=True)
    pathlib.Path(algn_dir+'/'+fltr).mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(fdir + '*3'+fltr+'.fits'))
    print('Total ',fltr ,' files:',len(files))
    for i in range (len(files)):
        img_map=sunpy.map.Map(files[i])
        obs_time=img_map.date
        #print((sunpy.time.parse_time(obs_time).datetime))
        if i >0:
            t1=sunpy.map.Map(files[i]).date
            t2=sunpy.map.Map(files[i-1]).date
            T1=(sunpy.time.parse_time(t1).datetime).timestamp()
            T2=(sunpy.time.parse_time(t2).datetime).timestamp()
            print(T1-T2)
        