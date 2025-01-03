
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

Filters=['NB03','NB04','NB08']#,'BB01','BB02','BB03']

search_fold='/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/'
fdir =search_fold 
cadence=[]
time_stmp=[]

for fltr in Filters:

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
            cadence.append(T1-T2)
            time_stmp.append(t1)
            print(T1-T2)
plt.plot(time_stmp,cadence)
plt.show()