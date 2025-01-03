
import os
import timeit
import numpy as np
import shutil
import matplotlib.pyplot as plt
import drms
import sunpy.map
from sunpy.net import attrs as a  
from datetime import datetime, timedelta
import requests
import re
from sunpy.net import Fido
from astropy import units as u
import sunpy.net.jsoc as jsoc
from sunpy.net import attrs as a
from itertools import repeat
import pandas as pd
from tqdm import tqdm
import requests
from bs4 import BeautifulSoup
import platform
from tqdm import tqdm

client = drms.Client()

save_dir = '/home/adithya/Research_Projects/Flare_studies/SUIT_Flares/May19_Flare/SHARP_Data/'
email_id = 'adithyabhattsringeri@gmail.com'

data = (pd.read_csv('/home/adithya/Research_Projects/Flare_studies/SUIT_Flares/May19_Flare/Suit_flare_Sharp_downloadList.csv', sep=',')).transpose()
Data=data.values
Ts=Data[2]
novaa_num=Data[1]
Harp_num=Data[0]

tot_files=len(Harp_num)
print('Total files to download',tot_files)

pars=['{magnetogram}','{Bp}','{Bt}','{Br}','{conf_disambig}','{bitmap}','{Bt_err}','{Bp_err}','{Br_err}']
loop_idx=0
for key in range(len(Ts)):

    time = str(Ts[key]) #Redundant line added for Uniformit
    noaa = novaa_num[key] 
    harp = Harp_num[key]
    loop_idx+=1

    #print(harp,noaa,time)
    if harp:
        now = datetime.now()
        
        print(f"[ {loop_idx} / {tot_files} ]Downloading files for NOAA: {noaa} | HARP number: {harp} \t [{now}]")
        
        noaa_dir = os.path.join(save_dir,str(noaa)) 

        if not os.path.exists(noaa_dir):
            os.mkdir(noaa_dir)
        

        
        ts_str = time.replace(" ","_").replace("-","").replace(":","")
        print('time:',ts_str)
        ts_dir = os.path.join(noaa_dir,ts_str)
        #print(ts_dir)

        if not os.path.exists(ts_dir):
            os.mkdir(ts_dir)

        if len(os.listdir(ts_dir)) < 9 :    

            qry_date=time[:4]+'.'+time[5:7]+'.'+time[8:10]+'_'+time[11:]+'_TAI'
            for parameters in tqdm(pars):
                query_string = 'hmi.sharp_cea_720s'+'['+str(harp)+']'+'['+qry_date+']' +parameters
                #print(query_string)

                export_request = client.export(query_string, method='url', protocol='fits',email=email_id)
                export_request.wait()  
                export_request.status 
                export_request.request_url  
                export_request.download(ts_dir)

        else:
            print(f"Files for {noaa} at {time} already exist ") 

        
        now = datetime.now()
        print(f"   Downloaded all FITS files.  \t\t {now}")



print(f"Task Completed")

#Remove From here if running over SSH or if Beep annoys you. ;)  
