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
from astropy.io import fits
import pathlib


start = timeit.default_timer()
now = datetime.datetime.now()-timedelta(days=1)
fol_nm='/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/ReShaped'
jpg_fold=fol_nm+'/'+'Coloured_ROI_imgs'
algn_dir=fol_nm+'/'+'Aligned_images'

#Filter='NB03'
Filters=['NB03','NB04','NB08']
pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
search_fold='/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Flare_data/' #'/scratch/suit_data/level1fits/2024/'+str(now.month).zfill(2)+'/'+str(now.day).zfill(2)+'/'+'normal_2k'+'/'
fdir =search_fold 


for fltr in Filters:
    pathlib.Path(jpg_fold+'/'+fltr).mkdir(parents=True, exist_ok=True)
    pathlib.Path(algn_dir+'/'+fltr).mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(fdir + '*3'+fltr+'.fits'))
    print('Total ',fltr ,' files:',len(files))
    aln_imgs=[]
    #Sequence = sunpy.map.Map(files, sequence=True)
    for i in range(len(files)):
        Img=fits.open(files[i])
        Data=Img[0].data
        print(Data.shape,os.path.basename(files[i]))
        if Data.shape[0]==704:
            Data=Img[0].data[:640,64:]
            Head=Img[0].header
            fld=algn_dir+'/'+fltr
            fnm=os.path.join(fld,os.path.basename(files[i]))
            print(fnm)
            hdu=fits.PrimaryHDU(Data)
            hdu.header=Head
            hdu.writeto(fnm,overwrite=True)


