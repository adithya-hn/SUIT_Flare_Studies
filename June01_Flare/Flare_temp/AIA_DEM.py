#-----------------------------------------------------------------
#Author: Prakhar
#Purpose: Calculate DEM and generate temperature Map from AIA Images
#Method: Using established Hannah & Kontar (2012) 
#Data of creation: 8/12/2023

#Modification History (Adithya)
# - Optimised and added the folder reading format for batch mode
# - Added Tempetrature error file
# - Corrected the Tempareture error data shape in demreg/dem2pos.py
# - Added Modified headers for Temperature and error fits file
#
#-----------------------------------------------------------------

#--Package importing---
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import scipy.io as io
from demregpy import dn2dem
import demregpy
import glob
from sunpy.net import Fido, attrs as a
from astropy.visualization import wcsaxes
from sys import path as sys_path
import astropy.time as atime
from astropy.coordinates import SkyCoord
from astropy import units as u
import sunpy.map
import sys
from sunpy.net import Fido, attrs
from sunpy.map import Map
import pathlib
from astropy import units as u, time as time
from astropy.io import fits
from aiapy.calibrate import degradation
from aiapy.calibrate.util import get_correction_table
from aiapy.calibrate import register, update_pointing
from aiapy.calibrate import degradation, register, update_pointing, correct_degradation

import os
import warnings
import timeit
warnings.simplefilter('ignore')

# Change to your local copy's location...
sys_path.append('/Analysis/Research_Projects/AIA_temp/DEM_Package/demreg-master/python/')
from dn2dem_pos import dn2dem_pos

#--------$$$-----------

#Required files

trin=io.readsav('/Analysis/Research_Projects/AIA_temp/DEM_Package/aia_tresp_en.dat')#AIA Response file
fl_list = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June01_Flare/Flare_temp/List_flare_frames.dat', dtype='str').transpose()

xrt_fls=fl_list
xrt_fls.sort()
Length = len(xrt_fls)
pathlib.Path("Temp_Fits").mkdir(parents=True, exist_ok=True) #Make folders
pathlib.Path("Temp_maps").mkdir(parents=True, exist_ok=True)
pathlib.Path("TempEr_Fits").mkdir(parents=True, exist_ok=True)

#print(1001,xrt_fls)
tot_n_fls=len(xrt_fls)
print('Total images set',len(xrt_fls))
rej=[]

startTime = timeit.default_timer()
prvTime=startTime
totelIm=0

for k in range(tot_n_fls):
    try:
        image=(os.path.split(xrt_fls[k]))[1]
        print('Processing image:',image)
        f_name=image

        fdir = xrt_fls[k]+'/'#/home/adithya/Desktop/Flare_temp/AIA_ROIs/ROI_Sample4_full_res/' +image+'/'
        #print(fdir)
        
        Temp_Bins=31

        #--------$$$-------------

        files = sorted(glob.glob(fdir + '*.fits'))
        #print(files)
        if len(files) ==0:
            print('No AIA files in folder')
            continue

        if len(files) <6:
            print('Insufficient files')
            continue

        if len(files) >6:
            print('More than 6 images')
            continue

        maps = [sunpy.map.Map(file) for file in files]
        maps = [Map(file) for file in files]
        maps = sorted(maps, key=lambda x: x.wavelength)
        #maps = [aiaprep(m) for m in maps]
        maps = [correct_degradation(m)/m.exposure_time for m in maps]
        

        wvn0 = [m.meta['wavelnth'] for m in maps]
        srt_id = sorted(range(len(wvn0)), key=wvn0.__getitem__)

        maps = [maps[i] for i in srt_id]
        #print([m.meta['wavelnth'] for m in maps])
    
        wvn = [m.meta['wavelnth'] for m in maps]
        worder=np.argsort(wvn)

        durs = [m.meta['exptime'] for m in maps]

        durs=np.array(durs)
        wvn=np.array(wvn)

        #print('duration -',durs)
        

        # Get the temperature response functions in the correct form for demreg
        tresp_logt=np.array(trin['logt'])
        nt=len(tresp_logt)
        nf=len(trin['tr'][:])
        trmatrix=np.zeros((nt,nf))
        #print('trmtx shape',trmatrix.shape)
        for i in range(0,nf):
            trmatrix[:,i]=trin['tr'][i]

        

        #print(nt,nf)
        temps=np.logspace(5.7,7.0,num=Temp_Bins)
        # Temperature bin mid-points for DEM plotting
        mlogt=([np.mean([(np.log10(temps[i])),np.log10((temps[i+1]))]) \
                for i in np.arange(0,len(temps)-1)])

        mtemps=([np.mean([(temps[i]),(temps[i+1])]) \
                    for i in np.arange(0,len(temps)-1)])
        log_temps = np.log10(temps)

        nx = int(maps[0].dimensions.x.value)
        ny = int(maps[0].dimensions.y.value)
    
        nf = len(files)
        data = np.zeros([ny, nx, nf])
        #print(data.shape)
        #convert from our list to an array of data
        for j in np.arange(nf):
            #print('shape',(maps[j].data).shape)
            data[:,:,j]=maps[j].data
            
        data[data < 0]=0
        serr_per=10.0
       
        npix=1 #nx*ny
        edata=np.zeros([ny,nx,nf])
        gains=np.array([18.3,17.6,17.7,18.3,18.3,17.6])
        dn2ph=gains*[94,131,171,193,211,335]/3397.0
        rdnse=1.15*np.sqrt(npix)/npix
        drknse=0.17
        qntnse=0.288819*np.sqrt(npix)/npix
        #print('---1')
        
        for j in np.arange(nf):
            etemp=np.sqrt(rdnse**2.+drknse**2.+qntnse**2.+(dn2ph[j]*abs(data[:,:,j]))/(npix*dn2ph[j]**2))
            esys=serr_per*data[:,:,j]/100.
            edata[:,:,j]=np.sqrt(etemp**2. + esys**2.)
        
       
        
        #print('Calculating DEM:')
        
        dem0, edem0, elogt0, chisq0, dn_reg0 = dn2dem_pos(data, edata, trmatrix, tresp_logt, temps,max_iter=100, rgt_fact=1.5)
        

        em0 = np.zeros_like(dem0)
        em0_temp= np.zeros_like(dem0)
        ertem1=np.zeros_like(edem0)
        erDem=np.zeros_like(edem0)
        

        for j in range(0,30):
            em0[:, :, j] = dem0[:, :, j] *(temps[j + 1] - temps[j])
            em0_temp[:, :, j] = em0[:, :, j]*(mtemps[j])
            ertem1[:,:,j]=pow((mtemps[j]),2)*pow((temps[j + 1] - temps[j]),2)*pow(edem0[:,:,j],2)
            erDem[:,:,j]=pow(edem0[:,:,j],2)

                
        total_em = np.sum(em0[:, :,0:30], axis=2)
        total_weighted_em = np.sum(em0_temp[:, :,0:30],axis=2)
        mean_temp = (total_weighted_em)/(total_em)

        Header_part=maps[1].fits_header
        

        hdu=fits.PrimaryHDU(mean_temp)
     
        hdu.header['DATE-OBS']=maps[0].meta.get('DATE-OBS')
        hdu.header['CRPIX1']=maps[0].meta.get('CRPIX1')
        hdu.header['CRPIX2']=maps[0].meta.get('CRPIX2')
        hdu.header['CRVAL1']=maps[0].meta.get('CRVAL1')
        hdu.header['CRVAL2']=maps[0].meta.get('CRVAL2')
        hdu.header['CTYPE1']=maps[0].meta.get('CTYPE1')
        hdu.header['CTYPE2']=maps[0].meta.get('CTYPE2')
        hdu.header['CUNIT1']=maps[0].meta.get('CUNIT1')
        hdu.header['CUNIT2']=maps[0].meta.get('CUNIT2')
        hdu.header['CDELT1']=maps[0].meta.get('CDELT1')
        hdu.header['CDELT2']=maps[0].meta.get('CDELT2')
        hdu.header['CROTA2']=maps[0].meta.get('CROTA2')
        hdu.header['DSUN_REF']=maps[0].meta.get('DSUN_REF')
        hdu.header['DSUN_OBS']=maps[0].meta.get('DSUN_OBS')
        hdu.header['RSUN_REF']=maps[0].meta.get('RSUN_REF')
        hdu.header['RSUN_OBS']=maps[0].meta.get('RSUN_OBS')
        hdu.header['CRLN_OBS']=maps[0].meta.get('CRLN_OBS')
        hdu.header['CRLT_OBS']=maps[0].meta.get('CRLT_OBS')
        hdu.header['HGLN_OBS']=maps[0].meta.get('HGLN_OBS')
        hdu.header['HGLT_OBS']=maps[0].meta.get('HGLT_OBS')
        hdu.header['HGLN_OBS']=maps[0].meta.get('HGLN_OBS')
        

        
        hdu.writeto('Temp_Fits/Temp_{}.fits'.format(f_name),overwrite=True)

        plt.imshow(mean_temp, 'jet',origin='lower',vmin=1000000,vmax=8000000)
        #plt.imshow(np.log10(mean_temp), 'hinodexrt', origin='lower')
        plt.title('Mean Temperature Map')
        plt.axis('off')
        plt.colorbar()
        plt.savefig('Temp_maps/Temp_{}.png'.format(f_name), dpi=1000, bbox_inches='tight')
        #plt.show(block=False)
        #plt.pause(1)
        plt.close()
        #plt.imshow()

        

        
        #error estimation
        ErdmSmsq=np.sum((pow(edem0[:,:,0:30],2)),axis=2)
        ErdmSm=np.sqrt(np.sum((erDem[:,:,0:30]),axis=2))
        demSmsq=pow(np.sum(em0[:, :,0:30], axis=2),2)
        Tsq=pow(mean_temp,2)
        erT=10**elogt0
        #print(demSmsq.shape,Tsq.shape,ErdmSmsq.shape,erT.shape,em0.shape)
        demTsq=np.sum((pow((em0*erT),2))[:,:,0:30],axis=2)

        #print(demSmsq.shape,Tsq.shape,ErdmSmsq.shape)
        #print((np.sum(ertem1[:,:,0:30],axis=2).shape))
        ErMtem=np.sqrt((np.sum(ertem1[:,:,0:30], axis=2)/demSmsq)+(Tsq/demSmsq)*ErdmSmsq+demTsq/demSmsq)
        
        hdu1=fits.PrimaryHDU(ErMtem)
        hdu1.header=Header_part
        hdu1.writeto('TempEr_Fits/Temp_Er_{}.fits'.format(f_name),overwrite=True)
        
        tempstopTime = timeit.default_timer()
        filetime = np.round((tempstopTime - prvTime),2)
        prvTime=tempstopTime
        print('[ ',k,' /',tot_n_fls,' ] ',f_name,'Time:',round((filetime/60),1),'Min') 


    except:
        print('skipped the files')
        print('error is:', sys.exc_info()[0])
        rej.append(image)

np.savetxt('Skipped_files.dat',rej,fmt='%s')
