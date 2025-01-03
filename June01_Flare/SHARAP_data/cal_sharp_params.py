'''
Purpose: Create Currenrt Helicity, Unsigned Current and Free energy image fits
Method: Calling the modified Monica Bobara's calculate sharpkeys 
Source Code: Monica Bobra @mbobra
Author: Adithya

'''

import numpy as np
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import math
import sympy as sp
from tqdm import tqdm
import os
import sunpy.map
from skimage.measure import block_reduce
import scipy.ndimage
import imageio
from skimage import exposure
from datetime import datetime, timedelta
import json

import warnings
warnings.simplefilter('ignore')

import mod_calculate_sharpkeys as calculate_sharpkeys



base_dir = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June01_Flare/SHARAP_data/' # where all downloaded data available


radsindeg = np.pi/180.
munaught  = 0.0000012566370614

noaas =['13697'] # os.listdir(base_dir)
#noaas.sort()
tot_files=len(noaas)
loop_idx=0
print('Total Active regions folder: ',tot_files)
skip_list=[]

wrk_dir=os.getcwd()
with open("skip_files.txt", "a") as f: #clear previous entries
	f.truncate(0)

def update_skipFile(file,miss_par): # accounts for skipped file
	with open("skip_files.txt", "a") as f:

		f.write( file+' '+miss_par+"\n")

file_pth = os.path.join(wrk_dir,'skip_files.txt') 
#print(file_pth)

for noaa in noaas:
	loop_idx+=1
	if noaa.isdigit():

		times = os.listdir(os.path.join(base_dir,noaa))
		times.sort()
		tot_time_stamps=len(times)
		print()
		print(f"Calculating parameters for NOAA {noaa} [ {loop_idx} / {tot_files} ]")
			

		for tstamp in tqdm(times):

			ts_dir = os.path.join(os.path.join(base_dir,noaa),tstamp)
			
			if os.path.isdir(ts_dir):
				files = os.listdir(ts_dir)
				files.sort()

				file_by=file_bx=file_bz=file_by_err=file_bx_err=file_bz_err=file_los=file_bitmap=file_conf_disambig = 0
				f_count=0

				for f in files:
					f = os.path.join(ts_dir,f)
					if f.endswith('Bt.fits'):
						file_by = f
					elif f.endswith('Bp.fits'):
						file_bx = f
					elif f.endswith('Br.fits'):
						file_bz = f
					elif f.endswith('conf_disambig.fits'):
						file_conf_disambig = f
					elif f.endswith('bitmap.fits'):
						file_bitmap = f
					elif f.endswith('Bt_err.fits'):
						file_by_err = f
					elif f.endswith('Bp_err.fits'):
						file_bx_err = f
					elif f.endswith('Br_err.fits'):
						file_bz_err = f
					elif f.endswith('magnetogram.fits'):
						file_los = f

					elif f.endswith('.fits.1'):
						os.remove(f)
						
					elif f.endswith('calculated_params.json'):
						
						f_count=1

						break

				#print(inc,mgm,fld_err,azi_err,inc_err,bmp,cds,fld,azi)
				if f_count ==1:
					print('Already files exist')
					continue

				if file_by ==0:
					print('Missing B theta file')
					print('Skipping it')
					skip_list.append(ts_dir)
					update_skipFile(ts_dir,'B theta')
					continue
				if file_bz ==0:
					print('Missing B r file')
					print('Skipping it')
					skip_list.append(ts_dir)
					update_skipFile(ts_dir,' Br')
					continue

				if  file_bx ==0:
					print('Missing B phi file')
					print('Skipping it')
					skip_list.append(ts_dir)
					update_skipFile(ts_dir,'Bp')
					continue

				if file_by_err ==0:
					print('Missing Bt error file')
					print('Skipping it')
					skip_list.append(ts_dir)
					update_skipFile(ts_dir,'Bt error')
					continue

				if file_bz_err ==0:
					print('Missing Br error file')
					print('Skipping it')
					skip_list.append(ts_dir)
					update_skipFile(ts_dir,'Br error')
					continue

				if file_bx_err ==0:
					print('Missing B phi file')
					print('Skipping it')
					skip_list.append(ts_dir)
					update_skipFile(ts_dir,'Bphi')
					continue

				if file_conf_disambig ==0:
					print('Missing conf disambig file')
					print('Skipping it')
					skip_list.append(ts_dir)
					update_skipFile(ts_dir,'conf disambig')
					continue

				if file_los ==0:
					print('Missing mag file')
					print('Skipping it')
					skip_list.append(ts_dir)
					update_skipFile(ts_dir,'mag')
					continue

				if file_bitmap ==0:
					print('Missing bitmap file')
					print('Skipping it')
					skip_list.append(ts_dir)
					update_skipFile(ts_dir,'bitmap')
					continue
				
				#print('Getting the data.')
				bz, by, bx, bz_err, by_err, bx_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, los, los_err = calculate_sharpkeys.get_data(file_bz, file_by, file_bx, file_bz_err, file_by_err, file_bx_err, file_conf_disambig, file_bitmap, file_los)
		
                # Creating common header file for all fits file
				with fits.open(file_bz) as hdul_bz:
					#print('Bz file:', os.path.basename(file_bz))
					bz = hdul_bz[0].data if hdul_bz[0].data is not None else hdul_bz[1].data
					bz_header = hdul_bz[1].header if hdul_bz[1].header is not None else hdul_bz[0].header

					#--For proper Header information--
					req_header_keys=['DATE-OBS','T_REC ','T_OBS ','CTYPE1','CTYPE2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','CDELT1','CDELT2','CUNIT1','CUNIT2','IMCRPIX1','IMCRPIX2','IMCRVAL1','IMCRVAL2','CROTA2','WCSNAME','DSUN_OBS','DSUN_REF','RSUN_REF','CRLN_OBS','CRLT_OBS','CAR_ROT','RSUN_OBS','LAT_MIN','LON_MIN','LAT_MAX','LAT_MAX','LON_MAX','N_PATCH','N_PATCH1','N_PATCHM','HARPNUM','NOAA_AR','NOAA_NUM','NOAA_ARS',]

					emt_ar=np.empty_like(bz[1].data)
					empty_fits=fits.PrimaryHDU(emt_ar)
					for  key_ in req_header_keys:
						empty_fits.header[key_]=bz_header[key_]
					
					file_headers=empty_fits.header # this header you can use for all files
				
				#---------'''------#
				
				
				# compute the vertical current and associated errors
				current                =calculate_sharpkeys.computeJz(bx, by, bx_err, by_err, conf_disambig, bitmap, nx, ny)
				jz, jz_err, derx, dery = current[0], current[1], current[2], current[3]
				# compute the moments of the current helicity and associated errors 
				mean_ih, mean_ih_err, total_us_ih, total_us_ih_err, total_abs_ih, total_abs_ih_err,US_JH_img,USJH_err_img=calculate_sharpkeys.computeHelicity(jz, jz_err, bz, bz_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec)
				potential =calculate_sharpkeys.greenpot(bz, nx, ny)
				bpx, bpy  = potential[0], potential[1]
				# compute the energy stored in the magnetic field and its associated errors
				meanpot, meanpot_err, totpot, totpot_err,tot_pot_img,tot_pot_img_er=calculate_sharpkeys.computeFreeEnergy(bx_err, by_err, bx, by, bpx, bpy, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, conf_disambig, bitmap)
				# compute the moments of the vertical current density and associated errors
				mean_jz, mean_jz_err, us_i, us_i_err,USJz_Img,USJz_er_Img = calculate_sharpkeys.computeJzmoments(jz, jz_err, derx, dery, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, munaught)
				Rparam, Rparam_err,pmap = calculate_sharpkeys.computeR(los, los_err, nx, ny, cdelt1_arcsec)

				print('R_VALUE ', Rparam,'Mx')
				print('The error in R_VALUE is', Rparam_err)
				
				# compute the degree to which the observed field is sheared and its associated errors
				meanshear_angle, meanshear_angle_err, area_w_shear_gt_45,shear_img,shear_img_er = calculate_sharpkeys.computeShearAngle(bx_err, by_err, bz_err, bx, by, bz, bpx, bpy, nx, ny, conf_disambig, bitmap)
				print('MEANSHR ',meanshear_angle,'degree')
				print('ERRMSHA ',meanshear_angle_err,'degree')
				print('SHRGT45 ',area_w_shear_gt_45,'as a percentage')
				
				'''---Verifying---
				print('Total Unsigned Current Helicity: ',abs(US_JH_img[US_JH_img!=0]).sum(),total_us_ih)
				print('Total Free energy: ',tot_pot_img[tot_pot_img!=0].sum(),totpot)
				print('Total Unsigned Current: ',USJz_Img[USJz_Img!=0].sum(),us_i)
				'''

				f_pattern = str(file_bz.replace("Br.fits",""))
				jz_fits = f_pattern+'jz'+'.fits'
				jz_err_fits = f_pattern+'jz_err'+'.fits'

				us_fits = f_pattern+'USJH'+'.fits'
				bh_err_fits = f_pattern+'USJH_err'+'.fits'

				JH_hdu = fits.PrimaryHDU(US_JH_img, header=file_headers)
				JH_hdu.writeto(us_fits, overwrite=True)

				JH_er_hdu = fits.PrimaryHDU(USJH_err_img, header=file_headers)
				JH_er_hdu.writeto(bh_err_fits, overwrite=True)

				totpot_fits = f_pattern+'TOTPOT'+'.fits'
				totpot_err_fits = f_pattern+'TOTPOT_err'+'.fits'

				TOTPOT_hdu = fits.PrimaryHDU(tot_pot_img, header=file_headers)
				TOTPOT_hdu.writeto(totpot_fits, overwrite=True)

				TOTPOT_er_hdu = fits.PrimaryHDU(tot_pot_img_er, header=file_headers)
				TOTPOT_er_hdu.writeto(totpot_err_fits, overwrite=True)

				shear_fits = f_pattern+'Shear'+'.fits'
				shear_err_fits = f_pattern+'Shear_err'+'.fits'

				Shear_hdu = fits.PrimaryHDU(shear_img, header=file_headers)
				Shear_hdu.writeto(shear_fits, overwrite=True)

				Shear_hdu_er = fits.PrimaryHDU(shear_img_er, header=file_headers)
				Shear_hdu_er.writeto(shear_err_fits, overwrite=True)

				pil_fits = f_pattern+'PIL'+'.fits'
				pil_err_fits = f_pattern+'PIL_err'+'.fits'

				PIL_hdu = fits.PrimaryHDU(pmap, header=file_headers)
				PIL_hdu.writeto(pil_fits, overwrite=True)


				json_data = {
				"mean_current_helicity" : mean_ih,
				"total_unsigned_current_helicity" : total_us_ih,
				"absolute_current_helicity" : total_abs_ih,
				"mean_current_helicity_err" : mean_ih_err,
				"total_unsigned_current_helicity_err" : total_us_ih_err,
				"total_unsigned_current_helicity_err" : total_abs_ih_err,
				"mean_photospheric_magnetic_free_energy" : meanpot,
				"mean_photospheric_magnetic_free_energy_error" : meanpot_err,
				"total_photospheric_magnetic_free_energy" : totpot,
				"total_photospheric_magnetic_free_energy_error": totpot_err
				}

				
				json_path = os.path.join(ts_dir,tstamp)
				json_path = os.path.join(json_path,f_pattern+"calculated_params.json")
				with open(json_path, "w") as outfile: 
					json.dump(json_data, outfile, sort_keys=False, indent=4)

			
			
			
			else:
				print(f" Some or all Files of NOAA {noaa} are Missing for timestamp {tstamp}.")

	



				


				




				
                

