'''
Created 30 Oct 2025
Author: adithya-hn

Purpose: running version of alin_suit_fltr_to_sdo

'''

import os
import glob
import datetime
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
import align_suit_fltr_to_aia


Filters=['6173']
suit_filters=['NB02','NB05','NB06','NB07']
suit_raw_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/Processed/contam_corr_data/'
#ref_fd_img_pth='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/aia/hmi.ic_45s.20240710_032915_TAI.2.continuum.fits'
ref_fd_img_pth='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/aia/hmi.ic_45s.20240710_040000_TAI.2.continuum.fits'
#Template cutout
tx1,ty1=-200,-260 
tx2,ty2=-60,-120

shift_x=0
shift_y=0

save_aligned_fits='yes'
save_pngs='no'     #aligned pngs
draw_contours='yes'
fol_nm=os.getcwd() #Custom folder to save contour images

jpg_fold=fol_nm+'/'+'Contour_imgs'


for fltr in suit_filters:
    print(f'Starting with {fltr}')
    fltr_fl = glob.glob(suit_raw_files + '*3'+f'{fltr}.fits')
    fltr_fl=sorted(fltr_fl, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    print(f'Total SUIT {fltr} files:',len(fltr_fl))
    if len(fltr_fl)==0:
        print(f'No {fltr} filter files')
        continue
    align_suit_fltr_to_aia.co_align_maps(suit_raw_files,jpg_fold,ref_fd_img_pth,fltr,fltr_fl,tx1,tx2,ty1,ty2,save_pngs=None,save_fits='yes',draw_contours='yes',cor_x=shift_x,cor_y=shift_y)
    print('-------------')