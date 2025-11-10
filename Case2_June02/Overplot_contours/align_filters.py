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
# suit_filters=['NB01','NB02','NB05','NB06','NB07','NB09']
suit_filters=['NB03','NB04']
suit_raw_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case2_June02/data/raw/'
#ref_fd_img_pth='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case2_June02/data/aia/hmi.ic_45s.20240602_022915_TAI.2.continuum.fits'
ref_fd_img_pth='/media/adithya/Adi_disk4/SUIT_flare_work/case2_Jun02/data/aia/aia_fd/aia.lev1_uv_24s.2024-06-02T022952Z.1600.image_lev1.fits'
tx1,ty1=-330,-350 
tx2,ty2=-200,-230

shift_x=5
shift_y=5

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