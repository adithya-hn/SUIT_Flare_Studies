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



suit_raw_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case3_June02/data/raw/'

# Filters=['6173']
# suit_filters=['NB01','NB02','NB05','NB06','NB07','NB09']
# ref_fd_img_pth='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case3_June02/data/aia/hmi.ic_45s.20240602_062915_TAI.2.continuum.fits'

Filters=['1600']
suit_filters=['NB04','NB03','NB08']
ref_fd_img_pth='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case3_June02/data/aia/aia.lev1_uv_24s.2024-06-02T062952Z.1600.image_lev1.fits'

#Template cutout
tx1,ty1=-330,-350 
tx2,ty2=-180,-230

save_aligned_fits='yes'
save_pngs='no'     #aligned pngs
draw_contours='yes'
fol_nm=os.getcwd() #Custom folder to save contour images

jpg_fold=fol_nm+'/'+'Contour_imgs'
"""nb1_fl = glob.glob(suit_raw_files + '*3'+'NB01.fits')
nb2_fl = glob.glob(suit_raw_files + '*3'+'NB02.fits')
nb5_fl = glob.glob(suit_raw_files + '*3'+'NB05.fits')
nb6_fl = glob.glob(suit_raw_files + '*3'+'NB06.fits')
nb7_fl = glob.glob(suit_raw_files + '*3'+'NB07.fits')

print('Total SUIT NB01 files:',len(nb1_fl))
print('Total SUIT NB02 files:',len(nb2_fl))
print('Total SUIT NB05 files:',len(nb5_fl))
print('Total SUIT NB06 files:',len(nb6_fl))
print('Total SUIT NB07 files:',len(nb7_fl))

print('---------------')"""


for fltr in suit_filters:
    print(f'Starting with {fltr}')
    fltr_fl = glob.glob(suit_raw_files + '*3'+f'{fltr}.fits')
    fltr_fl=sorted(fltr_fl, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    print(f'Total SUIT {fltr} files:',len(fltr_fl))
    if len(fltr_fl)==0:
        print(f'No {fltr} filter files')
        continue
    align_suit_fltr_to_aia.co_align_maps(suit_raw_files,jpg_fold,ref_fd_img_pth,fltr,fltr_fl,tx1,tx2,ty1,ty2,save_pngs=None,save_fits='yes',draw_contours='yes',cor_x=5,cor_y=5)
