'''
Created 25 Oct 2025
Author: adithya-hn

Purpose: running version of co_alin_suit_aia

'''

import os
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
import co_align_suit_aia

Filters=['1600']
suit_raw_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/raw/'
aia_imgs_pth='/media/adithya/Adi_disk4/SUIT_flare_work/case1_Jun01/data/aia/cut_outs/'
hmi_imgs_pth='/media/adithya/Adi_disk4/SUIT_flare_work/case1_Jun01/data/hmi/HMI_cutouts/'
suit_filters=['NB03','NB08','NB04']
ref_fd_1600_pth='/media/adithya/Adi_disk4/SUIT_flare_work/case1_Jun01/data/aia/aia_fd/aia.lev1_uv_24s.2024-06-01T070952Z.1600.image_lev1.fits'
tx1,ty1=-475,-330
tx2,ty2=-350,-230

alin_fltr='NB04'
save_aligned='yes'
save_pngs='no'
draw_contours='yes'
fol_nm=os.getcwd() #Custom folder to save contour images

co_align_suit_aia.start_coaligning(Filters,suit_raw_files,aia_imgs_pth,hmi_imgs_pth,ref_fd_1600_pth,tx1,tx2,ty1,ty2,alin_fltr,fol_nm,draw_contours,save_aligned_fits,save_pngs)
