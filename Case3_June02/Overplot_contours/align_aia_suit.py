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
suit_raw_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case3_June02/data/raw/'
aia_imgs_pth='/media/adithya/Adi_disk4/SUIT_flare_work/case3_Jun02/data/aia/cut_outs/'
hmi_imgs_pth='/media/adithya/Adi_disk4/SUIT_flare_work/case3_Jun02/data/hmi/HMI_cutouts/'
ref_fd_1600_pth='/media/adithya/Adi_disk4/SUIT_flare_work/case3_Jun02/data/aia/aia_fd_data/aia.lev1_uv_24s.2024-06-02T062952Z.1600.image_lev1.fits'

#Template cutout
tx1,ty1=-330,-350 
tx2,ty2=-180,-230


alin_fltr='NB04'   #Filter to align other SUIT filters 
save_aligned_fits='yes'
save_pngs='no'     #aligned pngs
draw_contours='yes'
fol_nm=os.getcwd() #Custom folder to save contour images

co_align_suit_aia.start_coaligning(Filters,suit_raw_files,aia_imgs_pth,hmi_imgs_pth,ref_fd_1600_pth,tx1,tx2,ty1,ty2,alin_fltr,fol_nm,draw_contours,save_aligned_fits,save_pngs)
