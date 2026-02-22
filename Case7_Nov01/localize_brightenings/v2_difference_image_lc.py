
'''
Created on 1 Oct 2025
@author: adithya-hn

- Running code for difference image
'''

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
#import suit_difference_image
import v2_diff_img

p1,p2,p3,p4=40,670,40,680
flt='NB04' #NB04,NB08
data_folder=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/aligned_crop_fits/*{flt}.fits'
clip_bar=False
get_histograms=False
v2_diff_img.diff_img(data_folder,flt,thresh_sig=5)

