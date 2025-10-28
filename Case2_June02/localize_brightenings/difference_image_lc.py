
'''
Created on 1 Oct 2025
@author: adithya-hn

- Running code for difference image
'''

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
import suit_difference_image

flt='NB04'
p1,p2,p3,p4=0,550,0,650
data_folder=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case2_June02/data/aligned_crop/*{flt}.fits'
clip_bar=False  #to remove the ccd quadrant
c1,c2=80,120
thresh_sig=4

get_histograms=False
suit_difference_image.diff_img(data_folder,flt,c1,c2,thresh_sig,clip_bar,get_histograms)

