
'''
Created on 1 Oct 2025
@author: adithya-hn

- Running code for difference image
'''

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
import suit_difference_image

flt='NB04'
p1,p2,p3,p4=130,690,20,560
data_folder=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/1600_aligned/*{flt}.fits'
clip_bar=False  #to remove the ccd quadrant
c1,c2=80,120

get_histograms=False
suit_difference_image.diff_img(data_folder,flt,p1,p2,p3,p4,c1,c2,clip_bar)

