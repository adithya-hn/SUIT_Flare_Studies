
'''
Created on 1 Oct 2025
@author: adithya-hn

- Running code for difference image
'''

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
import suit_difference_image

p1,p2,p3,p4=40,670,40,680
flt='NB08' #NB04,NB08
data_folder=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/1600_aligned/*{flt}.fits'
clip_bar=True
c1=70
c2=130
get_histograms=False
suit_difference_image.diff_img(data_folder,flt,p1,p2,p3,p4,c1,c2,clip_bar)

