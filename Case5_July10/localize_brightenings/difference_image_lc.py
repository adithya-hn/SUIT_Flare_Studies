
'''
Created on 1 Oct 2025
@author: adithya-hn

- Running code for difference image
'''

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
import suit_difference_image

flt='NB04'
p1,p2,p3,p4=50,670,10,520
data_folder=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/1600_aligned/*{flt}.fits'

clip_bar=False
c1=70
c2=130
get_histograms=False
thresh_sig=4
suit_difference_image.diff_img(data_folder,flt,c1,c2,thresh_sig,clip_bar)

