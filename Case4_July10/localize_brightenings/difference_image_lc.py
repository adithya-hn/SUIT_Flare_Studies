
'''
Created on 1 Oct 2025
@author: adithya-hn

- Running code for difference image
'''

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
import suit_difference_image
#import v2_diff_img

p1,p2,p3,p4=40,670,40,680
flt='NB04' #NB04,NB08
data_folder=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/aligned_crop/*{flt}.fits'
clip_bar=False
c1=70
c2=130
thresh_sig=5
get_histograms=False
save_fits=True

#suit_difference_image.diff_img(data_folder,flt,thresh_sig=thresh_sig,save_fits=False,text_col='k')

suit_difference_image.diff_img(data_folder,flt,thresh_sig=thresh_sig,get_histograms=False,text_col='k') #,get_histograms=True,save_fits=True