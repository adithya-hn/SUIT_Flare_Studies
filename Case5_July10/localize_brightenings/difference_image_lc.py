
'''
Created on 1 Oct 2025
@author: adithya-hn

- Running code for difference image
'''

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
import suit_difference_image
import suit_running_difference_image

flt='NB04' #NB04,NB08
data_folder=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop_fits/*{flt}.fits'
clip_bar=False  #to remove the ccd quadrant
c1,c2=80,120
thresh_sig=5
get_histograms=False
#suit_running_difference_image.running_diff_img(data_folder,flt,thresh_sig=thresh_sig,get_histograms=True,save_fits=False,text_col='k')

suit_difference_image.diff_img(data_folder,flt,thresh_sig=thresh_sig,get_histograms=False,text_col='k') #,get_histograms=True,save_fits=True
