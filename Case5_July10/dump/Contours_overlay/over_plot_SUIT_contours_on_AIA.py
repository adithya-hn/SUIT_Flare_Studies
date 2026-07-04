import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
from sunkit_image.coalignment import calculate_match_template_shift,apply_shifts
from datetime import timedelta
import timeit
import pathlib
from astropy.coordinates import SkyCoord
import numpy as np

from sunpy.time import parse_time
from astropy.time import Time
from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord


fltr='NB08'
fltr2='171'
search_fold=f'/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed_2/Aligned_images/{fltr}/' #Custom Folder
base_fold=f'/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/AIA_Data/{fltr2}/'    
print(f'Searching for {fltr} images in {search_fold} folder')

fdir =search_fold 
files = glob.glob(fdir + '*3'+fltr+'.fits')
files=sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
files=files
b_files=glob.glob(base_fold + '*.fits')
print('Total files:',len(files))
base_time_array=[]

#--- Time using file name---
for f in range(len(b_files)):
    print(os.path.basename(b_files[f])[24:46])
    base_time_array.append(datetime.datetime.strptime(os.path.basename(b_files[f])[24:46], "%Y-%m-%dT%H:%M:%S.%f"))
base_time_array=Time(parse_time(base_time_array)) #np.array(base_time_array)
#---

# create folder for images
fol_nm=os.getcwd() 
jpg_fold=fol_nm+'/'+'Contour_imgs'

pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold+f'/{fltr2}').mkdir(parents=True, exist_ok=True)

plot_data=[]
tot_count=[]
dates=[]

for i in range(len(files)):
    suitMap=sunpy.map.Map(files[i])
    base_time=Time(parse_time(suitMap.date))
    idx=np.argmin(np.abs(base_time_array - base_time))# find the closesed image timestamp 

    suit_pos = get_horizons_coord(-21, suitMap.date)
    suitMap.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun)) #header correctin based on SOHO pos

    BaseMap=sunpy.map.Map(b_files[idx])
    base_data=BaseMap.data#/int(BaseMap.meta.get('EXPTIME'))
    Base_img=sunpy.map.Map(base_data,BaseMap.fits_header)

    img_head=suitMap.fits_header
    suit_data=suitMap.data #gaussian_filter(suitMap.data,sigma=sigma)
    alned_data=suit_data*1000/int(suitMap.meta.get('MEAS_EXP'))
    
    #Contour levels
    th_lvs=[3900,4100,4300]
    th_lvs2=[12000,13000]
    Thresh_val=3600
    
    

    fl_nm=jpg_fold+f'/{fltr2}'+'/'+os.path.basename(files[i])[:-4]+'jpg'
    fig=plt.figure(figsize=(10,10))

    ax = fig.add_subplot(111, projection=suitMap) 
    Base_img.plot(cmap='gray',autoalign=True)

    suitMap.draw_contours(axes=ax, levels=th_lvs,zorder=1,colors=['blue','red', 'yellow'])
   
    #suitMap.draw_limb(axes=ax)
    #suitMap.draw_grid(axes=ax)
    plot_str='Ca II h: '+str(base_time_array[idx])
    ax.text(50,50, plot_str, color='white', fontsize=10)
    plt.draw()
    plt.colorbar()
    #plt.title('NB03 diff: '+str(suitMap_.date))
    #plt.axis('off')
    #plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    #plt.tight_layout()
    plt.savefig(fl_nm)
    plt.close()    
    #print(i,' / ',len(files))
plot_data=np.array(plot_data)
tot_count=np.array(tot_count)
np.savetxt(f'{fltr2}_threshold_count.csv',np.c_[dates,plot_data,tot_count],delimiter=',',fmt='%s')


