from glob import glob
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import sunpy.map
import numpy as np
import astropy.units as u
import sys
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from sunpy.coordinates import RotatedSunFrame
import pathlib
import os
import datetime
from sunpy.time import parse_time
from astropy.time import Time
from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord
from skimage.morphology import disk, closing,opening,dilation
from sunpy.coordinates import Helioprojective, SphericalScreen, propagate_with_solar_surface


channel= 'HMI'  # Channel name for the output file

hmi_files = sorted(glob(f"/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case2_June02/data/hmi/HMI/HMI_cutouts/*.fits")) 
suit_files=sorted(glob(f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case2_June02/data/1600_aligned/*NB03.fits'))

#--------------------------------------------------------------------------
base_time_array=[]
mg_map_time=[]
for f in range(len(hmi_files)):
    base_time_array.append(datetime.datetime.strptime(os.path.basename(hmi_files[f])[10:25], "%Y%m%d_%H%M%S"))

for j in range(len(suit_files)):
    mg_map_time.append(datetime.datetime.strptime(os.path.basename(suit_files[j]).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    mg_map_time_array=Time(parse_time(mg_map_time))

jpg_fold=os.getcwd()
pathlib.Path('HMI_box/er_box').mkdir(parents=True, exist_ok=True) 

times = []
intensities = []
print("Number of files found:", len(hmi_files)) 


#sys.exit(1)

plot_data=[]
tot_count=[]
qs=[]
dates=[]
test_point=[]
i=0
ref_mag_map=sunpy.map.Map(hmi_files[0])
for file in tqdm(hmi_files):
    hmi_map=sunpy.map.Map(file)
    hmi_time=Time(parse_time(hmi_map.date))
    idx=np.argmin(np.abs(mg_map_time_array - hmi_time))
    suitMap=sunpy.map.Map(suit_files[idx]) 
    suit_pos = get_horizons_coord(-21, suitMap.date)
    suitMap.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
    MgII_data=suitMap.data*1000/int(suitMap.meta.get('CMD_EXPT'))

    flt_th_lvs=[-100,100]
        
    '''
    Thresh1_data=np.sum(np.where(abs(base_data)>flt_th_lvs[0],abs(base_data),0))
    Thresh2_data=np.sum(np.where(abs(base_data)>flt_th_lvs[1],abs(base_data),0))
    c1_data.append(Thresh1_data)
    c2_data.append(Thresh2_data)
    dates.append(base_time)
    '''

   
    th_lvs2=[-100,100]*u.G
    mask = np.abs(hmi_map.data)>50 
    kernel=disk(3)
    mask_=dilation(mask,kernel)
    mask_=dilation(mask_,kernel)
    mask_=dilation(mask_,kernel)
    masked_hmi = np.zeros_like(hmi_map.data)
    masked_hmi[mask_] = 1  # binary mask

    mask_map=sunpy.map.Map(masked_hmi,hmi_map.meta)
    with propagate_with_solar_surface():
            mask_map=mask_map.reproject_to(ref_mag_map.wcs) 


    fl_nm=jpg_fold+f'/'+os.path.basename(suit_files[idx])[:-4]+'jpg'
    fig=plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection=suitMap)
    
    suitMap.plot(cmap='gray',autoalign=True)
    mask_map.plot(axes=ax,cmap='viridis',alpha=0.5)
    
    #hmi_map.draw_contours(axes=ax, levels=th_lvs2,zorder=2,colors=['pink','green'],alpha=0.4)
    #ax.imshow(masked_hmi, origin="lower", alpha=0.3, cmap="Reds")

    plt.colorbar()
        
    plot_str='Ca II h: '+str(suitMap.date)
    ax.text(50,50, plot_str, color='white', fontsize=10)
    plt.draw()
    #plt.colorbar()
    plt.savefig(fl_nm)
    plt.close()
        

