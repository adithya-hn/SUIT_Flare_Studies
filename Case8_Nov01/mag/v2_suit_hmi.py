from glob import glob
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')
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
from skimage.measure import label, regionprops
from skimage import filters, measure
from scipy import stats
from scipy.optimize import curve_fit
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

import warnings
import logging
warnings.simplefilter('ignore')
log_ = logging.getLogger('sunpy')
log_.setLevel('WARNING')
logging.getLogger('reproject').setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.ERROR)


channel= 'HMI'  # Channel name for the output file

hmi_files = sorted(glob(f"/media/adithya/Adi_disk4/SUIT_flare_work/case8_nov01/data/HMI/HMI_cutouts/*.fits")) 
suit_files=sorted(glob(f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/data/aligned_crop/*NB04.fits'))

#--------------------------------------------------------------------------

def get_suit_scale_rebined_map(ref_fd_1600,ref_mg_rot):
    scale=ref_fd_1600.scale[0].value/ref_mg_rot.scale[0].value
    fd_new_dem=[ref_fd_1600.data.shape[1]*scale,ref_fd_1600.data.shape[0]*scale]*u.pixel
    ref_aia_resmp=ref_fd_1600.resample(fd_new_dem)
    blo = SkyCoord(ref_mg_rot.bottom_left_coord.Tx, ref_mg_rot.bottom_left_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    tro = SkyCoord(ref_mg_rot.top_right_coord.Tx, ref_mg_rot.top_right_coord.Ty, frame=ref_aia_resmp.coordinate_frame)
    return ref_aia_resmp.submap(blo,top_right=tro)  


hmi_time=[]
mg_map_time=[]
for f in range(len(hmi_files)):
    hmi_time.append(datetime.datetime.strptime(os.path.basename(hmi_files[f])[10:25], "%Y%m%d_%H%M%S"))
    hmi_time_array=Time(parse_time(hmi_time))

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




suit_ref_map=sunpy.map.Map(suit_files[0])
print('ref file: ',suit_files[0])
img_mode=[]
qs_mode=[]
qs_median=[]
qs_mean=[]
qs_std=[]
Date_time=[]
ar_count=[]
ar_mode=[]
ar_1_count=[]
ar_1_mode=[]
ar_sigma=[]
ar1_sigma=[]
for file in tqdm(suit_files):
    suitMap=sunpy.map.Map(file)
    suit_time=Time(parse_time(suitMap.date))
    idx=np.argmin(np.abs(hmi_time_array - suit_time))
    ref_mg_rot=suitMap.rotate(angle=suitMap.meta["P_ANGLE"] * u.deg,missing=0)
    hmi_map_raw=sunpy.map.Map(hmi_files[idx])
    hmi_map=get_suit_scale_rebined_map(hmi_map_raw,ref_mg_rot)
    blo = SkyCoord(suitMap.bottom_left_coord.Tx, suitMap.bottom_left_coord.Ty, frame=suitMap.coordinate_frame)
    tro = SkyCoord(suitMap.top_right_coord.Tx, suitMap.top_right_coord.Ty, frame=suitMap.coordinate_frame)
    suit_pos = get_horizons_coord(-21, suitMap.date)
    suitMap.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
    #ref_aia_resmp.submap(blo,top_right=tro)  
    #suitMap=suitMap.submap(coords)

    MgII_data=suitMap.data*1000/int(suitMap.meta.get('CMD_EXPT'))
    print('img mode',stats.mode(MgII_data.astype(int), axis=None, keepdims=True)[0][0][0])
    img_mode.append(stats.mode(MgII_data.astype(int), axis=None, keepdims=True)[0][0][0])
    th_lvs2=[-100,100]*u.G
    mask = np.abs(hmi_map.data)>50
    ar_mask=np.abs(hmi_map.data)>200 
    kernel=disk(8) #np.array([[1,1,1,1,1,1],[1,1,1,1,1,1],[1,1,1,1,1,1]])#
    #print(kernel)
    mask_=dilation(mask,kernel)
    #mask_=dilation(mask_,kernel)
    # mask_=dilation(mask_,kernel)
    ar_mask_=dilation(ar_mask,kernel)
    #ar_mask_=dilation(ar_mask_,kernel)
    #ar_mask_=dilation(ar_mask_,kernel)
    masked_hmi = np.zeros_like(hmi_map.data)
    ar_masked_hmi = np.zeros_like(hmi_map.data)
    masked_hmi[mask_] = 1  # binary mask
    ar_masked_hmi[ar_mask_]=1

    mask_map=sunpy.map.Map(masked_hmi,hmi_map.meta)
    ar_mask_map=sunpy.map.Map(ar_masked_hmi,hmi_map.meta)
    with propagate_with_solar_surface():
            mask_map=mask_map.reproject_to(suit_ref_map.wcs) 
            ar_mask_map=ar_mask_map.reproject_to(suit_ref_map.wcs) 

    fl_nm=jpg_fold+f'/qs_mask/'+os.path.basename(file)[:-4]+'jpg'
    fig=plt.figure()
    ax = fig.add_subplot(111, projection=suitMap)
    #suitMap.plot(cmap='gray',autoalign=True)
    #mask_map.plot(axes=ax,cmap='viridis',alpha=0.5)
    binary_img=(mask_map.data).astype(np.bool)
    binary_ar=(ar_mask_map.data).astype(np.bool)
    inv=np.invert(binary_img)
    inv_ar=np.invert(binary_ar)

    plt.imshow(MgII_data*inv,vmin=0,vmax=12000)
    contours = measure.find_contours(mask_map.data, level=0.5)
    for contour in contours:
        ax.plot(contour[:, 1], contour[:, 0], linewidth=1.5, color='red',alpha=0.5)
    #plt.imshow(mask_map.data,cmap='gray',alpha=1)      
    #plt.imshow(ar_mask_map.data,cmap='Oranges',alpha=0.3)
    
    plt.colorbar()
    #plt.grid()
    plt.savefig(fl_nm)
    plt.close()
    bins=np.arange(0,30000,5)
    hist_data=np.array((MgII_data*inv).flatten(),dtype=np.int64)
    hist_data=hist_data[hist_data>0]

    
    qs_mode.append(stats.mode((hist_data))[0])
    qs_median.append(np.median(hist_data))
    qs_mean.append(np.mean(hist_data))
    qs_std.append(np.std(hist_data))
    Date_time.append(suitMap.date.datetime)
    print('--->',stats.mode((hist_data))[0],np.median(hist_data),np.mean(hist_data),np.std(hist_data))
    
    masked_ar1=np.array((MgII_data*binary_img).flatten(),dtype=np.int64)
    masked_ar=np.array((MgII_data*binary_ar).flatten(),dtype=np.int64)
    masked_ar1=masked_ar1[masked_ar1>0] #avoid masked zeros
    masked_ar =masked_ar[masked_ar>0]
    ar_1_mode.append(stats.mode(masked_ar1)[0])
    ar_mode.append(stats.mode(masked_ar)[0])
    ar_count.append(np.mean(masked_ar))
    ar_1_count.append(np.mean(masked_ar1))
    ar_sigma.append(np.std(masked_ar))
    ar1_sigma.append(np.std(masked_ar1))
    
    #print(qs_count)
        
   
    bins=np.arange(0,20000,10)
        # n, bins, patches = plt.hist(hist_data, bins=Bins, color='gray', alpha=0.7, edgecolor='black')
        # bin_centers = (bins[:-1] + bins[1:]) / 2
        # max_bin_index = np.argmax(n)
        # qs_guess = bin_centers[max_bin_index]

        # def double_gaussian(x, amp1, mean1, sigma1, amp2, mean2, sigma2):
        #     gauss1 = amp1 * np.exp(-(x - mean1)**2 / (2 * sigma1**2))
        #     gauss2 = amp2 * np.exp(-(x - mean2)**2 / (2 * sigma2**2))
        #     return gauss1+gauss2
        # p0 = [np.max(n), qs_guess, 100,  # First Gaussian (QS)
        #     np.max(n)/5, qs_guess+500, 200]  # Second Gaussian (active region)
        
        # params, cov = curve_fit(double_gaussian, bin_centers, n, p0=p0,maxfev=5000)
        # amp1, mean1, sigma1, amp2, mean2, sigma2 = params


        # print(f"QS (Gaussian 1): Mean = {mean1:.2f}, Sigma = {sigma1:.2f}")
        # print(f"Bright Features (Gaussian 2): Mean = {mean2:.2f}, Sigma = {sigma2:.2f}")
    plt.hist(hist_data,bins=bins)
    plt.axvline(stats.mode((hist_data))[0],label='mode',color='k',lw=0.5)
    plt.axvline(np.median(hist_data),label='median',color='y',lw=0.5)
    plt.axvline(np.mean(hist_data),label='mean',color='r',lw=0.5)
    plt.axvline((np.mean(hist_data)+3*np.std(hist_data)),label='mean+3std')
    plt.savefig(f'histograms/{hmi_map.date}.png')
    plt.close()
    
    fig=plt.figure()
    ax = fig.add_subplot(111, projection=suitMap)
    plt.imshow(MgII_data*binary_ar)
    plt.colorbar()
    plt.savefig(f'ar_mask/Ar_mask_{suitMap.date}.png',dpi=300)
    plt.close()

    
         
np.savetxt('QS_data.csv',np.c_[Date_time,img_mode,qs_mode,qs_median,qs_mean,qs_std,ar_1_mode,ar_1_count,ar1_sigma,ar_mode,ar_count,ar_sigma],fmt='%s',delimiter=',',comments='',header='Date,img_mode,QS_mode,QS_median,QS_mean,QS_std,AR_mode,AR mean,Ar_sigma,ARmode_200G,ARmean200G,AR_sigma200G')

# qs_count=np.array(qs_count)
# ar_count=np.array(ar_count)
#Date_time=np.array(Date_time, dtype='datetime64')
'''fig, axs = plt.subplots(1, 1, figsize=(5,8))
ax1=ax.twinx()
ax.plot(Date_time,qs_count,label='QS_value')
ax1.plot(Date_time,ar_count,label='AR_count')
ax.set_xlabel('Time')


ax1.set_ylabel('Quiet sun value')
ax.set_ylabel('AR count')
ax.legend()
ax1.legend()'''
plt.figure(figsize=(12,6))
plt.plot(Date_time,qs_mode,label='qs mode')
plt.plot(Date_time,qs_mean,label='qs mean')
plt.plot(Date_time,qs_median,label='qs median')
plt.legend()
plt.savefig('QS Count.png',dpi=300)
plt.show()
        

