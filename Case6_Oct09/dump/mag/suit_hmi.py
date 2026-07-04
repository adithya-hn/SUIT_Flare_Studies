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
from reproject import reproject_exact

channel= 'HMI'  # Channel name for the output file
fltr='NB08'
hmi_files = sorted(glob(f"/media/adithya/Adi_disk4/SUIT_flare_work/case6_oct09/data/HMI/HMI_cutouts/*.fits")) 
suit_files=sorted(glob(f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/aligned_crop/*{fltr}.fits'))

#--------------------------------------------------------------------------
hmi_time=[]
mg_map_time=[]
for f in range(len(hmi_files)):
    hmi_time.append(datetime.datetime.strptime(os.path.basename(hmi_files[f])[10:25], "%Y%m%d_%H%M%S"))
    hmi_time_array=Time(parse_time(hmi_time))

for j in range(len(suit_files)):
    mg_map_time.append(datetime.datetime.strptime(os.path.basename(suit_files[j]).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    mg_map_time_array=Time(parse_time(mg_map_time))

jpg_fold=os.getcwd()
pathlib.Path('qs_mask').mkdir(parents=True, exist_ok=True) 
pathlib.Path('ar_mask').mkdir(parents=True, exist_ok=True) 

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

cTx1=100
cTy1=100
arW=320
arH=360

ref_mag_map=sunpy.map.Map(hmi_files[0])
hmi_time=Time(parse_time(ref_mag_map.date))
idx1=np.argmin(np.abs(mg_map_time_array - hmi_time))
suit_ref_map=sunpy.map.Map(suit_files[idx1])
#print(suit_ref_map.dsun)
rotation_angle=suit_ref_map.meta["P_ANGLE"]
fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(projection=suit_ref_map)
cen_cord      = SkyCoord(cTx1 * u.arcsec, cTy1 * u.arcsec, frame=suit_ref_map.coordinate_frame)
offset_frame1 = SkyOffsetFrame(origin=cen_cord, rotation=-rotation_angle*u.deg)
width1  = arW * u.arcsec
height1 = arH * u.arcsec
coords = SkyCoord(lon=[-1/2, 1/2] * width1, lat=[-1/2, 1/2] * height1, frame=offset_frame1)
suit_ref_map.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
suit_ref_map.plot()
plt.close()
suit_ref_map = suit_ref_map.submap(coords)

qs_count=[]
Date_time=[]
ar_count=[]
ar_mode=[]
ar_1_count=[]
ar_1_mode=[]
ar_sigma=[]
ar1_sigma=[]
suit_pos = get_horizons_coord(-21, ref_mag_map.date)
for file in tqdm(suit_files):
    suitMap=sunpy.map.Map(file)
    suit_time=Time(parse_time(suitMap.date))
    idx=np.argmin(np.abs(hmi_time_array - suit_time))
    hmi_map=sunpy.map.Map(hmi_files[idx])
    
    cen_cord      = SkyCoord(cTx1 * u.arcsec, cTy1 * u.arcsec, frame=suitMap.coordinate_frame)
    offset_frame1 = SkyOffsetFrame(origin=cen_cord, rotation=-rotation_angle*u.deg)
    width1  = arW * u.arcsec
    height1 = arH * u.arcsec
    coords = SkyCoord(lon=[-1/2, 1/2] * width1, lat=[-1/2, 1/2] * height1, frame=offset_frame1)
    
    suitMap.meta.update(get_observer_meta(suit_pos, rsun=suit_pos.rsun))
    suitMap=suitMap.submap(coords)
    MgII_data=suitMap.data*1000/int(suitMap.meta.get('CMD_EXPT'))
    
    th_lvs2=[-100,100]*u.G
    mask = np.abs(hmi_map.data)>50 
    ar_mask=np.abs(hmi_map.data)>200 
    kernel=np.array([[1,1,1],[1,1,1],[1,1,1]])#disk(3) #
    #print(kernel)
    mask_=dilation(mask,kernel)
    mask_=dilation(mask_,kernel)
    mask_=dilation(mask_,kernel)
    '''mask_=dilation(mask_,kernel)
    mask_=dilation(mask_,kernel)
    mask_=dilation(mask_,kernel)
    mask_=dilation(mask_,kernel)'''
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
            out_data_qs,fp1=reproject_exact((mask_map.data,mask_map.wcs),suit_ref_map.wcs)
            out_data_ar,fp2=reproject_exact((ar_mask_map.data,ar_mask_map.wcs),suit_ref_map.wcs)
            #mask_map=mask_map.reproject_to(suit_ref_map.wcs) 
            #ar_mask_map=ar_mask_map.reproject_to(suit_ref_map.wcs) 

    mask_map=sunpy.map.Map(out_data_qs,suit_ref_map.meta)
    ar_mask_map=sunpy.map.Map(out_data_ar,suit_ref_map.meta)
    fl_nm=jpg_fold+f'/qs_mask/'+os.path.basename(file)[:-4]+'jpg'
    fig=plt.figure()
    ax = fig.add_subplot(111, projection=suitMap)
    #suitMap.plot(cmap='gray',autoalign=True)
    #mask_map.plot(axes=ax,cmap='viridis',alpha=0.5)
    binary_img=(mask_map.data).astype(np.bool)
    binary_ar=(ar_mask_map.data).astype(np.bool)
    inv=np.invert(binary_img)
    inv_ar=np.invert(binary_ar)

    plt.imshow(MgII_data*inv)
    contours = measure.find_contours(mask_map.data, level=0.5)
    for contour in contours:
        ax.plot(contour[:, 1], contour[:, 0], linewidth=1.5, color='red',alpha=0.5)
    #plt.imshow(mask_map.data,cmap='gray',alpha=1)      
    #plt.imshow(ar_mask_map.data,cmap='Oranges',alpha=0.3)
    
    plt.colorbar()
    plt.savefig(fl_nm)
    plt.close()
    bins=np.arange(0,40000,5)
    hist_data=np.array((MgII_data*inv).flatten(),dtype=np.int64)
    hist_data=hist_data[hist_data>0]

    
    qs_count.append(stats.mode(hist_data)[0])
    Date_time.append(suitMap.date.datetime)
    
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
    
    '''
    plt.hist(masked_ar,bins=bins)
    plt.savefig(f'histograms/{hmi_map.date}.png')
    plt.close()'''
    
    fig=plt.figure()
    ax = fig.add_subplot(111, projection=suitMap)
    plt.imshow(MgII_data*binary_ar)
    plt.colorbar()
    plt.savefig(f'ar_mask/Ar_mask_{suitMap.date}.png',dpi=300)
    plt.close()
np.savetxt(f'{fltr}QS_data.csv',np.c_[Date_time,qs_count,ar_1_mode,ar_1_count,ar1_sigma,ar_mode,ar_count,ar_sigma],fmt='%s',delimiter=',',comments='',header='Date,QS_mode,AR_mode,AR mean,Ar_sigma,ARmode_200G,ARmean200G,AR_sigma200G')

qs_count=np.array(qs_count)
ar_count=np.array(ar_count)
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
plt.plot(Date_time,qs_count)
plt.savefig('Count.png',dpi=300)
plt.show()
        

