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


channel= 'HMI'  # Channel name for the output file

files = sorted(glob(f"/media/adithya/Adi_disk4/SUIT_flare_work/case4_jul10/data/hmi/HMI_cutouts/*.fits")) 

pathlib.Path('HMI_box/er_box').mkdir(parents=True, exist_ok=True) 

times = []
intensities = []
print("Number of files found:", len(files)) 
test_map = sunpy.map.Map(files[0])
test_map.plot_settings["norm"]
vmin = -100
vmax = 100

# Set a custom normalization
test_map.plot_settings['norm'] = simple_norm(test_map.data,  'linear', vmin=vmin, vmax=vmax)
test_map.peek()

#sys.exit(1)

plot_data=[]
tot_count=[]
qs=[]
dates=[]
test_point=[]
i=0
for file in tqdm(files):
    suit_map=sunpy.map.Map(file)
    i+=1
    

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(projection=suit_map)
    suit_map.plot(axes=ax)

    rotation_angle=7.15 #suit_map.meta["CROTA2"]

    #AR sub box
    cen_cord      = SkyCoord(-255 * u.arcsec, -310 * u.arcsec, frame=suit_map.coordinate_frame)
    diffrot_point = SkyCoord(RotatedSunFrame(base=cen_cord, duration=45*(i)*u.second))
    transformed_diffrot_point = diffrot_point.transform_to(suit_map.coordinate_frame)

    offset_frame1 = SkyOffsetFrame(origin=transformed_diffrot_point, rotation=-rotation_angle*u.deg)
    width1  = 365 * u.arcsec
    height1 = 240 * u.arcsec
    coords = SkyCoord(lon=[-1/2, 1/2] * width1, lat=[-1/2, 1/2] * height1, frame=offset_frame1)

    #Quiet sub box
    #center_coord4 = SkyCoord(-170 * u.arcsec, -500 * u.arcsec, frame=suit_map.coordinate_frame)
    center_coord4 = SkyCoord(-275 * u.arcsec, -540 * u.arcsec, frame=suit_map.coordinate_frame) # point left, belwoe th ar box
    #center_coord4 = SkyCoord(-350 * u.arcsec, -570 * u.arcsec, frame=suit_map.coordinate_frame)

    #center_coord4 = SkyCoord(-170 * u.arcsec, -500 * u.arcsec, frame=suit_map.coordinate_frame)
    diffrot_point = SkyCoord(RotatedSunFrame(base=center_coord4, duration=45*(i)*u.second))
    transformed_diffrot_point_1 = diffrot_point.transform_to(suit_map.coordinate_frame)

    offset_frame4 = SkyOffsetFrame(origin=transformed_diffrot_point_1, rotation=-rotation_angle*u.deg)
    width4  = 40 * u.arcsec
    height4 = 40 * u.arcsec
    er_coords=SkyCoord(lon=[-1/2, 1/2] * width4, lat=[-1/2, 1/2] * height4, frame=offset_frame4)
    suit_map.draw_quadrangle(er_coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='Background')
    suit_map.draw_quadrangle(coords,axes=ax,edgecolor="red",linestyle="-",linewidth=2,label='Region of interest')
    
    er_box=suit_map.submap(er_coords)
    test_box=suit_map.submap(coords)

    Thresh_val= np.mean(er_box.data)*3
    fl_nm=f'HMI_box/'+os.path.basename(file)[:-4]+'jpg'
    #print(fl_nm)
    plt.savefig(fl_nm,dpi=300)
    plt.close()

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(projection=er_box)
    er_box.plot_settings['norm'] = simple_norm(er_box.data,  'linear', vmin=vmin, vmax=vmax)
    er_box.plot(axes=ax)
    plt.colorbar(ax=ax)
    plt.savefig(f'HMI_box/er_box/'+os.path.basename(file)[:-4]+'jpg',dpi=300)
    plt.close()


    qs.append(np.mean(abs(er_box.data)))
    dates.append(suit_map.date.datetime)
    test_point.append(np.mean(abs(test_box.data))) #test_box.meta.get("CMD_EXPT") test_box.meta.get("SHTR_STR"))#

dates=np.array(dates)
qs=np.array(qs)
test_point=np.array(test_point)
np.savetxt(f'HMI_qs1_count.csv',np.c_[dates,qs,test_point],comments='',header='Time,QS_mean,AR_mean',delimiter=',',fmt='%s')
plt.figure(figsize=(10, 5))
ax=plt.subplot(111)
ax1=ax.twinx()

date_stamp=dates[0].strftime('%Y-%m-%d')
plt.title(f'HMI AR and QS Box Light curve ({date_stamp})')
ax.plot(dates, test_point, marker='o',markersize=0.5,label='AR box Intensity')
ax1.plot(dates, qs,'r', marker='o',markersize=0.5,label='QS2 box Intensity')
plt.gcf().autofmt_xdate()
plt.xlabel("Time")
ax1.set_ylabel("QS Intensity")
ax.set_ylabel("AR Intensity")
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig(f'HMI_box_QS2_intensity_.png', dpi=300, bbox_inches='tight')
plt.close()

'''np.savetxt(f'{channel}_roi_lc.csv',np.c_[times,intensities],delimiter=',',fmt='%s',header='Date,Intensity')
plt.figure(figsize=(10, 5))
plt.plot(times, intensities, marker='o',markersize=0.5)
plt.gcf().autofmt_xdate()
plt.xlabel("Time")
plt.ylabel("Total Intensity")
plt.title(f"AIA {channel} Å Light Curve (roi)")
plt.grid(True)
plt.yscale('log')
plt.savefig(f'AIA_{channel}_roi_lc.png', dpi=300, bbox_inches='tight')
plt.show()'''