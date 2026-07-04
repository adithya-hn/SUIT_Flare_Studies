
#==========================================
#Differential rotation implimented
#==========================================

import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
from datetime import timedelta
import timeit
import pathlib
from astropy.coordinates import SkyCoord
#from colormap import filterColor
import numpy as np
from sunpy.coordinates import RotatedSunFrame
from tqdm import tqdm


start = timeit.default_timer()
now = datetime.datetime.now()-timedelta(days=1)

search_fold='/media/adithya/Adi_disk4/SUIT_flare_work/case4_jul10/data/HMI/HMI_cutouts/'

fdir =os.listdir(search_fold)
fdir.sort()

param='magnetogram'
plot_data=[]
dates=[]

fol_nm=os.getcwd()
box_pth=fol_nm+'/'+param
pathlib.Path(box_pth+'/Box').mkdir(parents=True, exist_ok=True)
files = sorted(glob.glob(f'{search_fold}*.fits'))
print('Total files:',len(files))
Sequence = sunpy.map.Map(files, sequence=True)
fltr_count=[]
date_array=[]
for i in tqdm(range (len(Sequence))):
    
    hmi_map=Sequence[i]
    #hmi_map.peek()
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(projection=hmi_map)
    #hmi_map.plot(axes=ax)

    point1 = SkyCoord(-173*u.arcsec, -181*u.arcsec, frame=hmi_map.coordinate_frame) #flull region
    diffrot_point1 = SkyCoord(RotatedSunFrame(base=point1, duration=45*(i+1)*u.second))
    transformed_diffrot_point1 = diffrot_point1.transform_to(hmi_map.coordinate_frame)

    point2 = SkyCoord(-213*u.arcsec, -229*u.arcsec, frame=hmi_map.coordinate_frame)
    diffrot_point2 = SkyCoord(RotatedSunFrame(base=point2, duration=45*(i+1)*u.second))
    transformed_diffrot_point2 = diffrot_point2.transform_to(hmi_map.coordinate_frame)
    
    coords = SkyCoord(Tx=(transformed_diffrot_point2.Tx.value, transformed_diffrot_point1.Tx.value) * u.arcsec,Ty=(transformed_diffrot_point2.Ty.value, transformed_diffrot_point1.Ty.value) * u.arcsec,frame=hmi_map.coordinate_frame,)
    hmi_map.draw_quadrangle(coords,axes=ax,edgecolor="blue",linestyle="-",linewidth=2,label='2-element SkyCoord')
    fnm=(os.path.basename(files[i]))[:-5]
    suit_box=hmi_map.submap(coords)
    suit_box.plot(axes=ax)

    #print('-------',f'{fol_nm}/{fnm}.jpg')
    plt.savefig(f'{box_pth}/Box/{fnm}.jpg',dpi=300)
    plt.close()
    
    fltr_count.append(np.sum(abs(suit_box.data)))
    date_array.append(suit_box.date.datetime)
np.savetxt('hmi_heating_box.csv',np.c_[date_array,fltr_count],header='time,HMI_count',comments='',fmt='%s')

plt.plot(date_array,fltr_count)
plt.savefig('HMI_box_lc.png',dpi=300)
plt.show()