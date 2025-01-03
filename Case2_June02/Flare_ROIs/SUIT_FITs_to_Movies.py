
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
#from colormap import filterColor


start = timeit.default_timer()
now = datetime.datetime.now()-timedelta(days=1)
fol_nm='/data/sreejith/MCNS_POC/Flare_ROIs/May_19'



Filters=['NB03','NB04']

print('Checking Data availabilty')
search_fold='/scratch/suit_data/level1fits/2024/06/05/normal_roi/' 

for flt in Filters:
    files = sorted(glob.glob(search_fold + '*3'+flt+'.fits'))
    print('Total ',flt ,' files:',len(files))


print('------------------------------')
fltr='NB04'
files = sorted(glob.glob(search_fold + '*3'+fltr+'.fits'))
files=files[100:110]
print('Making movie of ',fltr ,' files:',len(files))

aln_imgs=[]
Sequence = sunpy.map.Map(files, sequence=True)  
#for l in range(len(Sequence)):
#    Sequence[l].meta.update({'CROTA2':0})
mv_nm='June_05_'+fltr+'.mp4'
fig = plt.figure()
ax = fig.add_subplot(projection=Sequence.maps[0])
ani = Sequence.plot(axes=ax,cmap='gray')
#plt.axis('off')
ani.save(mv_nm)
plt.close()
print('Done')
    