


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
from colormap import filterColor


start = timeit.default_timer()
now = datetime.datetime.now()-timedelta(days=1)
fol_nm='/data/sreejith/MCNS_POC/Flare_ROIs/'

jpg_fold=fol_nm+'/'+'Coloured_ROI_imgs_01'
algn_dir=fol_nm+'/'+'Aligned_images_01'

#Filter='NB03'
Filters=['NB01','NB02','NB03','NB04','NB05','NB06','NB07','NB08']#,'BB01','BB02','BB03']
pathlib.Path(algn_dir).mkdir(parents=True, exist_ok=True)
pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
search_fold='/scratch/suit_data/level1.1fits/2024/06/01/normal_roi/' #'/scratch/suit_data/level1fits/2024/'+str(now.month).zfill(2)+'/'+str(now.day).zfill(2)+'/'+'normal_2k'+'/'
fdir =search_fold 
for fltr in Filters:
    pathlib.Path(jpg_fold+'/'+fltr).mkdir(parents=True, exist_ok=True)
    pathlib.Path(algn_dir+'/'+fltr).mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(fdir + '*3'+fltr+'.fits'))
    print('Total ',fltr ,' files:',len(files))
   
    aln_imgs=[]
    Sequence = sunpy.map.Map(files, sequence=True)
    mv_nm=fltr+'.mp4'
    #fig = plt.figure()
    #ax = fig.add_subplot(projection=Sequence.maps[80])
    #ani = Sequence.plot(axes=ax,cmap=filterColor[fltr])
    #plt.axis('off')
    #ani.save(mv_nm)
    #plt.close()
    
    aligned_maps=mc_coalign(Sequence,layer_index=70,clip=False)
    print('Aligned, Saving now...')
    for j in range(len(aligned_maps)):  #account for incerted image, j=0 is refference image
        aligned_img=sunpy.map.Map(aligned_maps[j].data,aligned_maps[j].fits_header)
        aligned_img.save(algn_dir+'/'+fltr+'/'+os.path.basename(files[j]),overwrite=True) #need this for alignement refference
        aln_imgs.append(algn_dir+'/'+fltr+'/'+os.path.basename(files[j]))
        fl_nm=jpg_fold+'/'+fltr+'/'+os.path.basename(files[j])[:-4]+'jpg'
        #Title_=aligned_img.meta.get('FTR_NAME') +' Filter: '+aligned_img.date
        fig=plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection=aligned_img)
        aligned_img.plot(title=False,cmap='gray',clip_interval=(2, 99.95) * u.percent)
        plt.figtext(0.04, 0.04,aligned_img.date, color='white',weight='bold', fontsize=18)
        plt.draw()
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(fl_nm,bbox_inches='tight', transparent="True", pad_inches=0)
        #plt.show()
        plt.close()
 

stop = timeit.default_timer()

print('Run Time: ', (stop - start)/60,'Mins') 