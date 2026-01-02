from collections import defaultdict
import sunpy.map as smap
import matplotlib.pyplot as plt 
import glob 
from sunpy.coordinates import Helioprojective, SphericalScreen, propagate_with_solar_surface
from sunpy.map import Map, MapSequence
from astropy.coordinates import SkyCoord
import astropy.units as u
from sunkit_image.coalignment import apply_shifts, calculate_match_template_shift
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
import os
from matplotlib.widgets import RectangleSelector
import numpy as np
from astropy.io import fits
import datetime
import roi_correction
import pathlib

raw_files_pth='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case9_Nov13/data/raw/'
filters=['NB02','NB05','NB06','NB07']

#------------------------------------------------------------------------------------------

pathlib.Path(raw_files_pth[:-4]+'Processed/contam_corr_data/').mkdir(parents=True,exist_ok=True)
pathlib.Path(raw_files_pth[:-4]+f'contam_correct').mkdir(parents=True, exist_ok=True)

for flt in filters:
    print(f' Correcting {flt}')
    raw_files=glob.glob(raw_files_pth+f'*3{flt}.fits')
    raw_files=sorted(raw_files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    MODE='median' # options: 'median' or 'max'
    PLOT= True # Plot preview?
    SAVE=True # Save o/p fits?
    print('No of files: ', len(raw_files))
    groups = defaultdict(list)
    maps=Map(raw_files,sequence=True)

    for i in range(len(maps)):
        mp= maps[i]
        x1 = mp.meta['X1'] #col
        y1 = mp.meta['Y1'] #row
        groups[(x1, y1)].append(raw_files[i])
    
    i=0

    for (x1, y1), f_seq in groups.items():
        print
    
        if len(f_seq) < 1:
            print('Not enough files to stcak')
            continue
        i+=1
        print(i,') ',x1,y1,'-- No. of files in batch: ',len(f_seq))
        if i==1: #taking referenceposition
            row=y1
            col=x1
        
        flat_shift=(-(y1-row),-(x1-col))
        print('Shift applied for flat ',flat_shift)
        savepath= raw_files_pth[:-4]+f'contam_correct/roi_flat_frame_{i}.fits'
        seq = Map(f_seq, sequence=True)
        
        print('  ')
        
        if i==1:
            aligned_sequence= roi_correction.align_maps(seq[:10]) # Generate flat using first 10 images of sequence
            flat_frame_= roi_correction.generate_flat(aligned_sequence,savepath, SAVE)
    
        flat_frame=np.roll(flat_frame_,shift=flat_shift,axis=(0,1))

        ref_map= aligned_sequence[0] 
        corrected_map_ls=[]
        for m in seq: #Multiprocess to be implemented here
            corrected_map= Map(m.data/flat_frame, m.meta)
            img_savepath=  raw_files_pth[:-4]+'Processed/contam_corr_data/'+m.meta['F_NAME']
            corrected_map.save(img_savepath, overwrite=True)
            