import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import pandas as pd
import datetime
import timeit
import pathlib
from astropy.coordinates import SkyCoord
import numpy as np
import glob
from scipy import stats

start = timeit.default_timer()

# Parameters

#csv_path = 'Flare_files_Nov11_M1.4_case28f.dat'  # <- change this to your actual file path
fdir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop_fits'
fol_nm = os.getcwd() + '/lc_images_/'
Filters = ['NB04']
def get_mg_threshold(ref_mg_rot):
    valid = ref_mg_rot.data[(ref_mg_rot.data > 100)]
    valid_int = valid.astype(int)  
    mode_val = stats.mode(valid_int, keepdims=True).mode[0]/ref_mg_rot.meta.get('CMD_EXPT')*1000  #normalised mode value    
    th_lvs=np.array([mode_val*0.7,mode_val*1.5,mode_val*3] )
    return mode_val,th_lvs

# Loop over filters 
for fltr in Filters:
    # Filter files for the current filter
    #files = sorted([f for f in all_files if f'*3{fltr}.fits' in f or f.endswith(f'3{fltr}.fits')])
    
    files = sorted(glob.glob(fdir+ '/*3'+fltr+'.fits'))
    if not files:
        print(f"No files found for filter {fltr}")
        continue

    print('Processing filter:', fltr)
    plot_data = []
    Plot_data_er = []
    dates = []
    box_pth = f'{fol_nm}/{fltr}/Box'
    pathlib.Path(f'{fol_nm}/{fltr}').mkdir(parents=True, exist_ok=True)
    pathlib.Path(box_pth).mkdir(parents=True, exist_ok=True)

    fltr_count = []
    date_array = []
    fltr_count_err = []
    bx_area = []
    er_bx_area = []

    Sequence = sunpy.map.Map(files, sequence=True)
    aligned_maps = Sequence
    ref_map=Sequence[0]
    thresh_lvs=get_mg_threshold(ref_map)

    for i, suit_map in enumerate(Sequence):
        fnm = os.path.basename(files[i])[:-5]
        F_name = f'{fol_nm}/{fltr}/{fnm}.jpg'
        Box_fnm = f'{box_pth}/{fnm}.jpg'
        Headr_data = Sequence[0].fits_header
        Headr_data['DATE-OBS'] = str(aligned_maps[i].date)

        # Intensity scaling
        Imn, Imx = 1000, 30000
        if fltr == 'NB08': Imx = 12000
        elif fltr == 'NB04': Imx = 33000
        elif fltr == 'BB01': Imx = 21000
        elif fltr in ['BB02', 'BB03']: Imx = 26000

        # Plot full map
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection=suit_map)
        suit_map.plot(axes=ax,cmap='suit_nb04', vmin=Imn, vmax=Imx)
        ref_map.draw_contours(axes=ax,levels=thresh_lvs[1],colors=['red','pink','green'])
        plt.colorbar()
        plt.savefig(F_name, dpi=300)
        plt.close()

stop = timeit.default_timer()
print('Run Time: ', (stop - start)/60, 'Mins')
