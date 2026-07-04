import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
from sunkit_image.coalignment import calculate_match_template_shift, apply_shifts
from datetime import timedelta
import timeit
import pathlib
from colormap import filterColor
from astropy.coordinates import SkyCoord
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
from PIL import Image
import pandas as pd
from scipy.ndimage import shift
from scipy.ndimage import gaussian_filter

# Threshold values:
nb3T = 12000
nb4T = 12500
nb8T = 4000
nb3Mx = 14000
nb4Mx = 15000
nb8Mx = 4300

Filters = ['NB03', 'NB04', 'NB08']

for fltr in Filters:
    plot_data = []
    tot_count = []
    dates = []
    areas_level1 = []  # To store the area under the first contour level
    areas_level2 = []  # To store the area under the second contour level
    areas_level3 = []  # To store the area under the third contour level
    actual_counts_level1 = []  # To store the actual count under the first contour level
    actual_counts_level2 = []  # To store the actual count under the second contour level
    actual_counts_level3 = []  # To store the actual count under the third contour level

    search_fold = f'/Analysis/Projects_Data/Flare_Data/July10_Flare_Data2/Processed/Aligned_images/{fltr}/'  # Custom Folder

    print(f'Searching for {fltr} images in {search_fold} folder')
    fdir = search_fold
    files = glob.glob(fdir + '*3' + fltr + '.fits')
    files = sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    print('Total files:', len(files))

    ref_img = sunpy.map.Map(files[0])
    fl_date = datetime.datetime.fromisoformat(str(ref_img.date))
    fol_nm = os.getcwd()
    print(fol_nm)

    jpg_fold = fol_nm + '/' + 'Threshold_imgs'
    algn_imgs = fol_nm + '/' + 'Contours'

    pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)
    pathlib.Path(jpg_fold + f'/{fltr}').mkdir(parents=True, exist_ok=True)
    pathlib.Path(algn_imgs + f'/{fltr}').mkdir(parents=True, exist_ok=True)

    sigma = 1
    if fltr == 'NB03':
        Thresh_val = nb3T
        Tmax = nb3Mx
        th_lvs = [11000, 11500, 13000]

    if fltr == 'NB04':
        Thresh_val = nb4T
        Tmax = nb4Mx
        th_lvs = [11000, 11500, 13000]

    if fltr == 'NB08':
        Thresh_val = nb8T
        Tmax = nb8Mx
        th_lvs = [3900, 4100, 4300]

    for i in range(len(files)):
        suitMap = sunpy.map.Map(files[i])
        img_head = suitMap.fits_header
        suit_data = suitMap.data  # gaussian_filter(suitMap.data,sigma=sigma)
        alned_data = suit_data * 1000 / suitMap.meta.get('MEAS_EXP')
        alignedMap = sunpy.map.Map(alned_data, img_head)
        dates.append(suitMap.date)

        # Calculate area and actual count for each contour level
        area_level1 = np.count_nonzero(alned_data > th_lvs[0])
        area_level2 = np.count_nonzero(alned_data > th_lvs[1]) 
        area_level3 = np.count_nonzero(alned_data > th_lvs[2])

        actual_count_level1 = np.sum(alned_data[alned_data > th_lvs[0]])
        actual_count_level2 = np.sum(alned_data[alned_data > th_lvs[1]])
        actual_count_level3 = np.sum(alned_data[alned_data > th_lvs[2]])

        areas_level1.append(area_level1)
        areas_level2.append(area_level2)
        areas_level3.append(area_level3)
        actual_counts_level1.append(actual_count_level1)
        actual_counts_level2.append(actual_count_level2)
        actual_counts_level3.append(actual_count_level3)

        #fl_nm = jpg_fold + f'/{fltr}/' + os.path.basename(files[i])[:-4] + 'jpg'
        
        fig2 = plt.figure(figsize=(10, 10))
        ax2 = fig2.add_subplot(111, projection=alignedMap)

        fl_nm2 = algn_imgs + f'/{fltr}/' + os.path.basename(files[i])[:-4] + 'jpg'
        print(f'Saving {fl_nm2}')
        alignedMap.draw_contours(axes=ax2, levels=th_lvs, zorder=1, colors=['red', 'yellow', 'blue'])
        alignedMap.plot(cmap='gray', vmin=0, vmax=Tmax)
        plt.draw()
        plt.colorbar()
        plt.savefig(fl_nm2)
        plt.close()

    #plot_data = np.array(plot_data)
    tot_count = np.array(tot_count)
    areas_level1 = np.array(areas_level1)
    areas_level2 = np.array(areas_level2)
    areas_level3 = np.array(areas_level3)
    actual_counts_level1 = np.array(actual_counts_level1)
    actual_counts_level2 = np.array(actual_counts_level2)
    actual_counts_level3 = np.array(actual_counts_level3)

    # Save the results to a CSV file
    np.savetxt(f'{fltr}_threshold_count.csv', np.c_[dates, areas_level1, areas_level2, areas_level3, actual_counts_level1, actual_counts_level2, actual_counts_level3], delimiter=',', fmt='%s')