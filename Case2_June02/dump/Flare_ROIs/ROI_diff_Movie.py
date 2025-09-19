import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import datetime
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
from scipy.ndimage import gaussian_filter
from astropy.coordinates import SkyCoord
import numpy as np
import pathlib
import map_to_image
from colormap import filterColor

def analyze_image_difference(dat_file, Filter='NB03', d=2):
    """Analyze differences in FITS images over time.

    Parameters:
        dat_file (str): The path to the .dat file containing the list of FITS files.
        Filter (str): The filter to be applied. Default is 'NB08'.
        d (int): The difference interval. Default is 2.

    Returns:
        None
    """
    print(f'Searching for {Filter} images in the provided .dat file')

    # Read the .dat file to get the list of FITS files
    with open(dat_file, 'r') as file:
        files = file.read().splitlines()

    # Ensure the files are sorted based on the timestamp in the filename
    print(files[0])
    files = sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    print('Total files:', len(files))

    if len(files) < d:
        print('Not enough files to analyze differences.')
        return

    ref_img = sunpy.map.Map(files[0])
    
    fl_date = datetime.datetime.fromisoformat(str(ref_img.date))
    fol_nm = f'{fl_date.day:02d}_{fl_date.month:02d}_{fl_date.year:02d}'

    jpg_fold = os.path.join(fol_nm, 'ROI_Diff_imgs', Filter)
    pathlib.Path(jpg_fold).mkdir(parents=True, exist_ok=True)

    sigma = 1
    for i in range(len(files) - d):
        suitMap = sunpy.map.Map(files[i])
        suitMap_ = sunpy.map.Map(files[i + d])
        
        fl_pth = os.path.join('Flare_imgs', Filter, f'{os.path.basename(files[i + d])[:-4]}.jpg')
        pathlib.Path(os.path.dirname(fl_pth)).mkdir(parents=True, exist_ok=True)
        
        Sequence = sunpy.map.Map([files[i], files[i + d]], sequence=True)
        align_sq = mc_coalign(Sequence, layer_index=1, clip=False)
        coords =  SkyCoord(Tx=(-370, -250) * u.arcsec, Ty=(-230, -330) * u.arcsec, frame=suitMap.coordinate_frame)
        
        diff_img = gaussian_filter((align_sq[1].data / align_sq[1].meta.get('MEAS_EXP')), sigma=2) - gaussian_filter((align_sq[0].data / align_sq[0].meta.get('MEAS_EXP')), sigma=2)
        map_to_image.sunpy_map_to_jpg(align_sq[0], fl_pth)
        
        diff_map = sunpy.map.Map(diff_img, Sequence[1].fits_header)
        Diff_Map = diff_map.submap(coords)
        fl_nm = os.path.join(jpg_fold, f'{os.path.basename(files[i + d])[:-4]}.jpg')
        
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection=Diff_Map)
        Diff_Map.plot(cmap=filterColor[Filter], vmin=0, vmax=8)
        ax.plot_coord(SkyCoord(-320 * u.arcsec, -310 * u.arcsec, frame=Diff_Map.coordinate_frame), color='white', marker='X', markersize=12)
        plot_str = f'{Diff_Map.date} - {suitMap.date}'
        ax.text(50, 50, plot_str, color='white', fontsize=10)
        plt.xticks([])
        plt.yticks([])
        plt.tight_layout()
        plt.savefig(fl_nm)
        plt.close()

        print(i, ' / ', len(files))

    print('Image difference analysis completed.')

# Example usage
analyze_image_difference('Flare_files_Jun2_M1.2.dat', Filter='NB08', d=2)
