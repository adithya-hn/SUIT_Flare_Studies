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

start = timeit.default_timer()

# Parameters
csv_path = 'Flare_files_Jul11_M1.2_case25f.dat'  # <- change this to your actual file path
fol_nm = os.getcwd() + '/lc_images/'
Filters = ['NB03', 'NB04', 'NB08', 'BB01', 'BB02', 'BB03']

# ROI and Background box coordinates
#cTx1, cTy1, cTx2, cTy2 = -413, -178, -285, -93
#Tx_er1, Ty_er1, Tx_er2, Ty_er2 = -174, -85, -100, -11

cTx1=66
cTy1=-383
cTx2=191
cTy2=-200

Tx_er1=157
Ty_er1=-131
Tx_er2=235
Ty_er2=-52

# Read CSV of image paths
with open(csv_path, 'r') as f:
    all_files = [line.strip() for line in f if line.strip().endswith('.fits')]

# Loop over filters
for fltr in Filters:
    # Filter files for the current filter
    files = sorted([f for f in all_files if f'*3{fltr}.fits' in f or f.endswith(f'3{fltr}.fits')])

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

    Sequence = sunpy.map.Map(files, sequence=True)
    aligned_maps = Sequence

    for i, suit_map in enumerate(Sequence):
        fnm = os.path.basename(files[i])[:-5]
        F_name = f'{fol_nm}/{fltr}/{fnm}.jpg'
        Box_fnm = f'{box_pth}/{fnm}.jpg'
        Headr_data = Sequence[0].fits_header
        Headr_data['DATE-OBS'] = str(aligned_maps[i].date)

        # Intensity scaling
        Imn, Imx = 1000, 30000
        if fltr == 'NB08': Imx = 15000
        elif fltr == 'BB01': Imx = 21000
        elif fltr in ['BB02', 'BB03']: Imx = 26000

        # Plot full map
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection=suit_map)
        suit_map.plot(axes=ax, vmin=Imn, vmax=Imx)
        coords = SkyCoord(Tx=(cTx1, cTx2) * u.arcsec, Ty=(cTy1, cTy2) * u.arcsec, frame=suit_map.coordinate_frame)
        er_coords = SkyCoord(Tx=(Tx_er1, Tx_er2) * u.arcsec, Ty=(Ty_er1, Ty_er2) * u.arcsec, frame=suit_map.coordinate_frame)
        suit_map.draw_quadrangle(coords, axes=ax, edgecolor="red", linestyle="-", linewidth=2, label='Region of interest')
        suit_map.draw_quadrangle(er_coords, axes=ax, edgecolor="blue", linestyle="-", linewidth=2, label='Background')
        plt.colorbar()
        plt.savefig(F_name, dpi=300)
        plt.close()

        # Plot submaps
        suit_box = suit_map.submap(coords)
        er_box = suit_map.submap(er_coords)
        er_area = np.prod(er_box.data.shape)
        L, H = suit_box.data.shape
        bx_area.append(L * H)

        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(projection=suit_box)
        suit_box.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
        plt.savefig(Box_fnm, dpi=300)
        plt.close()

        # Light curve values
        exposure = Sequence[i].meta.get('CMD_EXPT')
        fltr_count.append(np.sum(suit_box.data * 1000 / exposure))
        fltr_count_err.append(np.sum(er_box.data * 1000 / exposure) / er_area)
        date_array.append(Sequence[i].date)

    # Save light curve
    np.savetxt(f'{fltr}_c25_lc_data.csv',
               np.c_[date_array, fltr_count, fltr_count_err, bx_area],
               delimiter=',', fmt='%s')

stop = timeit.default_timer()
print('Run Time: ', (stop - start)/60, 'Mins')
