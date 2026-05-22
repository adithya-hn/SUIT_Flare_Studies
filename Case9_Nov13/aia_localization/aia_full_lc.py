
#--------------------------------------------------------------------------------

#Date created: 2025-07-10
#Author: @adithya-hn
#Description: This script computes the full-disk and region-of-interest (ROI) light curves from AIA data, saving the results to a CSV file and generating plots. It also optionally creates cutout images of the ROI.

#History of changes:
#2025-07-10: Initial creation based on previous scripts.
#2025-09-23: Addded exposure time to the output csv for removing saturated images.

#--------------------------------------------------------------------------------


from glob import glob
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import Helioprojective
import pathlib
import os
import warnings
warnings.simplefilter('ignore')
import gc

# ----------------------------
# CONFIGURABLE PARAMETERS
# ----------------------------

File_Source = '/media/adithya/Adi_disk4/SUIT_flare_work/case9_nov13/data/aia/aia_fd_data/1'
tx1, tx2 = -310, 220
ty1, ty2 = -450, -10
make_cutouts = False

for channel in ['131']:
    # File discovery
    files = sorted(glob(f"{File_Source}/aia.lev1_*{channel}.image_lev1.fits"))

    print(f'No. of files found: {len(files)}')
    if len(files) ==0:
        print(f"No files found for channel {channel}. Skipping...")
        continue

    # Output folders
    cutout_dir = pathlib.Path(File_Source) / f'cut_outs/{channel}_cutouts'
    cutout_img_dir = pathlib.Path(File_Source) / f'cut_outs/{channel}_cutout_imgs'
    cutout_dir.mkdir(parents=True, exist_ok=True)
    cutout_img_dir.mkdir(parents=True, exist_ok=True)

    # Output containers
    times = []
    intensities_full = []
    intensities_roi = []
    intensities_outside_roi = []
    exposure=[]

    bad_count = 0

    out_csv = f'aia_{channel}_all_lc.csv'

    # Write header only if file doesn't exist
    if not os.path.exists(out_csv):
        with open(out_csv, 'w') as f:
            f.write("Date,FullDiskIntensity,ROI_Intensity,OutsideROI_Intensity,Exposure\n")


    # ----------------------------
    # MAIN LOOP
    # ----------------------------
    yy, xx = np.mgrid[0:4096, 0:4096]
    for file in tqdm(files, desc="Processing files"):
        try:
            m = sunpy.map.Map(file)
            fnm = os.path.basename(file)

            # Define observer-correct frame
            hpc_frame = Helioprojective(observer=m.observer_coordinate, obstime=m.date)
            cutout_coords = SkyCoord(Tx=(tx1, tx2) * u.arcsec, Ty=(ty1, ty2) * u.arcsec, frame=hpc_frame)

            # Radius mask (within 1.1 solar radii)
            cenX, cenY = int(m.meta["CRPIX1"]), int(m.meta["CRPIX2"])
            r_sun = int(m.meta['R_SUN'])
            r = np.sqrt((xx - cenX)**2 + (yy - cenY)**2)
            mask_95 = r < (1.1 * r_sun)

            # Full-disk data within 1.1 R_sun
            full_disk_data = np.where(mask_95, m.data, np.nan)
            intensity_full = np.nansum(full_disk_data) / m.exposure_time.value

            # Submap = ROI
            suit_box = m.submap(cutout_coords)
            roi_intensity = np.nansum(suit_box.data) / m.exposure_time.value

            # Mask ROI region from full disk
            bottom_left = suit_box.bottom_left_coord.to_pixel(m.wcs)
            top_right = suit_box.top_right_coord.to_pixel(m.wcs)
            x1, y1 = np.round(bottom_left).astype(int)
            x2, y2 = np.round(top_right).astype(int)

            # Clip to image bounds
            x1, x2 = np.clip([x1, x2], 0, m.data.shape[1] - 1)
            y1, y2 = np.clip([y1, y2], 0, m.data.shape[0] - 1)

            # Remove ROI from full disk
            full_disk_data_masked = full_disk_data.copy()
            full_disk_data_masked[y1:y2, x1:x2] = np.nan
            intensity_outside = np.nansum(full_disk_data_masked) / m.exposure_time.value
            
 
            # Store results
            '''
            expo= m.exposure_time.value
            dt=m.date.datetime
            times.append(m.date.datetime)
            intensities_full.append(intensity_full)
            intensities_roi.append(roi_intensity)
            intensities_outside_roi.append(intensity_outside)
            exposure.append( m.exposure_time.value)#np.nansum(m.data[0:100,0:100]))
            '''


            # Optional: Save cutout and image
            if make_cutouts:
                fig = plt.figure(figsize=(5, 5))
                ax = fig.add_subplot(projection=m)
                m.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
                m.draw_quadrangle(cutout_coords, axes=ax, edgecolor="red", linewidth=2)
                plt.savefig(cutout_img_dir / f'{fnm[:-5]}.png', dpi=300)
                plt.close()
                suit_box.save(cutout_dir / fnm, overwrite=True)
            with open(out_csv, 'a') as f:
                f.write(f"{m.date.isot},{intensity_full},{roi_intensity},{intensity_outside},{m.exposure_time.value}\n")

            # Free memory
            del m, full_disk_data, suit_box, mask_95
            gc.collect()

        except Exception as e:
            print(f"Could not process {os.path.basename(file)}: {e}")
            bad_count += 1
            continue

    # ----------------------------
    # PLOT LIGHT CURVES
    # ----------------------------
    plt.figure(figsize=(10, 5))
    plt.plot(times, intensities_full, label='Full Disk', marker='o', markersize=0.5)
    plt.plot(times, intensities_roi, label='ROI Only', marker='s', markersize=0.5)
    plt.plot(times, intensities_outside_roi, label='Outside ROI', marker='x', markersize=0.5)
    plt.gcf().autofmt_xdate()
    plt.xlabel("Time")
    plt.ylabel("Intensity (DN/s)")
    plt.title(f"AIA {channel} Å Light Curves")
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    #plt.tight_layout()
    plt.savefig(f'AIA_{channel}_fd_roi_lc.png', dpi=300)
    plt.close()

    print("Done. Images with issues:", bad_count)
