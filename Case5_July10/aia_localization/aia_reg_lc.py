"""
Created on Oct 9 2025
@Author: Adithya-hn

- to cretae the olight curve of the aia active regions
"""

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
from pathlib import Path
import pathlib
from astropy.coordinates import SkyCoord
#from colormap import filterColor
import numpy as np
from sunpy.coordinates import RotatedSunFrame
from sunpy.physics.differential_rotation import solar_rotate_coordinate
from sunpy.coordinates import Helioprojective
from sunpy.map.maputils import coordinate_is_on_solar_disk
from sunpy.time import parse_time
from tqdm import tqdm
import pandas as pd
import imageio.v3
import warnings
warnings.simplefilter('ignore')



#--------------------
data_folder=Path('/media/adithya/Adi_disk4/SUIT_flare_work/case5_jul10/data/aia/aia_fd')
regions_csv='aia131_regions.csv'
out_csv='aia131_region_lc.csv'
out_dir='AR_patches'


#--------------------
regions_df = pd.read_csv(regions_csv)
#out_dir.mkdir(exist_ok=True, parents=True)
fits_files = sorted(data_folder.glob("*.131.image_lev1.fits"))
start_time =sunpy.map.Map(fits_files[0]).date #parse_time('2024-06-01 06:58:56')
print('Total number AIA files: ', len(fits_files))

records = []
exposures=[]
for f in tqdm(fits_files):
    m = sunpy.map.Map(f)
    row_data = {"time": m.date.datetime}
    exp = m.meta.get("exptime", np.nan)
    exposures.append(exp)

    #for region_id in regions_df["region_id"].unique():
    for _, reg in regions_df.iterrows():
        try:
            region_id = int(reg["region_id"])
            reg = regions_df[regions_df["region_id"] == region_id].iloc[0]
            width=int(reg['width_pix'])
            height=int(reg['height_pix'])
            ctx=int(reg['Centre_hpc_Tx'])
            cty=int(reg['Centre_hpc_Ty'])
            strt_time_str=str(reg['file'])
            
            #print(strt_time_str[17:34])
            dist=out_dir+f'/{region_id}'
            Path(dist).mkdir(parents=True, exist_ok=True)
            center_coord=SkyCoord(ctx*u.arcsec, cty*u.arcsec, obstime=start_time, observer="earth", frame=Helioprojective)
            inside = coordinate_is_on_solar_disk(center_coord)
            if inside==False:
                print(f'Region id : {region_id} is outside the disk, off setting both x and y by 10')
                if ctx<0:
                    ctx=ctx+10
                if ctx>0:
                    ctx=ctx-10
                if cty<0:
                    cty=cty+10
                if cty>0:
                    cty=cty-10
                center_coord=SkyCoord(ctx*u.arcsec, cty*u.arcsec, obstime=start_time, observer="earth", frame=Helioprojective)

            rot_coord = solar_rotate_coordinate(center_coord,time=m.date)
            px = m.world_to_pixel(rot_coord)
            cx, cy = int(px.x.value), int(px.y.value)
            #print(cx,cy)
            # crop a fixed-size box around this shifted center
            x1 = cx - width//2 #floor divison
            y1 = cy - height//2
            x2 = cx + width//2 
            y2 = cy + height//2
            sub = m.data[y1:y2, x1:x2]
            sub8=sub.astype(np.uint8)
            
            imageio.v3.imwrite(dist+f'/{region_id}_{m.date}.png',sub8)
            intensity = np.nansum(sub)
            #rows.append({"time": m.date.datetime, "intensity": intensity})
            row_data[f"region_{region_id}"] = intensity
        except Exception as e:
            print(f"Skipping {f.name}: {e}")
            continue

    
    records.append(row_data)

# --- Convert to DataFrame in wide format ---
df_wide = pd.DataFrame(records).sort_values("time")
df_wide.insert(1, "exposure", exposures)
df_wide.to_csv(out_csv, index=False)
print("Saved:", out_csv)
# --- Plot all curves ---
plt.figure(figsize=(10,6))
for col in df_wide.columns[2:]:
    plt.plot(df_wide["time"], df_wide[col], marker=".", label=col)
plt.xlabel("Time (UTC)")
plt.ylabel("Mean Intensity (DN)")
plt.title("AIA 131 Light Curves (all detected regions)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('All_lc.png', dpi=150)
plt.show()