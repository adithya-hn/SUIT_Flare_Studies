from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from sunpy.map import Map
from sunpy.physics.differential_rotation import solar_rotate_coordinate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec    
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib.patches import Rectangle
from astropy.wcs.utils import skycoord_to_pixel
from pathlib import Path
from sunpy.coordinates import RotatedSunFrame
from sunpy.physics.differential_rotation import solar_rotate_coordinate
from sunpy.coordinates import Helioprojective
from sunpy.map.maputils import coordinate_is_on_solar_disk
from datetime import datetime
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

start_time= Time("2024-10-31T23:30:08")
data_folder=Path('/media/adithya/Adi_disk4/SUIT_flare_work/case8_nov01/data/aia/aia_fd_data')
f='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/aia/aia.lev1_euv_12s.2024-11-01T002320Z.131.image_lev1.fits'
regions_df=pd.read_csv('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/aia_localization/aia131_regions.csv')
region_lc=''
exposures=[]
pad_pixels = 60              # pad each bbox by this many pixels
max_regions_to_plot = 20     # safety cap
props_to_use=20
m = Map(f)
row_data = {"time": m.date.datetime}
exp = m.meta.get("exptime", np.nan)
exposures.append(exp)
# fig, ax = plt.subplots(figsize=(8,8))
# ax.imshow(m.data, origin='lower',cmap='sdoaia131', vmin=np.nanpercentile(m.data, 1), vmax=np.nanpercentile(m.data, 99.8))
fig = plt.figure(figsize=(18,8))
gs  = GridSpec(1,2,
               width_ratios=[1,2],
               hspace=0)
ax = fig.add_subplot(gs[:,0])

ax.imshow(m.data, origin='lower',cmap='sdoaia131', vmin=np.nanpercentile(m.data, 1), vmax=np.nanpercentile(m.data, 99.8))
for _, reg in regions_df.iterrows():
    region_id = int(reg["region_id"])
    reg = regions_df[regions_df["region_id"] == region_id].iloc[0]
    width=int(reg['width_pix'])
    height=int(reg['height_pix'])
    ctx=int(reg['Centre_hpc_Tx'])
    cty=int(reg['Centre_hpc_Ty'])
    strt_time_str=str(reg['file'])

    center_coord=SkyCoord(ctx*u.arcsec, cty*u.arcsec, obstime=start_time, observer="earth", frame=Helioprojective)
    inside = coordinate_is_on_solar_disk(center_coord)
    if inside==False:
        print(f'Region id : {region_id} is outside the disk, off setting both x and y by 40')
        if ctx<0:
            ctx=ctx+40
        if ctx>0:
            ctx=ctx-40
        if cty<0:
            cty=cty+40
        if cty>0:
            cty=cty-40
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
    rect = Rectangle((x1, y1), width, height, linewidth=1.6, edgecolor='red', facecolor='none')
    ax.add_patch(rect)
    ax.text(x1+2, y1+10, f"R{region_id}", color='yellow', fontsize=8, bbox=dict(facecolor='black', alpha=0.4, pad=1))
# plt.show()
ax.set_title(fr'AIA 131$~\mathrm{{\AA}}$ (2024-Nov-01 00:23:18)')
ax.set_ylabel('Y(pixels)')
ax.set_xlabel('X(pixels)')
cols_to_plot = ["region_1", "region_2", "region_3","region_4","region_5","region_6","region_7","region_8","region_9","region_10","region_11"]  # change names as needed
# --- Load your AIA light curves ---
aia_df = pd.read_csv("aia131_region_lc.csv", parse_dates=["time"])
# --- Bin AIA to 30s cadence (matching HELIOS) ---
aia_df["time"] = pd.to_datetime(aia_df["time"])
aia_df = aia_df.set_index("time")
aia_filtered = aia_df[aia_df["exposure"] >= 0.1]
aia_norm = aia_filtered.copy()
region_cols = [c for c in aia_filtered.columns if c.startswith("region_")]
aia_norm[region_cols] = aia_norm[region_cols].div(aia_norm["exposure"].replace(0, np.nan), axis=0)
aia_binned = aia_norm.resample("60s").mean().reset_index()

Helios=(np.load("cdte_data_flare_7.npy", allow_pickle=True)).transpose()
datetime_objects = pd.to_datetime(Helios[0])
helios_time=helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]

helios_rate=np.array((Helios[1]+Helios[2]),dtype=float)
helios_df = pd.DataFrame({"time": helios_time, "helios_rate": helios_rate},index=helios_time)
helios_df["time"] = pd.to_datetime(helios_df["time"])
helios_df = helios_df.set_index("time")

helios_df_binned = helios_df#.resample("1min").mean()

combined_df = pd.merge_asof(aia_binned.sort_values("time"),
                   helios_df_binned.sort_values("time"),
                   on="time", direction="nearest")

ax1 = fig.add_subplot(gs[0,1])
ax1.plot(combined_df["time"], combined_df["helios_rate"], label="HELIOS", color="black", lw=2)
ax2=ax1.twinx()

for col in cols_to_plot:
    valid = combined_df[[col, "helios_rate"]].dropna()
    # r, _ = pearsonr(valid[col], valid["helios_rate"])
    ax2.plot(combined_df["time"], combined_df[col], label=f'{col}')


ax2.legend(loc='upper right',fontsize=18)
ax1.legend(loc='upper left',fontsize=18)
ax1.set_xlabel("Time (UT)")
ax1.set_ylabel("HEL1OS Count rate (counts/s)")
ax2.set_ylabel("AIA 131 Region counts (DN/s)")
ax1.set_yscale("log")
ax2.set_yscale("log")
plt.title("HELIOS vs AIA Region Light Curves (2024-Nov-01)")

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)

plt.tight_layout()
plt.savefig('c7_AIA-HEL1OS-lc.pdf',dpi=300)
plt.close()

