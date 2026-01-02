import os
import numpy as np
import astropy.units as u
from sunpy.map import Map
#from sunpy.physics.differential_rotation import propagate_with_solar_surface
from sunpy.coordinates import Helioprojective, SphericalScreen, propagate_with_solar_surface
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
import glob
import datetime
from tqdm import tqdm

# -------------------------------
# CONFIG
# -------------------------------
NPROC = max(1, cpu_count() - 6)  # leave 1 core free
PNG_DIR = "../data/aligned_crop_jpgs"
FITS_DIR = "../data/derot_fits"
DEROT_PNG_DIR = "../data/derotated_jpgs"

os.makedirs(PNG_DIR, exist_ok=True)
os.makedirs(FITS_DIR, exist_ok=True)
os.makedirs(DEROT_PNG_DIR, exist_ok=True)

# # Precomputed once
# x1, y1, x2, y2 = (
#     max(blX) + 30,
#     max(blY) + 30,
#     min(trX) - 30,
#     min(trY) - 30
# )

# -------------------------------
# WORKER FUNCTION
# -------------------------------
def process_one_file(fits_path):
    
    suit_map = Map(fits_path)
    #print('ok')

    with propagate_with_solar_surface():
        smap = suit_map.reproject_to(ref_baseMap.wcs)

    smap.meta['f_name']=suit_map.meta['f_name']
    fname = smap.meta["F_NAME"]
    #print(fname)
    smap.save(os.path.join(FITS_DIR,fname),overwrite=True)

    # ---- Bounding box ----
    data = smap.data
    nz = np.argwhere(data > 100)
    ymin, xmin = nz.min(axis=0)
    ymax, xmax = nz.max(axis=0)

    # ---- Derotated PNG ----
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection=ref_baseMap)
    smap.plot(axes=ax)
    ax.scatter([xmin], [ymin], c="red", s=10)
    ax.scatter([xmax], [ymax], c="blue", s=10)
    plt.savefig(
        os.path.join(DEROT_PNG_DIR, fname[:-5] + ".jpg"),
        dpi=200)
    plt.close(fig)

    return {
        "file": fname,
        "xmin": xmin,
        "ymin": ymin,
        "xmax": xmax,
        "ymax": ymax,
        "date": suit_map.date.datetime
    }

    
    

if __name__ == "__main__":

    suit_raw_files= '../data/aligned_fits/'
    fltr='NB04'
    fltr_fl = glob.glob(suit_raw_files + '*3'+f'{fltr}.fits')
    fltr_fl=sorted(fltr_fl, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    ref_baseMap=Map(fltr_fl[1146])
    print('Reference derot frame: ',ref_baseMap.date) 
    # for i in range(10):
    #     res=process_one_file(fltr_fl[i])




    with Pool(NPROC) as pool:
        results = list(tqdm( pool.imap_unordered(process_one_file, fltr_fl), total=len(fltr_fl)))
        print(results)
    import pandas as pd

    df = pd.DataFrame(results)
    df.to_csv('test.csv',index=False)

    # Separate failures if needed
    df_ok = df[~df["file"].str.contains("error", na=False)]
    df_err = df[df.columns.intersection(["file", "error"])].dropna(how="all")

    # Write outputs
    df_ok.to_csv("../data/bbox_results.csv", index=False)
    df_err.to_csv("../data/bbox_errors.csv", index=False)
    #len(fltr_fl)
    # from joblib import Parallel, delayed

    # Parallel(n_jobs=NPROC, backend="loky")(
    # delayed(process_one_file)(f) for f in fltr_fl)