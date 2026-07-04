import pandas as pd
import os
import astropy.units as u
from sunpy.map import Map
from tqdm import tqdm
import matplotlib.pyplot as plt

df = pd.read_csv("../data/bbox_results.csv")

# Intersection of all valid bounding boxes
x1 = df["xmin"].max()+30
y1 = df["ymin"].max()+30
x2 = df["xmax"].min()-30
y2 = df["ymax"].min()-30

if (x2 <= x1) or (y2 <= y1):
    raise ValueError("No common overlapping region found")

print("Common crop box (pixels):", x1, y1, x2, y2)

IN_DIR  = "../data/derot_fits"
OUT_DIR = "../data/derotated_common_crop"
PNG_DIR = "../data/aligned_crop_jpgs"
os.makedirs(OUT_DIR, exist_ok=True)

for fname in tqdm(df["file"].unique()):
    fpath = os.path.join(IN_DIR, fname)

    if not os.path.exists(fpath):
        continue

    smap = Map(fpath)

    cropped = smap.submap(
        bottom_left=[x1, y1] * u.pix,
        top_right=[x2, y2] * u.pix)

    cropped.save(
        os.path.join(OUT_DIR, fname),
        overwrite=True
    )
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection=cropped)
    cropped.plot(axes=ax)
    plt.savefig(
        os.path.join(PNG_DIR, fname[:-5] + ".jpg"),
        dpi=200)
    plt.close(fig)