import asdf
import numpy as np
import sunpy.map as smap

# Create MapSequence
maps = smap.Map("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/aligned_crop_/*NB04.fits", sequence=True)

# Build ASDF tree
tree = {
    "data": np.stack([m.data for m in maps]),   # (nt, ny, nx)
    "meta": [dict(m.meta) for m in maps],       # per-frame metadata
}

# Write ASDF file
with asdf.AsdfFile(tree) as af:
    af.write_to("mapsequence.asdf")


with asdf.open("mapsequence.asdf") as af:
    data = af["data"]
    meta = af["meta"]

maps = smap.Map(data[0], meta[0])#, sequence=True)
maps.peek()