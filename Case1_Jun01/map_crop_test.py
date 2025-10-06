
import sunpy.map
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.wcs import WCS


def submap_with_ginput(sunpy_map):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=sunpy_map)
    sunpy_map.plot(axes=ax)
    ax.set_title("Click two corners of the ROI (bottom-left and top-right)")

    # Wait for 2 clicks
    pts = plt.ginput(2, timeout=0)
    plt.close(fig)

    if len(pts) != 2:
        raise RuntimeError("ROI selection cancelled or failed.")

    (x1, y1), (x2, y2) = pts
    bottom_left = (min(x1, x2), min(y1, y2)) * u.pixel
    top_right = (max(x1, x2), max(y1, y2)) * u.pixel
    X=(x1+x2)/2
    Y=(y1+y2)/2
    w=WCS(sunpy_map.fits_header) 
    sky = w.pixel_to_world(X,Y)
    bl=w.pixel_to_world(min(x1, x2),min(y1, y2))
    tr=w.pixel_to_world(max(x1, x2),max(y1, y2))
    print('Mid_points: ',sky.Tx,sky.Ty)
    print('Width and Height:',tr.Tx-bl.Tx,tr.Ty-bl.Ty)

    submap = sunpy_map.submap(bottom_left=bottom_left, top_right=top_right)
    return submap

img=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/1600_aligned/SUT_T24_0785_000396_Lev1.0_2024-06-01T07.34.27.945_0973NB03.fits')

croped_map=submap_with_ginput(img)
croped_map.peek()