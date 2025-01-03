import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from sunpy.net import Fido
import glob
import datetime
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
from sunkit_image.coalignment import calculate_match_template_shift,apply_shifts
from datetime import timedelta
import timeit
import pathlib
from astropy.coordinates import SkyCoord
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
#import ImagesToMovie_pkg
import matplotlib.image as mpimg
from PIL import Image
import pandas as pd
from scipy.ndimage import shift
from scipy.ndimage import gaussian_filter


fl1='/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Flare_data/SUT_T24_0785_000396_Lev1.0_2024-06-01T07.55.39.499_0973NB03.fits'
fl2='/Analysis/Projects_Data/Flare_Data/June01_Flare_Data/Flare_data/SUT_T24_0785_000396_Lev1.0_2024-06-01T07.54.14.730_0973NB03.fits'

map1=sunpy.map.Map(fl1)
map2=sunpy.map.Map(fl2)
map1.peek()
map2.peek()
Sequence = sunpy.map.Map([fl2,fl1], sequence=True)
align_shift = mc_coalign(Sequence, layer_index=0)
coords = SkyCoord(Tx=(-650, -300) * u.arcsec,Ty=(-140, -420) * u.arcsec,frame=map1.coordinate_frame,)

diff_img=gaussian_filter((Sequence[0].data/Sequence[0].meta.get('MEAS_EXP')),sigma=2)-gaussian_filter((Sequence[1].data/Sequence[1].meta.get('MEAS_EXP')),sigma=2)
diff_map=sunpy.map.Map(diff_img,Sequence[0].fits_header)
Diff_Map=diff_map.submap(coords)
fnm='7_55_diff.fits'
Diff_Map.peek()
Diff_Map.save(fnm,overwrite=True)
plt.imshow(Diff_Map.data,vmin=0,vmax=1,origin='lower')
plt.colorbar()
plt.show()
