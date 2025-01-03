

import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import glob
import datetime
from sunkit_image.coalignment import calculate_match_template_shift, apply_shifts
from datetime import timedelta
import timeit
import pathlib
from astropy.coordinates import SkyCoord, SkyOffsetFrame
import numpy as np
from PIL import Image


from sys import path as sys_path

sys_path.append('/home/adithya/Adithya_repos')
import suitkit

map = sunpy.map.Map('/Analysis/Research_Projects/SUIT_work/SUT_T24_0785_000404_Lev1.0_2024-06-08T02.31.39.797_0972NB03.fits.fits')
#save_jpg(map, 'output_image.jpg')
suitkit.view_map(map)