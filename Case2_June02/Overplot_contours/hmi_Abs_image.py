import os
import matplotlib
matplotlib.use('Agg')
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
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
#import ImagesToMovie_pkg
import matplotlib.image as mpimg
from PIL import Image
import pandas as pd
from scipy.ndimage import shift
from scipy.ndimage import gaussian_filter
from sunpy.time import parse_time
from astropy.time import Time
from sunpy.map import get_observer_meta
from sunpy.coordinates import frames, get_horizons_coord
from tqdm import tqdm
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from sunpy.visualization.animator import MapSequenceAnimator
from sunpy.coordinates import Helioprojective, SphericalScreen, propagate_with_solar_surface



hmi_map=sunpy.map.Map('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case2_June02/data/hmi/HMI/HMI_cutouts/hmi.m_45s.20240602_023000_TAI.2.magnetogram.fits')
hmi_map_data=abs(hmi_map.data)

abs_hmi_map=sunpy.map.Map(hmi_map_data,hmi_map.meta)
abs_hmi_map.save('Abs_hmi.m_45s.20240602_023000_TAI.2.magnetogram.fits')