import numpy as np
from scipy.interpolate import interp2d
import math
from PIL import Image
import math
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sunpy.map
import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
import timeit
import datetime


suit_map = sunpy.map.Map('/Analysis/SUIT_data/Flare_data/SUT_T24_0785_000396_Lev1.0_2024-06-01T08.49.22.903_0983NB03.fits')
#suit_map.peek()
'''center_coord = SkyCoord(-450 * u.arcsec, -250 * u.arcsec, frame=suit_map.coordinate_frame) #54,157
width = 170 * u.arcsec
height =170 * u.arcsec
rectangle = SkyCoord(lon=[-1/2, 1/2] * width, lat=[-1/2, 1/2] * height, frame=center_coord)'''

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(projection=suit_map)
suit_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
coords = SkyCoord(
    Tx=(-520, -360) * u.arcsec,
    Ty=(-160, -320) * u.arcsec,
    frame=suit_map.coordinate_frame,)
suit_map.draw_quadrangle(
    coords,
    axes=ax,
    edgecolor="blue",
    linestyle="-",
    linewidth=2,
    label='2-element SkyCoord')

box1=suit_map.submap(coords)
box1.peek()

