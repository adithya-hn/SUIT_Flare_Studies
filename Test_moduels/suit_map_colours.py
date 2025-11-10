import sunpy.map
import os
from astropy.io import fits
from colormap import filterColor
import matplotlib.pyplot as plt
import numpy as np

file='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Flare_stat/case12/data/raw/SUT_T24_1592_000646_Lev1.0_2024-11-08T02.59.26.996_0973NB04.fits'

image=sunpy.map.Map(file)
image.peek(cmap=filterColor['NB04'])

#Astropy
suit_img =fits.open(file)
suit_data = suit_img[0].data
suit_header = suit_img[0].header

plt.imshow(suit_data, cmap=filterColor['NB04'])
plt.colorbar()
plt.title("SUI Image")
plt.savefig("suit_astropy_img.png")
plt.show()