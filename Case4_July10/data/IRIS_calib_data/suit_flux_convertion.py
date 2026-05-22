from astropy.io import fits
import numpy as np
import glob
import sunpy.map
import astropy.units as u
import os
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

file_list = sorted(glob.glob("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/aligned_crop_fits/*NB04.fits"))

flux_maps = []

for f in tqdm(file_list):
    m = sunpy.map.Map(f)
    C=2.45*10**2
    # convert to flux
    flux = m.data * C+3.45*1e6
    # fig,(ax1,ax2)=plt.subplots(1,2)
    # ax1.imshow(m.data)
    # ax2.imshow(flux)
    # ax1.set_title('raw')
    # ax2.set_title('Flux map')
    # plt.show()
    m.meta["exposure"]=1
    m.meta['CMD_EXPT']=1000
    fm=sunpy.map.Map(flux,m.meta)
    fm.save(fm.meta["f_name"])

   
