
import os
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.units as u
from sunpy.map import Map
from sunpy.map import MapSequence
from tqdm import tqdm
import glob
import datetime
import numpy as np


suit_raw_files= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/raw/'
fltr_fl = glob.glob(suit_raw_files + '*3NB04.fits')
fltr_fl=sorted(fltr_fl, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))

for i in tqdm(range(len(fltr_fl))):
    suit_map=Map(fltr_fl[i])
    intial_count=np.sum(suit_map.data)
    suit_map_rot=suit_map.rotate(angle=suit_map.meta["CROTA2"] * u.deg,missing=0)
    after_rot=np.sum(suit_map_rot.data)
    suit_map_=suit_map.rotate(angle=-suit_map.meta["CROTA2"] * u.deg,missing=0)
    back_to_same=np.sum(suit_map_.data)
    print(suit_map.meta["CROTA2"])
    print('First rotation loss:',(intial_count-after_rot)*100/intial_count, 'Final loss= ',(intial_count-back_to_same)*100/intial_count )