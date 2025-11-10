import warnings

import matplotlib.pyplot as plt
import numpy as np
from numpy.exceptions import VisibleDeprecationWarning
from parfive import Downloader

from astropy.time import Time

from sunkit_spex.extern.rhessi import RhessiLoader
from sunkit_spex.legacy.fitting.fitter import Fitter

warnings.filterwarnings("ignore", category=RuntimeWarning)
try:
    warnings.filterwarnings("ignore", category=VisibleDeprecationWarning)
except AttributeError:
    warnings.filterwarnings("ignore", category=np.exceptions.VisibleDeprecationWarning)

'''dl = Downloader()
base_url = "https://sky.dias.ie/public.php/dav/files/BHW6y6aXiGGosM6/rhessi/"
file_names = ["20021005_103800_spec.fits", "20021005_103800_srm.fits"]

#https://umbra.nascom.nasa.gov/rhessi/hessidata/2002/10/

for fname in file_names:
    dl.enqueue_file(base_url + fname, path="./rhessi/")
files = dl.download()
'''

from sunpy.net import Fido, attrs as a
results = Fido.search(a.Time("2014/1/7", "2014/1/8"), a.Instrument.rhessi, a.Physobs.spe)
Path= '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/RHESSI_Flares/Case1/data'
files = Fido.fetch(results,path=Path)
#rhess_spec = RhessiLoader(spectrum_fn="./rhessi/hsi_20140107_092400_004.fits", srm_fn="./rhessi/20021005_103800_srm.fits")

time_profile_size = (9, 6)
spec_plot_size = (6, 6)
tol = 1e-20
spec_font_size = 18
default_text = 10
#xlims, ylims = [3, 100], [5e-4, 1e4]
plt.rcParams["font.size"] = spec_font_size
plt.figure(figsize=time_profile_size)

# the line that actually plots
#rhess_spec.lightcurve(energy_ranges=[[5, 10], [10, 30], [25, 50]])

plt.show()
plt.rcParams["font.size"] = default_text