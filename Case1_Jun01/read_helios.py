
import pandas as pd
data = pd.read_pickle("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/helios/for_adithya/cdte_corr_spgm_2024_06_01_07_33_to_2024_06_01_09_33.pkl")
#print(data.columns)

from astropy.io import fits

data2=fits.open('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/data/helios/2024/07/10/HLS_20240710_000051_43141sec_lev1_V211/cdte/hel1os_cdte_spectra_cdte1.fits')

print(data2[1].header)




