import numpy as np
import matplotlib.pyplot as plt
import scipy.io
# from scipy.integrate import simps
import os
kiris = scipy.io.readsav('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/IRIS_flux_calib/iris_images/iris_real_units_20241009.sav') # restoring DN to flux unit converted IRIS values
# kint = kiris['kint'] #converted intensity
# kerr = kiris['kerr'] # error in conversion
# kwvl = kiris['kwvl'] # iris wavelength array
# res=[0.698,0.698]
# kint = kint.copy()
# kwvl = kwvl.copy()
# kerr = kerr.copy()
# kint *= 1e4*1e-7 # converting ergs/cm^2 to SI unit
# kerr *= 1e4*1e-7
# kwvl /= 10. #converting angstrom to nm
print(kiris.keys())

data=kiris['iris_real_units_per_ang'][0]
print(data['SCALE'])
os._exit(1)


#reading SUIT effective area files
# wave_suit, nb3 = np.loadtxt(f'../../SUIT_eff_prof_Soumya_Feb13_2025/NB3_eff_outband_fix.txt',unpack=True)
# nb3_i = interp1d(wave_suit,nb3,fill_value='extrapolate')
wave_suit1, nb4 = np.loadtxt(f'NB4_eff_outband_fix.txt',unpack=True)
nb4_i = interp1d(wave_suit1,nb4,fill_value='extrapolate')
# nb3_area = nb3_i(kwvl)
nb4_area = nb4_i(kwvl)


data = np.reshape(kint,(kint.shape[0]*kint.shape[1],kint.shape[2]))  # converting iris intensity into array
data_err = np.reshape(kerr,(kerr.shape[0]*kerr.shape[1],kerr.shape[2]))
im_nb3 = np.zeros(data.shape[0]); overlap_region = np.zeros(data.shape[0])
im_nb4 = np.zeros(data.shape[0]); overlap_region = np.zeros(data.shape[0])

omega_pixel = (np.pi/3600/180)**2*res[0]*res[1] # We are getting rid of IRIS solid angle using it's plate scale. 

for i,j in enumerate(data):
 #Counts in DN unit
    # im_nb3[i] =  simps(j*nb3_area*omega_pixel/((cn.h*cn.c)/(kwvl*1e-9)/3.0),kwvl)  
    im_nb4[i] =  simps(j*nb4_area*omega_pixel/((cn.h*cn.c)/(kwvl*1e-9)/3.0),kwvl)
# The below final image must be in units of DN/sec/pix for IRIS which is same as of SUIT's which should be used for scaling factor
image_nb3 = np.reshape(im_nb3,(kint.shape[0],kint.shape[1]))
image_nb4 = np.reshape(im_nb4,(kint.shape[0],kint.shape[1]))




