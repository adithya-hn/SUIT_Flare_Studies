import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from plots_styl import set_pub_style
set_pub_style()


# Load response curve (wavelength vs transmission)
response_data = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/processed/Tot_NB08_transm.csv", delimiter=',')  # Columns: wavelength, transmission
wavelengths_response = response_data[:, 0]
transmission = response_data[:, 1]

# Load theoretical curve (wavelength vs intensity)
theoretical_data = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/raw/CaII_H_spectra.txt",skiprows=1,delimiter=' ')  # Columns: wavelength, intensity
wavelengths_theoretical = theoretical_data[:, 0]
intensity_theoretical = theoretical_data[:, 1]
wavln=np.arange(39500,39850,1)/100
print(wavln)
# Interpolate theoretical intensity at response curve wavelengths
interp_func1 = interp1d(wavelengths_theoretical, intensity_theoretical, kind='linear', fill_value="extrapolate")
interp_func2 = interp1d(wavelengths_response, transmission, kind='linear', fill_value="extrapolate")

interp_intensity = interp_func1(wavln)
intrep_transm=interp_func2(wavln)
tot_count=[]
for i in range (250):
    shifted_int=np.insert(interp_intensity,0,0)
    shifted_intensity=shifted_int[:-(i+1)]
    print(len(interp_intensity),len(shifted_intensity),len(shifted_int))
    expected_counts = shifted_intensity * intrep_transm
    tot_count.append(np.sum(expected_counts))
    print('Total counts: ',np.sum(expected_counts))
    interp_intensity=shifted_int
    

steps=np.arange(250)
plt.plot(steps,tot_count)
plt.show()

