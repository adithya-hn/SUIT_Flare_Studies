import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from plots_styl import set_pub_style
set_pub_style()


# Load theoretical curve (wavelength vs intensity)
theoretical_data = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/raw/CaII_H_spectra.txt",skiprows=1,delimiter=' ')  # Columns: wavelength, intensity
wavelengths_theoretical = theoretical_data[:, 0]
intensity_theoretical = theoretical_data[:, 1]


#plt.plot(wavelengths_response,transmission)
plt.plot(wavelengths_theoretical,intensity_theoretical)
plt.ylabel('Intensity')
plt.xlabel('Wavelength')
plt.title('Ca-II h theoritical spectra for $\mu$=0')
plt.savefig('Theoritical spectra')
plt.show()
