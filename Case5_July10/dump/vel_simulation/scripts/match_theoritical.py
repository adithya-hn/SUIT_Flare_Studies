import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Load response curve (wavelength vs transmission)
response_data = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/processed/Tot_NB08_transm.csv", delimiter=',')  # Columns: wavelength, transmission
wavelengths_response = response_data[:, 0]
transmission = response_data[:, 1]

# Load theoretical curve (wavelength vs intensity)
theoretical_data = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/raw/CaII_H_spectra.txt",skiprows=1,delimiter=' ')  # Columns: wavelength, intensity
wavelengths_theoretical = theoretical_data[:, 0]
intensity_theoretical = theoretical_data[:, 1]

# Find indices where response curve wavelengths match theoretical wavelengths
matching_indices = np.isin(wavelengths_theoretical, wavelengths_response)
print(matching_indices)
# Extract matching theoretical intensities
matched_intensities = intensity_theoretical[matching_indices]

# Ensure order of intensities matches response wavelengths
matched_intensities = np.array([intensity_theoretical[np.where(wavelengths_theoretical == wl)][0] if wl in wavelengths_theoretical else 0 for wl in wavelengths_response])

# Calculate expected counts
expected_counts = matched_intensities * transmission

#plt.plot(wavelengths_response,transmission)
plt.plot(wavelengths_response,matched_intensities)
plt.show()

# Save results
output_data = np.column_stack((wavelengths_response, transmission, matched_intensities, expected_counts))
np.savetxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/processed/expected_counts.csv", output_data, delimiter=",", header="Wavelength,Transmission,Matched_Intensity,Expected_Counts", comments='')

print("Expected counts saved in 'expected_counts.csv'")
