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

# Interpolate theoretical intensity at response curve wavelengths
interp_func = interp1d(wavelengths_theoretical, intensity_theoretical, kind='linear', fill_value="extrapolate")
interpolated_intensity = interp_func(wavelengths_response)

print('wavelength bin: ', np.diff(wavelengths_response))
plt.plot(np.diff(wavelengths_response))
plt.title('Trnasmission profile wavelgth bin')
plt.ylabel('Wavelength difference')
plt.xlabel('Steps')
plt.savefig('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/results/Transm_wavelength_bin.png')
plt.show()
# Calculate expected counts
expected_counts = interpolated_intensity * transmission
print('Total count: ',np.sum(expected_counts))
# Save results
output_data = np.column_stack((wavelengths_response, transmission, interpolated_intensity, expected_counts))
np.savetxt("expected_counts.csv", output_data, delimiter=",", header="Wavelength,Transmission,Interpolated_Intensity,Expected_Counts", comments='')

# Plot for visualization
fig,axs=plt.subplots(1,1, figsize=(11,5))

axs.plot(wavelengths_response, interpolated_intensity, label="Interpolated Intensity", color="red",linewidth=0.5)
axs.plot(wavelengths_theoretical, intensity_theoretical, label="Theoretical Curve", linestyle="dashed",linewidth=0.7)
ax2=axs.twinx()
ax2.plot(wavelengths_response,transmission,'k',linewidth=0.5,label='Filter transmission')
ax2.fill_between(wavelengths_response,transmission,color='gray',alpha=0.3)
axs.fill_between(wavelengths_response, interpolated_intensity,color='blue',alpha=0.2)
ax2.set_ylabel('$\%$ Transmission')
axs.set_ylabel('Intensity')
axs.set_xlabel("Wavelength")
plt.title('Theoritical spectra and Filter Transmission')
#plt.xlabel("Wavelength")
plt.figlegend(bbox_to_anchor=(0.0001, 0.35, 0.34, 0.5))
plt.savefig('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/results/area_theoritical_n_transm.png')
plt.show()

print("Expected counts saved in 'expected_counts.csv'")
