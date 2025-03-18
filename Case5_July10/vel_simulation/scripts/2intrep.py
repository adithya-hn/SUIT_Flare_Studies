import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from plots_styl import set_pub_style
set_pub_style()


# Load response curve (wavelength vs transmission)
response_data = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/processed/Shifted_transm_profile.csv", delimiter=',')  # Columns: wavelength, transmission
wavelengths_response = response_data[:, 0]
transmission = response_data[:, 1]

# Load theoretical curve (wavelength vs intensity)
theoretical_data = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/raw/CaII_H_spectra.txt",skiprows=1,delimiter=' ')  # Columns: wavelength, intensity
wavelengths_theoretical = theoretical_data[:, 0]
intensity_theoretical = theoretical_data[:, 1]
wavln=np.arange(3950000,3985000,1)/10000
print(wavln)
# Interpolate theoretical intensity at response curve wavelengths
interp_func1 = interp1d(wavelengths_theoretical, intensity_theoretical, kind='linear', fill_value="extrapolate")
interp_func2 = interp1d(wavelengths_response, transmission, kind='linear', fill_value="extrapolate")

interp_intensity = interp_func1(wavln)
interp_intensity_=np.copy(interp_intensity)
theo_intrep_int=np.copy(interp_intensity)
intrep_transm=interp_func2(wavln)
tot_count=[]
shifted_array_b=[]
shifted_array_r=[]

shifted_int_=np.concatenate((interp_intensity_,np.zeros(1000)))
for i in range (2000):
    
    if i<1000:

        shifted_intensity_=shifted_int_[(1000-i):len(shifted_int_)-i]
        expected_counts = shifted_intensity_ * intrep_transm
        tot_count.append(np.sum(expected_counts))
    else:

        shifted_int=np.insert(interp_intensity,0,0)
        shifted_intensity=shifted_int[:-(i+1-1000)]
        expected_counts = shifted_intensity * intrep_transm
        tot_count.append(np.sum(expected_counts))
        interp_intensity=shifted_int
    if i == 891:
        shifted_array_b.append(shifted_intensity_)
        print(shifted_array_b[0])

    
shifted_array=np.array(shifted_array_b)
tot_count=np.array(tot_count)
print(min(tot_count),np.where(tot_count==min(tot_count)))
steps=np.arange(2000)
plt.plot(steps,tot_count)

plt.xlabel('steps')
plt.ylabel('Counts')
plt.title('red blue shift count')
plt.savefig('step_counts.png')
plt.close()

fig,axs=plt.subplots(1,1, figsize=(11,5))
ax2=axs.twinx()
axs.plot(wavln,shifted_array[0],label='Shifted spectra at max count')
axs.plot(wavln,theo_intrep_int,label='Theoritecal spectra of ca II h')
ax2.plot(wavln,intrep_transm,'k',label='Transmission profile')
axs.set_xlabel('wavelength')
ax2.set_ylabel('$\%$ Transmission')
axs.set_ylabel('Intensity')
plt.title('Red and Blue shift')
plt.figlegend(bbox_to_anchor=(0.0001, 0.35, 0.34, 0.5))
plt.savefig('shift_count.png')
plt.close()

