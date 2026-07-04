

import numpy as np
import matplotlib.pyplot as plt

# Load the zDCF output
data = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/time_lags/helios/nb03/ccf.dcf")
tau, err_low, err_up, zdcf,decfer_low,decfer_up,n = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4],data[:, 5],data[:, 6]
# Compute asymmetric error bars
yerr = [err_low, err_up]
xerr = [decfer_low, decfer_up]

# Plotting
plt.figure(figsize=(8, 5))
plt.errorbar(tau, zdcf, yerr=yerr,xerr=xerr, fmt='o', ecolor='gray', capsize=3, label='zDCF')
plt.axhline(0, color='k', linestyle='--', lw=0.8)
plt.xlabel('Lag (s)', fontsize=12)
plt.ylabel('zDCF', fontsize=12)
plt.title('Z-transformed Discrete Correlation Function', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.5)
plt.ylim(-20,20)
plt.legend()
plt.tight_layout()

plt.show()
