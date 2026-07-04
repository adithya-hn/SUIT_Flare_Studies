# Description: This script merges two text files with the same first column into a single CSV file.
#

import numpy as np

# Load the two text files
data1 = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/raw/NB8_1_inband.txt")
data2 = np.loadtxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/raw/NB8_1_oob.txt")

# Ensure first columns are identical
if not np.array_equal(data1[:, 0], data2[:, 0]):
    raise ValueError("First columns do not match exactly!")


# Convert second columns to float for processing
col2_file1 = data1[:, 1].astype(float)
col2_file2 = data2[:, 1].astype(float)

# Set negative values to zero
col2_file1[col2_file1 < 0] = 0
col2_file2[col2_file2 < 0] = 0

tot_transmission = col2_file1 * col2_file2
# Stack the first column and modified second columns side by side
merged_data = np.column_stack((data1[:, 0], col2_file1, col2_file2))

# Save as CSV
np.savetxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/processed/NB08_1_2_trnsm.csv", merged_data, fmt="%s", delimiter=",")
np.savetxt("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/vel_simulation/data/processed/NB08_tot_trnsm.csv", np.c_[data1[:, 0],tot_transmission], fmt="%s", delimiter=",")

print("Merged file saved as merged_output.csv")
