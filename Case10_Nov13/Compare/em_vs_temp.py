

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.dates as mdates
# Replace with your actual file path and delimiter if not CSV
df = pd.read_csv('fit_results_AL1_SOLEXS_20241113_SDD2_L1_2411131500_2411131730_TEMP_EM.txt',sep='\t')  # Or use sep='\t' or delim_whitespace=True if needed

# Combine date and time into a datetime column
df['datetime'] = pd.to_datetime(df['Time'])

# Convert emission measure and temperature to numeric, if needed
df['emission_measure'] = pd.to_numeric(df['EM'], errors='coerce')
df['temperature'] = pd.to_numeric(df['Temperature'], errors='coerce')


# Define time range for the flare
start_time = '2024-11-13 16:50:00'
end_time = '2024-11-13 17:30:00'

mask = (df['datetime'] >= start_time) & (df['datetime'] <= end_time)
flare_df = df.loc[mask]

time_nums = flare_df['datetime'].map(mdates.date2num)

# Create scatter plot with color mapped to time
plt.figure(figsize=(8, 6))
sc = plt.scatter( flare_df['emission_measure'],flare_df['temperature'],
                 c=time_nums, cmap='rainbow', edgecolor='k', s=30)

# Add colorbar with time formatting
cbar = plt.colorbar(sc, label='Time (UTC)')
cbar.ax.yaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
plt.xscale('log')
plt.ylabel('Temperature (MK)')
plt.xlabel('Emission Measure ( 10$^{46}$cm$^{-3}$)')
plt.title('EM vs Temperature (2024-11-13)')
plt.grid(True)
plt.tight_layout()
plt.savefig('Em_vs_temp_c10_.png',dpi=300)
plt.show()
