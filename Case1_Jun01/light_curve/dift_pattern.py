import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()


# Load the catalogue
catalogue_path = "/Analysis/Research_Projects/SUIT_work/SUIT_Database/SUIT_Data_Jun_catalogue_2024.csv"  # your monthly CSV file
df = pd.read_csv(catalogue_path)

# Convert 'time' column to datetime
df['time'] = pd.to_datetime(df['time'])

# Define your desired time range
start_time = pd.Timestamp("2024-06-01 07:00:00")
end_time   = pd.Timestamp("2024-06-01 09:30:00")

# Query: Select only 2k images (NAXIS1 == 2048) within the time range
query_df = df[(df['NAXIS1'] == 2048) & (df['time'] >= start_time) & (df['time'] <= end_time)]

# Extract CRPIX values
crpix_values = query_df[['time', 'crpix1', 'crpix2']]

print(crpix_values)
crpix_values.to_csv("case1_crpix_values_2k.csv", index=False)

fig=plt.figure(figsize=(12,8))
plt.title('2K Drift pattern')
plt.plot(query_df['time'],query_df['crpix1'])
plt.ylabel('CRPIX1 value')
plt.xlabel('Time')
plt.ylim(1270,1290)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.savefig('dift_pattern')
plt.show()
