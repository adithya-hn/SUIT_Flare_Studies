import pandas as pd
import numpy as np

# === USER INPUTS ===
input_file  = "c4_NB04_lc_data.csv"  # your input file
in_helios   = "helios_CdTe_c4.csv"

start_time  = "2024-07-10T04:15:00"  # ISOT start time
end_time    = "2024-07-10T04:30:00"    # ISOT end time

output_file = "NB04_flare_segment_1.csv"  # output file
out_helios  = "helios_flare_segment_1.csv" 

# === LOAD AND FILTER ===
# Load CSV and convert 'time' to datetime

df = pd.read_csv(input_file)

df['Time'] = pd.to_datetime(df['Time'])
# Convert input times to datetime
start_dt = pd.to_datetime(start_time)
end_dt = pd.to_datetime(end_time)

# Slice the DataFrame
segment = df[(df['Time'] >= start_dt) & (df['Time'] <= end_dt)].copy()
segment["Time_sec"] = (segment["Time"].astype("int64") / 1e9)# seconds since 1970-01-01 , devided by 1e9 to scale down
segment["Error"]= np.sqrt(segment["Total"]*segment["Exposure"])/segment["Exposure"]
segment_out=segment[['Time_sec','Total','Error']]

# Save the result
segment_out.to_csv(output_file, index=False)
print(f"Saved {len(segment)} rows to {output_file}")


#------------------------------------------------
df = pd.read_csv(in_helios)
df['Time'] = pd.to_datetime(df['Time'])

# Convert input times to datetime
start_dt = pd.to_datetime(start_time)
end_dt = pd.to_datetime(end_time)

# Slice the DataFrame
segment = df[(df['Time'] >= start_dt) & (df['Time'] <= end_dt)].copy()
segment["Time_sec"] = (segment["Time"].astype("int64") / 1e9)# seconds since 1970-01-01 , devided by 1e9 to scale down


segment_out=segment[['Time_sec','Total','CdTe1+2Er']]
# Save the result
segment_out.to_csv(out_helios, index=False)
print(f"Saved {len(segment)} rows to {out_helios}")
