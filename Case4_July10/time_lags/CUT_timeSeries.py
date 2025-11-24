import pandas as pd

# === USER INPUTS ===
input_file  = "NB03_flare.csv"  # your input file
in_helios   = "helios_CdTe_c4.csv"

start_time  = "2023-11-28T10:15:00"  # ISOT start time
end_time    = "2023-11-28T10:20:00"    # ISOT end time

output_file = "NB03_flare_segment.csv"  # output file
out_helios  = "helios_flare_segment.csv" 

# === LOAD AND FILTER ===
# Load CSV and convert 'time' to datetime
df = pd.read_csv(input_file)
df['time'] = pd.to_datetime(df['time'])

# Convert input times to datetime
start_dt = pd.to_datetime(start_time)
end_dt = pd.to_datetime(end_time)

# Slice the DataFrame
segment = df[(df['time'] >= start_dt) & (df['time'] <= end_dt)]

# Save the result
segment.to_csv(output_file, index=False)
print(f"Saved {len(segment)} rows to {output_file}")
