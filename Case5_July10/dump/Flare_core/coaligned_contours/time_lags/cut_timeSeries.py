import pandas as pd
import numpy as np

for flt in ['nb03', 'nb04', 'nb08']:
    # === USER INPUTS ===
    input_file = f"{flt}_contours.csv"  # your input file
    start_time = "2024-07-10T15:20:00.000"  # ISOT start time
    end_time = "2024-07-10T16:00:00.000"    # ISOT end time
    output_file = f"{flt}_flare_segment.csv"  # output file



    # === LOAD AND FILTER ===
    # Load CSV and convert 'time' to datetime
    df = pd.read_csv(input_file)
    df['filter2_file'] = pd.to_datetime(df['filter2_file'])

    # Convert input times to datetime
    start_dt = pd.to_datetime(start_time)
    end_dt = pd.to_datetime(end_time)

    # Filter based on time
    segment = df[(df['filter2_file'] >= start_dt) & (df['filter2_file'] <= end_dt)].copy()

    # Convert time to seconds relative to start time
    segment['filter2_file'] = (segment['filter2_file']-start_dt).dt.total_seconds()
    print(f"Filtered {(segment['filter2_file'])} rows based on time")
    # Assume 2nd column is intensity/counts
    intensity_col = segment.columns[1]
    segment['error'] = np.sqrt(segment[intensity_col])

    # Reorder columns
    segment = segment[['filter2_file', intensity_col, 'error']]

    # Save to file
    segment.to_csv(output_file, index=False, header=False)
    print(f"Saved segment with seconds and Poisson errors to {output_file}")

