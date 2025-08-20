import pandas as pd
import numpy as np

for flt in ['NB03', 'NB04', 'NB08','flare_3']: #,
    # === USER INPUTS ===
    if flt == 'flare_3':

        numpy_array=np.load(f"cdte_data_{flt}.npy", allow_pickle=True)
        print(numpy_array.shape)
        df =pd.DataFrame(numpy_array) #, header=None, names=['time', 'intensity_cdte1','intensity_cdte2'])
        df.columns = ['time', 'intensity_cdte1','intensity_cdte2','er1','er2']
    else:
        input_file = f"{flt}_c3_lc_data.csv"  # your input file

        #/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/time_lags/_NB04_c8_lc_data.csv
        df = pd.read_csv(input_file, header=None, names=['time', 'intensity','QS_mean','AR_box_area','QS_area'])
    start_time = "2024-06-02T07:52:00.000"  # ISOT start time
    end_time   = "2024-06-02T08:22:00.000"    # ISOT end time
    output_file = f"{flt}_flare_segment.csv"  # output file



    # === LOAD AND FILTER ===
    # Load CSV and convert 'time' to datetime
    #df = pd.read_csv(input_file)
    
    df['time'] = pd.to_datetime(df['time'])

    #print(df['time'])

    # Convert input times to datetime
    start_dt = pd.to_datetime(start_time)
    end_dt = pd.to_datetime(end_time)

    # Filter based on time
    segment = df[(df['time'] >= start_dt) & (df['time'] <= end_dt)].copy()

    # Convert time to seconds relative to start time
    segment['time'] = (segment['time']-start_dt).dt.total_seconds()
    print(f"Filtered {(segment['time'])} rows based on time")
    # Assume 2nd column is intensity/counts
    if flt == 'flare_3':
        segment['intensity'] = segment['intensity_cdte1']+segment['intensity_cdte2']
        print('__________________________________________________')
        print(segment['intensity'],segment['intensity_cdte1'],segment['intensity_cdte2'])
        print('-----------/////\\\\\----')
    intensity_col = segment.columns[1]
    segment['error'] = np.sqrt((segment['intensity']).astype(int))

    # Reorder columns
    segment = segment[['time','intensity', 'error']]

    # Save to file
    segment.to_csv(output_file, index=False, header=False)
    print(f"Saved segment with seconds and Poisson errors to {output_file}")

