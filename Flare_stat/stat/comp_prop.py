import numpy as np
import pandas as pd
from scipy.integrate import trapz
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os

# Load flare cases
flare_cases = pd.read_csv("Flare_prop_time.csv")

# Initialize list to store results
all_results = []

# Process each flare case
for _, row in flare_cases.iterrows():
    file_name = row["filename"]  # First column contains the file name
    flare_class = file_name.split("_")[1]  # Extract flare class from filename
    filter_name = file_name[:4]  # First four characters are the filter
    start_time = Time(row["start_time"], format='isot').unix
    peak_time = Time(row["peak_time"], format='isot').unix
    end_time = Time(row["end_time"], format='isot').unix
    
    # Load light curve data
    df = pd.read_csv(file_name, names=["time", "total_count", "qs_mean", "box_area"], delimiter=',')
    df['time'] = df['time'].apply(lambda t: Time(t, format='isot').unix)
    df['datetime'] = df['time'].apply(lambda t: Time(t, format='unix').datetime)
    box_area = df['box_area'].iloc[0]
    
    # Compute properties
    preflare_intensity = df[(df['time'] >= start_time - 120) & (df['time'] < start_time)]['total_count'].mean()
    flare_peak_intensity = df[(df['time'] >= peak_time - 120) & (df['time'] < peak_time)]['total_count'].mean()
    quiet_sun_intensity = df[(df['time'] >= start_time) & (df['time'] <= end_time)]['qs_mean'].mean()
    fluence = trapz(df[(df['time'] >= start_time) & (df['time'] <= end_time)]['total_count'] - preflare_intensity, df[(df['time'] >= start_time) & (df['time'] <= end_time)]['time'])
    
    # Store results
    results = {
        "File": file_name,
        "Filter": filter_name,
        "Flare_Class": flare_class,
        "Start_Time": row["start_time"],
        "Peak_Time": row["peak_time"],
        "End_Time": row["end_time"],
        "Preflare_Intensity": preflare_intensity,
        "Flare_Peak_Intensity": flare_peak_intensity,
        "Quiet_Sun_Intensity": quiet_sun_intensity,
        "Fluence": fluence,
    }
    all_results.append(results)
    
    # Plot
    plt.figure(figsize=(10, 5))
    plt.plot(df['datetime'], df['total_count'], label='Flare Light Curve', color='r')
    plt.axvline(Time(start_time, format='unix').datetime, color='g', linestyle='--', label='Start Time')
    plt.axvline(Time(peak_time, format='unix').datetime, color='b', linestyle='--', label='Peak Time')
    plt.axvline(Time(end_time, format='unix').datetime, color='k', linestyle='--', label='End Time')
    plt.xlabel('Time (HH:MM)')
    plt.ylabel('Total Count')
    plt.title(f'Flare Light Curve - {filter_name} {flare_class}')
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.xticks(rotation=45)
    plt.legend()
    plot_filename = f"flare_plot_{file_name.split('.')[0]}.png"
    plt.savefig(plot_filename)
    plt.close()
    
# Save combined report
report_df = pd.DataFrame(all_results)
report_df.to_csv("combined_flare_report.csv", index=False)
print("Combined flare report saved.")
