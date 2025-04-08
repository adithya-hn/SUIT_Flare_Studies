import numpy as np
import pandas as pd
from scipy.integrate import trapz
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()


# Load light curve data (CSV without headers, columns: time (ISOT), total count, QS mean count, box area)
def load_light_curve(file):
    return pd.read_csv(file, names=["time", "total_count", "qs_mean", "box_area"], delimiter=',')

# Convert ISOT time to timestamp
def isot_to_timestamp(isot_time):
    return Time(isot_time, format='isot').unix

# Convert timestamp to HH:MM format
def timestamp_to_hhmm(timestamp):
    return Time(timestamp, format='unix').datetime.strftime('%H:%M')

# Compute average intensity over a given time range
def average_intensity(df, start_time, end_time, column):
    subset = df[(df['time'] >= start_time) & (df['time'] <= end_time)]
    return subset[column].mean()

def peak_intensity(df, start_time, end_time, column):
    subset = df[(df['time'] >= start_time) & (df['time'] <= end_time)]
    return subset[column].max()

# Compute FWHM
def compute_fwhm(df, peak_intensity):
    half_max = peak_intensity / 2
    above_half_max = df[df['total_count'] >= half_max]
    if len(above_half_max) > 1:
        fwhm = above_half_max['time'].iloc[-1] - above_half_max['time'].iloc[0]
    else:
        fwhm = np.nan
    return fwhm

# Compute fluence
def compute_fluence(df, start_time, end_time, preflare_intensity):
    flare_data = df[(df['time'] >= start_time) & (df['time'] <= end_time)]
    excess_intensity = flare_data['total_count'] - preflare_intensity
    return trapz(excess_intensity, flare_data['time'])

# Convert timestamp to datetime object
def timestamp_to_datetime(timestamp):
    return Time(timestamp, format='unix').datetime


#----------------------------------------------------------------

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
    noaa_region=row['noaa_region']
    goes_class=row['goes_class']

#----------------------------------------------------------------
# User-defined parameters

    
    df['datetime'] = df['time'].apply(timestamp_to_datetime)
  # Constant box area

    # Compute intensities
    preflare_intensity  = average_intensity(df, start_time - 120, start_time, 'total_count')
    flare_peak_intensity= peak_intensity(df, start_time, peak_time+60, 'total_count')
    quiet_sun_intensity = average_intensity(df, start_time, end_time, 'qs_mean')
    preflare_level      = average_intensity(df, start_time - 3200, start_time-1200, 'total_count')

    # Compute properties
    intensity_ratio = flare_peak_intensity / preflare_intensity
    flare_contrast = (flare_peak_intensity - preflare_intensity) / preflare_intensity
    qs_peak_ratio=   (flare_peak_intensity/box_area) / quiet_sun_intensity
    enhancement_qs = ((flare_peak_intensity/box_area) - quiet_sun_intensity) / quiet_sun_intensity
    fluence = compute_fluence(df, start_time, end_time, preflare_intensity)
    fwhm = compute_fwhm(df, flare_peak_intensity)

    # Save results
    results = {
        "Flare_ID": f"{Time(start_time, format='unix').isot.replace(':', '').replace('-', '').replace('T', '_')}",
        "Filter": filter_name,
        "NOAA_Region": row["noaa_region"],
        "GOES_Class": row["goes_class"],
        "Start_Time_UTC": row["start_time"],
        "Peak_Time_UTC": row["peak_time"],
        "End_Time_UTC": row["end_time"],
        "Preflare_Intensity": preflare_intensity,
        "Flare_Peak_Intensity": flare_peak_intensity,
        "Quiet_Sun_Intensity": quiet_sun_intensity,
        "Intensity_Ratio": intensity_ratio,
        "Peak_Ratio": qs_peak_ratio,
        "Flare_Contrast": flare_contrast,
        "Enhancement_QS": enhancement_qs,
        "Fluence": fluence,
        "FWHM": fwhm,
        "Box_Area": box_area,
        "Preflare_level": preflare_level/(box_area*quiet_sun_intensity)
    }
    all_results.append(results)
    # Convert to DataFrame and save
    report_df = pd.DataFrame([results])
  
    # Generate plot
    plt.figure(figsize=(10, 5))
    plt.plot(df['datetime'], df['total_count'], label='Flare Light Curve', color='r')
    plt.axvline(timestamp_to_datetime(start_time), color='g', linestyle='--', label=f'Start Time')
    plt.axvline(timestamp_to_datetime(peak_time), color='b', linestyle='--', label=f'Peak Time')
    plt.axvline(timestamp_to_datetime(end_time), color='k', linestyle='--', label=f'End Time')
    plt.axhline(preflare_intensity, color='purple', linestyle='--', label='Preflare Intensity')
    plt.axhline(flare_peak_intensity, color='orange', linestyle='--', label='Flare Peak Intensity')
    plt.axhline(preflare_level, color='aqua', linestyle='--', label='Pre flare level')

    # Annotate properties
    props_text = (
                  f"Intensity Ratio: {intensity_ratio:.2f}\n"
                  f"Peak_QS Ratio: {qs_peak_ratio:.2f}\n"
                  f"Flare Contrast: {flare_contrast:.2f}\n"
                  f"Enhancement (QS): {enhancement_qs:.2f}\n"
                  f"Fluence: {fluence:.2f}\n"
                  f"Box Area: {box_area}\n"
                  f"Pre flare level QS ratio: {preflare_level/(box_area*quiet_sun_intensity)}")

    
    if filter_name=='NB08':
        plt.text(df['datetime'].min(), df['total_count'].max() * 0.99, props_text, fontsize=10,bbox=dict(facecolor='white', alpha=0.7))
    else:
        plt.text(df['datetime'].min(), df['total_count'].max() * 0.9, props_text, fontsize=10,
             bbox=dict(facecolor='white', alpha=0.7))

    plt.xlabel('Time (HH:MM)')
    plt.ylabel('Total Count')
    plt.title(f'Flare Light Curve - NOAA {noaa_region} - {goes_class}')
    #plt.legend()
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    #plt.legend()
    plt.xticks(rotation=45)
    plot_filename = f"flare_plot_{file_name.split('.')[0]}.png"
    plt.savefig(plot_filename)
    plt.close()

# Save combined report
report_df = pd.DataFrame(all_results)
report_df.to_csv("combined_flare_report.csv", index=False)
print("Combined flare report saved.")
