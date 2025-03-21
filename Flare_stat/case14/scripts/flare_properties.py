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
# User-defined parameters
#flare_light_curve_file = 'NB08_X1.4_Light_curve_data.csv'
flare_light_curve_files = ['NB03_c14_lc_data.csv', 'NB04_c14_lc_data.csv', 'NB08_c14_lc_data.csv']

start_time_isot = '2024-11-10T11:51:00'
peak_time_isot = '2024-11-10T12:06:00'
end_time_isot = '2024-11-10T12:14:00'
noaa_region = 13889
goes_class = 'M9.4'

#----------------------------------------------------------------

# Convert ISOT times to timestamps
start_time = isot_to_timestamp(start_time_isot)
peak_time = isot_to_timestamp(peak_time_isot)
end_time = isot_to_timestamp(end_time_isot)
filters_results=[]
for file in flare_light_curve_files:
    # Load light curve
    df = load_light_curve(file)
    df['time'] = df['time'].apply(isot_to_timestamp)
    df['datetime'] = df['time'].apply(timestamp_to_datetime)
    box_area = df['box_area'].iloc[0]  # Constant box area

    # Compute intensities
    preflare_intensity  = average_intensity(df, start_time - 120, start_time, 'total_count')
    flare_peak_intensity= average_intensity(df, peak_time - 120, peak_time, 'total_count')
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
        "NOAA_Region": noaa_region,
        "GOES_Class": goes_class,
        "Start_Time_UTC": start_time_isot,
        "Peak_Time_UTC": peak_time_isot,
        "End_Time_UTC": end_time_isot,
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

    # Convert to DataFrame and save
    report_df = pd.DataFrame([results])
    report_filename = f"flare_report_{file.split('.')[0]}.csv"
    report_df.to_csv(report_filename, index=False)
    print(f"Flare analysis complete. Report saved as {report_filename}.")

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

    
    if file=='NB08_c14_lc_data.csv':
        plt.text(df['datetime'].min(), df['total_count'].max() * 0.99, props_text, fontsize=10,bbox=dict(facecolor='white', alpha=0.7))
    else:
        plt.text(df['datetime'].min(), df['total_count'].max() * 0.9, props_text, fontsize=10,
             bbox=dict(facecolor='white', alpha=0.7))

    plt.xlabel('Time (HH:MM)')
    plt.ylabel('Total Count')
    plt.title(f'Flare Light Curve - NOAA {noaa_region} - {goes_class} ({file})')
    #plt.legend()
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.xticks(rotation=45)
    plot_filename = f"flare_report_{file.split('.')[0]}.png"
    plt.savefig(plot_filename)
    plt.close()
    print(f"Plot saved as {plot_filename}.")