
''''
Created on 17 oct 2025
@author: adithya-hn

A simple example to load and plot STIX spectrum and lightcurve using SunPy and SunKit-Spex

'''


import matplotlib.pyplot as plt
from parfive import Downloader
import astropy.units as u
from astropy.time import Time
from sunkit_spex.extern.stix import STIXLoader
from sunkit_spex.legacy.fitting.fitter import Fitter, load
from datetime import datetime , timedelta
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
#---------------Input parameters----------------

spec_file="/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/stix/stx_spectrum_2410315184.fits"
srm_file="/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/stix/stx_srm_2410315184.fits"
Start_t="2024-11-01T00:00:30"
End_t="2024-11-01T02:30:30"
case='7_Nov01'

#-----------------------------------------------

time_profile_size = (9, 6)
spec_plot_size = (10, 10)
joint_spec_plot_size = (25, 10)
tol = 1e-20
spec_font_size = 18
xlims, ylims = [3, 100], [1e-1, 1e6]

plt.rcParams["font.size"] = spec_font_size
plt.figure(layout="tight",figsize=(12,8))
stix_spec = STIXLoader(spectrum_file=spec_file, srm_file=srm_file)
stix_spec.update_event_times(start=Time(Start_t), end=Time(End_t))
plot_range=[datetime.fromisoformat(Start_t)- timedelta(minutes=10),datetime.fromisoformat(End_t)+ timedelta(minutes=10)]
ax=stix_spec.lightcurve(energy_ranges=[[10,30]]) #
# plt.xlim(plot_range[0],plot_range[1])
# plt.ylim(ylims)
# plt.savefig(f"stix_lc_{case}_{Start_t}_{End_t}.png", dpi=300)
# plt.close()
# lc = stix_spec.lightcurve(energy_ranges=[[10, 30]]) 
# 3. Extract Data for CSV
times_for_csv = []
flux_for_csv = []
# lc_data = stix_spec.lightcurve(energy_ranges=[[10, 30]]) 
# 1. Generate the plot and capture the Axes object
#    This draws the 10-30 keV line onto the figure.
# ax = stix_spec.lightcurve(energy_ranges=[[10, 30]])

# 2. Extract data directly from the plot lines
#    Since 'ax' contains the plot, we can grab the x (time) and y (flux) arrays.
try:
    lines = ax.get_lines()
    
    # Assuming the first line is the one we want (usually true for single energy range)
    x_data = lines[0].get_xdata()
    y_data = lines[0].get_ydata()
    print(len(x_data))
    timestamps = mdates.num2date(x_data)
    timestamps_iso = [t.strftime('%Y-%m-%dT%H:%M:%S') for t in timestamps]

    # Create raw DataFrame
    df_raw = pd.DataFrame({
        'Time_UT': timestamps_iso,
        'Flux_10_30keV': y_data
    })

    # --- CLEANING STEP ---
    # The plot has duplicates for timestamps (vertical steps).
    # We keep the 'last' entry for each timestamp to get the flux STARTING at that time.
    df_clean = df_raw.drop_duplicates(subset='Time_UT', keep='last')

    # Save Clean CSV
    csv_filename = f"stix_lc_10-30keV_{case}_{Start_t}_{End_t}.csv"
    df_clean.to_csv(csv_filename, index=False)
    print(f"✅ Saved Cleaned CSV to: {csv_filename}")
    print(f"   (Rows reduced from {len(df_raw)} to {len(df_clean)})")

except IndexError:
    print("❌ Error: No plot lines found.")
except Exception as e:
    print(f"❌ Error extracting data: {e}")

# 5. Finalize Plot Visuals
ax.set_xlim(plot_range[0], plot_range[1])
ax.set_ylim(ylims)
ax.set_title("STIX Light Curve 10-30 keV")
plt.ylabel("Counts / s / keV / cm²") # Ensure units match your data

plt.savefig(f"stix_lc_{case}_{Start_t}_{End_t}.png", dpi=300)
plt.close()

fig_temp, ax_temp = plt.subplots()
stix_spec.lightcurve(energy_ranges=[[10, 30]]) 

# 3. Extract & Clean Raw Data
try:
    lines = ax_temp.get_lines()
    x_data = lines[0].get_xdata() # Matplotlib dates (floats)
    y_data = lines[0].get_ydata() # Flux / Count Rate

    # Convert Matplotlib floats to Pandas Datetime objects
    dates = mdates.num2date(x_data)
    
    # Create DataFrame
    df = pd.DataFrame({'Flux': y_data}, index=pd.to_datetime(dates))
    
    # Clean duplicates (handling the "step" plot style)
    # We keep the 'last' value for each exact timestamp
    df = df[~df.index.duplicated(keep='last')]

    # ---------------------------------------------------------
    # 4. RESAMPLE TO 1-MINUTE BINS
    # ---------------------------------------------------------
    # '1min' = 1 minute bin size.
    # .mean() calculates the average Flux/Rate during that minute.
    # Use .sum() ONLY if your unit is raw counts and you want total counts.
    df_1min = df.resample('1min').mean()
    
    # Remove any empty bins created by resampling (optional)
    df_1min = df_1min.dropna()

    # 5. Save 1-Min CSV
    csv_filename = f"stix_lc_10-30keV_1min_{case}_{Start_t}_{End_t}.csv"
    df_1min.to_csv(csv_filename, index_label='Time_UT')
    print(f"✅ Saved 1-minute binned CSV to: {csv_filename}")

except Exception as e:
    print(f"❌ Error processing data: {e}")
    df_1min = None

plt.close(fig_temp) # Close the temporary plot

# ---------------------------------------------------------
# 6. Plot the 1-Minute Binned Data
# ---------------------------------------------------------
if df_1min is not None:
    plt.figure(layout="tight", figsize=(12, 8))
    
    # Plotting the 1-min binned data
    # drawstyle="steps-mid" aligns the bin center nicely
    plt.plot(df_1min.index, df_1min['Flux'], 
             drawstyle="steps-mid", 
             color='blue', 
             linewidth=2, 
             label='10-30 keV (1-min bin)')

    # Formatting
    plot_range = [
        datetime.fromisoformat(Start_t) - timedelta(minutes=10),
        datetime.fromisoformat(End_t) + timedelta(minutes=10)
    ]
    plt.xlim(plot_range[0], plot_range[1])
    # plt.ylim(ylims)
    plt.yscale('log')
    plt.ylabel("Flux (counts/s/keV/cm²)")
    plt.title(f"STIX 10-30 keV Light Curve (1-minute bins)")
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Save Figure
    plt.savefig(f"stix_lc_1min_{case}_{Start_t}_{End_t}.png", dpi=300)
    plt.show()
    print("✅ Saved 1-minute binned Plot.")