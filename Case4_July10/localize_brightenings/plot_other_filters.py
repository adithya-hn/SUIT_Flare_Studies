import pandas as pd
import matplotlib.pyplot as plt
import glob

# Path to your CSV files (change accordingly)
files = glob.glob("Diff*.csv")   # e.g. "./data/*.csv"

plt.figure(figsize=(10,6))

f='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/localize_brightenings/NB05_sel_contour_data.csv'
df = pd.read_csv(f)
df['image datetime'] = pd.to_datetime(df['image datetime'])
# Detect gaps > 2 minutes and break by inserting NaN
dt = df['image datetime'].diff().dt.total_seconds()
df.loc[dt > 120, 'intensity_sum'] = None  # gap threshold = 120 sec

# Plot intensity with error bar (optional)
plt.plot(df['image datetime'], df['intensity_sum'], alpha=0.6, label=f.split('/')[-1])
#plt.plot(df['Date'], df['diff_count'], label=f.split('/')[-1])   # simple line plot
     

plt.xlabel("Time")
plt.ylabel("Intensity")
plt.title("Light Curves")
#plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('plot_nb5_selected_contour_lc.png',dpi=300)
plt.show()