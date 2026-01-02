import pandas as pd
import matplotlib.pyplot as plt
import glob
import numpy as np

# Path to your CSV files (change accordingly)
files = glob.glob("csv_files/*.csv")   # e.g. "./data/*.csv"
aia131_fl=pd.read_csv('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/localize_brightenings/other_filters/aia131_region_lc.csv')

aia_dt= pd.to_datetime(aia131_fl['time'])
reg=aia131_fl['region_1']/aia131_fl['exposure']
fig, ax = plt.subplots(figsize=(10,6))
ax1=ax.twinx()
for f in files:
    df = pd.read_csv(f)
    
    # If datetime column needs parsing
    df['Time'] = pd.to_datetime(df['Time'])
    # Detect gaps > 2 minutes and break by inserting NaN
    dt = df['Time'].diff().dt.total_seconds()
    df.loc[dt > 120, 'AR_total'] = None  # gap threshold = 120 sec
    
    # Plot intensity with error bar (optional)
    #plt.errorbar(df['Time'], df['diff_count'], yerr=df['Diff_error'], alpha=0.6, label=f.split('/')[-1])
    ax.plot(df['Time'], df['AR_total']/np.max(df['AR_total']), label=f.split('/')[-1])   # simple line plot
     
ax1.plot(aia_dt,reg,label='AIA 131',color='k')
plt.xlabel("Time")
plt.ylabel("Intensity")
plt.title("Light Curves")
#plt.yscale('log')
ax.legend()
ax1.legend(loc='upper right')
plt.tight_layout()
plt.savefig('plot_all_diff_img_lc.png',dpi=300)
plt.show()
