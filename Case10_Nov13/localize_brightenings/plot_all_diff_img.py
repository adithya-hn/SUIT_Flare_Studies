import pandas as pd
import matplotlib.pyplot as plt
import glob

# Path to your CSV files (change accordingly)
files = glob.glob("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case10_Nov13/localize_brightenings/Diff*.csv")   # e.g. "./data/*.csv"

plt.figure(figsize=(10,6))

for f in files:
    df = pd.read_csv(f)
    
    # If datetime column needs parsing
    df['Date'] = pd.to_datetime(df['Date'])
    # Detect gaps > 2 minutes and break by inserting NaN
    dt = df['Date'].diff().dt.total_seconds()
    df.loc[dt > 120, 'diff_count'] = None  # gap threshold = 120 sec
    
    # Plot intensity with error bar (optional)
    plt.errorbar(df['Date'], df['diff_count'], yerr=df['Diff_error'], alpha=0.6, label=f.split('/')[-1])
    #plt.plot(df['Date'], df['diff_count'], label=f.split('/')[-1])   # simple line plot
     

plt.xlabel("Time")
plt.ylabel("Intensity")
plt.title("Light Curves")
#plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.show()
