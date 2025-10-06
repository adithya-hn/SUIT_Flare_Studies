import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import numpy as np
from datetime import datetime, timedelta

# --- Load your AIA light curves ---
aia_df = pd.read_csv("aia131_region_lc.csv", parse_dates=["time"])
# --- Bin AIA to 30s cadence (matching HELIOS) ---
aia_df["time"] = pd.to_datetime(aia_df["time"])
aia_df = aia_df.set_index("time")
aia_filtered = aia_df[aia_df["exposure"] >= 0.1]

# Drop the exposure column if you don't want it in the output
#aia_filtered = aia_filtered.drop(columns=["exposure"])

# average values in 30s bins
aia_binned = aia_filtered.resample("60s").mean().reset_index()
# --- Normalize by exposure (elementwise); avoid divide-by-zero
aia_norm = aia_binned.copy()
region_cols = [c for c in aia_binned.columns if c.startswith("region_")]
aia_norm[region_cols] = aia_norm[region_cols].div(aia_norm["exposure"].replace(0, np.nan), axis=0)


#print(aia_binned)

# --- Load your HELIOS count rate file ---
# Suppose it has columns: time, helios_rate
Helios=(np.load("cdte_data_flare_5.npy", allow_pickle=True)).transpose()
datetime_objects = pd.to_datetime(Helios[0])
helios_time=helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]

helios_rate=np.array((Helios[1]+Helios[2]),dtype=float)
helios_df = pd.DataFrame({"time": helios_time, "helios_rate": helios_rate},index=helios_time)
helios_df["time"] = pd.to_datetime(helios_df["time"])
helios_df = helios_df.set_index("time")

print(helios_df.head())
#HELIOS_START = pd.Timestamp('2024-06-01 07:00:00') 


# --- Convert HELIOS to 1-minute bins: sum two 30s bins -> 1 minute counts
helios_df_binned = helios_df.resample("1min").mean()


print('Helios',helios_df_binned.iloc[0],'AIA',aia_binned.iloc[0])

#print(datetime_objects)

# --- Align times by merging ---
df = pd.merge_asof(aia_binned.sort_values("time"),
                   helios_df_binned.sort_values("time"),
                   on="time", direction="nearest")




# --- Plot all regions with HELIOS ---


correlations = {}

for col in df.columns:
    if col.startswith("region_"):
        valid = df[[col, "helios_rate"]].dropna()
        plt.figure(figsize=(12, 6))
        ax=plt.subplot(111)
        ax2=ax.twinx()
        ax.plot(df["time"], df["helios_rate"], label="HELIOS", color="black", lw=2)
        ax.set_yscale('log')
        if len(valid) > 2:
            r, _ = pearsonr(valid[col]/df["exposure"], valid["helios_rate"])
            correlations[col] = r
        ax2.plot(df["time"], df[col]/df["exposure"], label=f'{col} with r= {r:.02f}')  # normalized by exposure
        ax2.set_yscale('log')
        plt.legend()
        plt.xlabel("Time")
        ax2.set_ylabel("AIA reg Counts")
        ax.set_ylabel("AIA reg Counts")
        plt.title("HELIOS vs AIA Region Light Curves")
        plt.tight_layout()
        plt.savefig(f'aia_reg{col}_helios.png',dpi=300)
        plt.close()

'''# --- Compute correlation values ---

for col in df.columns:
    if col.startswith("region_"):
        valid = df[[col, "helios_rate"]].dropna()
        if len(valid) > 2:
            #print(valid[col]/df["exposure"],valid["helios_rate"])
            r, p = pearsonr(valid[col]/df["exposure"], valid["helios_rate"])
            correlations[col] = r
'''
print("Correlation values with HELIOS:")
for region, r in correlations.items():
    print(f"{region}: {r:.3f}")


plt.figure(figsize=(12, 6))
ax1=plt.subplot(111)
ax1.plot(df["time"], df["helios_rate"], label="HELIOS", color="black", lw=2)
ax2=ax1.twinx()
for col in df.columns:
    if col.startswith("region_"):
        ax2.plot(df["time"], df[col], label=f'{col}: {correlations[col]:.02f}')
    

plt.legend()
plt.xlabel("Time")
ax1.set_ylabel("HEL1OS Counts")
ax2.set_ylabel("aia reg counts")
ax1.set_yscale("log")
ax2.set_yscale("log")
plt.title("HELIOS vs AIA Region Light Curves (1min cadence)")
plt.tight_layout()
plt.savefig('AIA-HEL1OS-lightcurve.png',dpi=300)
plt.show()