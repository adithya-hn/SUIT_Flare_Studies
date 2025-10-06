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

# average values in 30s bins
aia_binned = aia_df.resample("60s").mean().reset_index()
#print(aia_binned)

# --- Load your HELIOS count rate file ---
# Suppose it has columns: time, helios_rate
Helios=(np.load("cdte_data_flare_1.npy", allow_pickle=True)).transpose()
datetime_objects = pd.to_datetime(Helios[0])
helios_time=helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]

helios_rate=np.array((Helios[1]+Helios[2]),dtype=float)
helios_df = pd.DataFrame({"time": helios_time, "helios_rate": helios_rate})

# Convert time to datetime
#aia_df["time"]    = pd.to_datetime(aia_df["time"])
#aia_df = aia_df.set_index("time").sort_index()

helios_df["time"] = pd.to_datetime(helios_df["time"])
helios_df = helios_df.set_index("time")
helios_df_binned = helios_df.resample("60s").mean().reset_index()
print(helios_df_binned)

#print(datetime_objects)

# --- Align times by merging ---
df = pd.merge_asof(aia_binned.sort_values("time"),
                   helios_df_binned.sort_values("time"),
                   on="time", direction="nearest")

# Now df has: time, exposure, region_1, region_2,..., helios_rate
# print(df.head())

# --- Plot all regions with HELIOS ---
plt.figure(figsize=(12, 6))
plt.plot(df["time"], df["helios_rate"], label="HELIOS", color="black", lw=2)
plt.yscale('log')

correlations = {}

for col in df.columns:
    if col.startswith("region_"):
        valid = df[[col, "helios_rate"]].dropna()
        if len(valid) > 2:
            r, _ = pearsonr(valid[col]/df["exposure"], valid["helios_rate"])
            correlations[col] = r
        plt.plot(df["time"], df[col]/df["exposure"], label=f'{col} with r= {r:.02f}')  # normalized by exposure

plt.legend()
plt.xlabel("Time")
plt.ylabel("Counts / Flux")
plt.title("HELIOS vs AIA Region Light Curves")
plt.tight_layout()
plt.show()

# --- Compute correlation values ---

for col in df.columns:
    if col.startswith("region_"):
        valid = df[[col, "helios_rate"]].dropna()
        if len(valid) > 2:
            #print(valid[col]/df["exposure"],valid["helios_rate"])
            r, p = pearsonr(valid[col]/df["exposure"], valid["helios_rate"])
            correlations[col] = r

print("Correlation values with HELIOS:")
for region, r in correlations.items():
    print(f"{region}: {r:.3f}")
