import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import numpy as np
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()



cols_to_plot = [ "region_2", "region_3","region_5"]  # change names as needed


# --- Load your AIA light curves ---
aia_df = pd.read_csv("aia131_region_lc.csv", parse_dates=["time"])
# --- Bin AIA to 30s cadence (matching HELIOS) ---
aia_df["time"] = pd.to_datetime(aia_df["time"])
aia_df = aia_df.set_index("time")
aia_filtered = aia_df[aia_df["exposure"] >= 0.1]
aia_norm = aia_filtered.copy()
region_cols = [c for c in aia_filtered.columns if c.startswith("region_")]
aia_norm[region_cols] = aia_norm[region_cols].div(aia_norm["exposure"].replace(0, np.nan), axis=0)
aia_binned = aia_norm.resample("60s").mean().reset_index()

Helios=(np.load("cdte_data_flare_7.npy", allow_pickle=True)).transpose()
datetime_objects = pd.to_datetime(Helios[0])
helios_time=helio_time_array=[datetime.strptime(str(ts)[:26], "%Y-%m-%d %H:%M:%S.%f") for ts in datetime_objects]

helios_rate=np.array((Helios[1]+Helios[2]),dtype=float)
helios_df = pd.DataFrame({"time": helios_time, "helios_rate": helios_rate},index=helios_time)
helios_df["time"] = pd.to_datetime(helios_df["time"])
helios_df = helios_df.set_index("time")

helios_df_binned = helios_df#.resample("1min").mean()

combined_df = pd.merge_asof(aia_binned.sort_values("time"),
                   helios_df_binned.sort_values("time"),
                   on="time", direction="nearest")
# Plot
#combined_df[cols_to_plot]#.plot(figsize=(10, 5))


plt.figure(figsize=(12, 6))
ax1=plt.subplot(111)
ax1.plot(combined_df["time"], combined_df["helios_rate"], label="HELIOS", color="black", lw=2)
ax2=ax1.twinx()

for col in cols_to_plot:
    valid = combined_df[[col, "helios_rate"]].dropna()
    r, _ = pearsonr(valid[col], valid["helios_rate"])
    ax2.plot(combined_df["time"], combined_df[col], label=f'{col}: {r:.02f}')


ax2.legend()
plt.xlabel("Time")
ax1.set_ylabel("HEL1OS Counts")
ax2.set_ylabel("aia reg counts")
ax1.set_yscale("log")
ax2.set_yscale("log")
plt.title("HELIOS vs AIA Region Light Curves (1min cadence)")

time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)

plt.tight_layout()
plt.savefig('AIA-HEL1OS-lc.png',dpi=300)
plt.show()

# plt.xlabel("Time")
# plt.ylabel("Counts / Intensity")
# plt.title("Selected Light Curves")
# plt.legend()
# plt.tight_layout()
# plt.show()
