    # import pickle

    # with open("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case2_June02/data/helios/cdte_corr_spgm_2024_06_02_07_41_to_2024_06_02_09_01.pkl", "rb") as f:
    #     data = pickle.load(f)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_pickle("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case2_June02/data/helios/cdte_corr_spgm_2024_06_02_07_41_to_2024_06_02_09_01.pkl")
print(type(df))
print(df.keys())

df["time"] = pd.to_datetime(df["time"])

def band_counts(energy, counts, emin=22, emax=30):
    energy = np.asarray(energy)
    counts = np.asarray(counts)
    mask = (energy >= emin) & (energy <= emax)
    return np.nansum(counts[mask])
df["lc_22_30_keV"] = df.apply(
    lambda row: band_counts(
        row["cdte1_ene_vec"],
        row["cdte1_spgm_counts"],
        22, 30
    ),
    axis=1
)


plt.figure(figsize=(8,4))
plt.plot(df["time"], df["lc_22_30_keV"], drawstyle="steps-mid")
plt.xlabel("Time")
plt.ylabel("Counts (22–30 keV)")
plt.title("HELIOS CdTe1 Light Curve (22–30 keV)")
plt.tight_layout()
plt.yscale('log')
plt.show()
