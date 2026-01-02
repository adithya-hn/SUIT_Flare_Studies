import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


df=pd.read_csv('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/localize_brightenings/contour_stats.csv')


Bins=np.arange(df["area_px"].min(),df['area_px'].max(),50)
print(len(Bins))
plt.figure(figsize=(6,4))
plt.hist(df["area_px"], bins=Bins)
plt.xlabel("Contour area (pixels)")
plt.ylabel("Count")
plt.tight_layout()
plt.show()