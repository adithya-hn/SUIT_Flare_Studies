import numpy as np
import glob
import matplotlib.pyplot as plt

files = sorted(glob.glob("csv_files/*.csv"))

all_times = []
all_ints = []

for f in files:
    data = np.genfromtxt(f, delimiter=',', names=True, dtype=None, encoding='utf-8')
    
    # Convert time to datetime64
    t = np.array(data['Peak_time'], dtype='datetime64[ms]')
    x = np.array(data['Peak_total_count'], dtype=float)
    
    all_times.append(t)
    all_ints.append(x)

plt.figure(figsize=(8,5))

for t, x in zip(all_times, all_ints):
    plt.plot(t, x, alpha=0.7)

plt.xlabel("Time")
plt.ylabel("Intensity")
plt.tight_layout()
plt.close()

all_c = np.concatenate(all_ints)
print(len(all_c))
plt.figure(figsize=(6,5))
plt.hist(all_c, bins=20)
plt.xlabel("Peak Total Count")
plt.ylabel("Frequency")
plt.tight_layout()
plt.show()


# plt.figure(figsize=(6,5))

# for i, c in enumerate(all_ints):
#     plt.hist(c, bins=40, alpha=0.5, label=f"File {i}")

# plt.legend()
# plt.xlabel("Peak Total Count")
# plt.ylabel("Frequency")
# plt.tight_layout()
# plt.show()