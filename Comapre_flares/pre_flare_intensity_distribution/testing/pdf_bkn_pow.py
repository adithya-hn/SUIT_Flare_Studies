"""
@Author      : Adithya H N
@Created On  : 2026-04-12
@Last Updated: 2026-04-15
@Project     : Pre-flare study
@Version     : 1.0

@Description
-----------
Brief description: Trying power law fit on log bin data, power law varies based on the bin wodth.


"""



import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')

data=np.sort(data[data>0.33])

hist, edges = np.histogram(data, bins=np.logspace(np.log10(min(data)), np.log10(max(data)), 18),  density=True)

x = np.sqrt(edges[:-1] * edges[1:])   # bin centers (geometric mean)

# Remove zero bins
mask = hist > 0
xfit = x[mask]
yfit = hist[mask]

# Broken power-law model
def broken_powerlaw(x, A, xb, a1, a2):
    return np.where(
        x < xb,
        A * (x / xb)**(-a1),
        A * (x / xb)**(-a2)
    )

# Initial guesses
p0 = [1.0, 1.0, 1.0, 3.0]

# Fit
popt, pcov = curve_fit(
    broken_powerlaw,
    xfit,
    yfit,
    p0=p0,
    bounds=([0, min(xfit), 0, 0],
            [np.inf, max(xfit), 10, 10])
)

A_fit, xb_fit, a1_fit, a2_fit = popt

print("Break =", xb_fit)
print("Slope1 =", a1_fit)
print("Slope2 =", a2_fit)

# Plot
xx = np.logspace(np.log10(min(xfit)), np.log10(max(xfit)), 20)
bins = np.logspace(np.log10(min(data)), np.log10(max(data)), 20)

plt.hist(data, bins=bins, density=True, histtype='step')

# plt.loglog(xfit, yfit, 'o')
plt.loglog(xx, broken_powerlaw(xx, *popt), 'r-')
plt.xlabel("x")
plt.ylabel("PDF")
plt.show()