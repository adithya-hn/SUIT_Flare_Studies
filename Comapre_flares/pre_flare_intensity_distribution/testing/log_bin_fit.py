import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')

# data=data[(data>0.36) & (data<1.0)]
# data=data[data>1.0]
data=data[data>3.6]
data=np.sort(data)
def empirical_ccdf(E):
    n = len(E)
    return np.arange(n, 0, -1) / n
ccdf = empirical_ccdf(data)

# data = data[data > x_min]

# define bins
nbins = 18
bins = np.logspace(np.log10(data.min()), np.log10(data.max()), nbins)
counts, edges = np.histogram(data, bins=bins)
bin_widths = np.diff(edges)
bin_centers = np.sqrt(edges[:-1] * edges[1:])  # geometric mean

pdf = counts / (np.sum(counts) * bin_widths)
# mask = counts > 0
bin_centers = bin_centers#[mask]
pdf = pdf#[mask]
bin_widths = np.diff(edges)


dNdE = counts / bin_widths
dNdE = dNdE

mask = (bin_centers > 0) & (bin_centers < 3)

x = np.log10(bin_centers[mask])
y = np.log10(dNdE[mask])

coeffs = np.polyfit(x, y, 1)#, w=w)

slope = coeffs[0]
intercept = coeffs[1]
alpha = -slope
A = 10**intercept

print("alpha =", alpha)
plt.figure(figsize=(10,6))
plt.loglog(bin_centers, dNdE, drawstyle="steps-mid",)
xfit = np.linspace(min(bin_centers[mask]), max(bin_centers[mask]), 100)
yfit = A * xfit**(-alpha)
plt.loglog(xfit, yfit, '-',label=fr'Power law fit  $\alpha$= {alpha:.2f}'+'\n'+ f'number of nbins= {nbins}')
plt.xlabel('Peak Flux (erg x cm^-2)')
plt.ylabel('dN/dF')
plt.title('Differential Frequency Distribution')
plt.legend()
plt.savefig(f'histogram_{nbins}.png',dpi=300)
plt.show()

# plt.loglog(bin_centers, pdf, 'o')
# plt.xlabel('Peak Flux')
# plt.ylabel('PDF')
# plt.show()