import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')
x_min = 0.1
data=np.sort(data)
def empirical_ccdf(E):
    n = len(E)
    return np.arange(n, 0, -1) / n
ccdf = empirical_ccdf(data)

# data = data[data > x_min]

# define bins
nbins = 10
bins = np.logspace(np.log10(data.min()), np.log10(data.max()), nbins)
counts, edges = np.histogram(data, bins=bins)
bin_widths = np.diff(edges)
bin_centers = np.sqrt(edges[:-1] * edges[1:])  # geometric mean

pdf = counts / (np.sum(counts) * bin_widths)
mask = counts > 0
bin_centers = bin_centers[mask]
pdf = pdf[mask]
bin_widths = np.diff(edges)
dNdE = counts / bin_widths
dNdE = dNdE[mask]

mask = (bin_centers > 0.1) & (bin_centers < 2)

x = np.log10(bin_centers[mask])
y = np.log10(dNdE[mask])

# convert errors to log-space
# yerr = errors[mask] / (dNdE[mask] * np.log(10))

# weights
# w = 1 / yerr**2

# weighted linear fit
coeffs = np.polyfit(x, y, 1)#, w=w)

slope = coeffs[0]
intercept = coeffs[1]

alpha = -slope
A = 10**intercept

print("alpha =", alpha)



plt.loglog(bin_centers, dNdE, 'o')
plt.xlabel('Energy')
plt.ylabel('dN/dE')
plt.title('Differential Frequency Distribution')
plt.show()

plt.loglog(bin_centers, pdf, 'o')
plt.xlabel('Energy')
plt.ylabel('PDF')
plt.show()

logx = np.log10(bin_centers)
logy = np.log10(pdf)

# choose fitting range (IMPORTANT)
mask = (bin_centers > 0.1) & (bin_centers < 2.2)

coeffs = np.polyfit(logx[mask], logy[mask], 1)
slope = coeffs[0]

alpha_est = -slope
print(alpha_est)

xmin = 0.3
data_ = data[data > xmin]

alpha_mle = 1 + len(data_) / np.sum(np.log(data_ / xmin))
print(alpha_mle)

# def neg_log_likelihood(params):
#     alpha, x0 = params
    
#     if alpha <= 1 or x0 < 0:
#         return np.inf
    
#     xmin = x_min
#     xmax = data.max()
    
#     # normalization
#     C = ((xmax + x0)**(1 - alpha) - (xmin + x0)**(1 - alpha)) / (1 - alpha)
        
#     logL = -alpha * np.sum(np.log(data + x0)) - len(data) * np.log(C)
    
#     return -logL
# res = minimize(neg_log_likelihood, x0=[2.0, 0.05])
# alpha_opt, x0_opt = res.x
# print(alpha_opt, x0_opt)
alpha_final=alpha_mle
Emin_opt=xmin
Efit = data[data >= Emin_opt]
N0 = len(Efit) 
model_ccdf = .6 * (Efit / Emin_opt) ** (-(alpha_final - 1))


plt.figure(figsize=(10,6))
plt.loglog(data, ccdf, '+',color='k', label='SUIT transeint events')
plt.loglog(Efit, model_ccdf, '-',
           linewidth=2,
           label=f'Power-law fit (alpha = {alpha_final:.2f})')

plt.show()